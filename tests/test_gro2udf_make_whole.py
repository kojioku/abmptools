# -*- coding: utf-8 -*-
"""
gro2udf ``--make-whole``: 各フレームの分子を、.top の結合をたどってつなぎ直す。

粗視化の本番 MD では、mdrun の .gro / .xtc が**全フレームで**分子が箱の境界を
またいで割れていることがある。``-pbc nojump`` (``--trajectory`` の既定) は
最初のフレームの状態を引き継ぐだけで割れを直さず、出力の間隔が粗い (10 ns ごと
など) と小分子がフレーム間に半箱動いて、本当の移動を折り返しと取り違える。
``-pbc mol`` は .tpr が要り、.tpr は書いた GROMACS 以上の版でないと読めない。

ここでは [ bonds ] と [ constraints ] をたどって、各フレームを独立に直す。
"""
import numpy as np
import pytest

from abmptools.gro2udf.make_whole import make_frames_whole, unwrap_steps
from abmptools.gro2udf.top_model import GROFrameData
from abmptools.gro2udf.top_parser import TopParser

#: 結合でつながる鎖 (CHN, 4 粒子) と、拘束だけでつながる三角形 (TRI, 3 粒子)、
#: 1 粒子の溶媒 (SOL)。comb-rule 2 の普通の型。
TOP = """\
[ defaults ]
1  2  no  1.0  1.0

[ atomtypes ]
  PA   45.0 0.000 A 0.43 2.0
  W    72.0 0.000 A 0.47 3.0

[ moleculetype ]
  CHN  1

[ atoms ]
  1 PA 1 CHN C1 1 0.0 45.0
  2 PA 1 CHN C2 2 0.0 45.0
  3 PA 1 CHN C3 3 0.0 45.0
  4 PA 1 CHN C4 4 0.0 45.0

[ bonds ]
  1 2 1 0.40 7000
  2 3 1 0.40 7000
  3 4 1 0.40 7000

[ moleculetype ]
  SOL  1

[ atoms ]
  1 W 1 SOL W 1 0.0 72.0

[ moleculetype ]
  TRI  1

[ atoms ]
  1 PA 1 TRI T1 1 0.0 45.0
  2 PA 1 TRI T2 2 0.0 45.0
  3 PA 1 TRI T3 3 0.0 45.0

[ constraints ]
  1 2 1 0.30
  1 3 1 0.30
  2 3 1 0.30

[ system ]
  test

[ molecules ]
  CHN  1
  SOL  2
  TRI  1
"""
BOX = 3.0   # nm
#: 原子の並び: CHN 0-3 | SOL 4,5 | TRI 6-8   (計 9)


def _whole_coords(shift=0.0):
    """分子がつながった座標 (箱の中央付近)。"""
    chn = [[1.0 + 0.4 * k + shift, 1.5, 1.5] for k in range(4)]     # 1.0 .. 2.2
    sol = [[0.5, 0.5, 0.5], [2.5, 2.5, 2.5]]
    tri = [[1.5, 0.2, 1.5], [1.8, 0.2, 1.5], [1.65, 0.46, 1.5]]
    return chn + sol + tri


def _wrapped(coords):
    """各原子を箱 [0, BOX) に折り返す (mdrun の出力と同じ割れ方)。"""
    return [[c % BOX for c in xyz] for xyz in coords]


def _raw(tmp_path):
    p = tmp_path / "w.top"
    p.write_text(TOP)
    return TopParser().parse(str(p))


def _bond_lengths(coords, pairs):
    x = np.asarray(coords)
    return [float(np.linalg.norm(x[i] - x[j])) for i, j in pairs]


CHN_BONDS = [(0, 1), (1, 2), (2, 3)]
TRI_LINKS = [(6, 7), (6, 8), (7, 8)]


class TestWalk:
    def test_constraints_are_walked_even_without_constraints_as_bonds(self, tmp_path):
        """拘束だけの分子 (TRI) も、--constraints-as-bonds 無しでつなぎ直す。"""
        raw = _raw(tmp_path)
        assert raw.bondlist[2] == []                            # 結合としては 0 本
        assert raw.constraint_pairs["TRI"] == [(1, 2), (1, 3), (2, 3)]
        children = sorted(int(c) for _, cs in unwrap_steps(raw) for c in cs)
        assert children == [1, 2, 3, 7, 8]                      # CHN 1-3, TRI 7-8

    def test_split_molecules_are_put_back(self, tmp_path):
        # CHN を右端にずらして境界をまたがせ、TRI を下端にまたがせる
        coords = _whole_coords(shift=1.5)                       # CHN 2.5 .. 3.7
        coords[6][1] = coords[7][1] = -0.1                      # TRI がy=0をまたぐ
        coords[8][1] = 0.16
        split = _wrapped(coords)
        assert max(_bond_lengths(split, CHN_BONDS)) > 2.0       # 割れている
        assert max(_bond_lengths(split, TRI_LINKS)) > 2.0

        f = GROFrameData(step=0, time=0.0, coord_list=split, cell=[BOX] * 3)
        out = make_frames_whole([f], unwrap_steps(_raw(tmp_path)))[0]
        assert _bond_lengths(out.coord_list, CHN_BONDS) == pytest.approx([0.4] * 3)
        assert max(_bond_lengths(out.coord_list, TRI_LINKS)) < 0.31

    def test_only_whole_box_shifts_are_applied(self, tmp_path):
        """座標は箱の整数倍しか動かさない (値そのものは丸めない)。"""
        split = _wrapped(_whole_coords(shift=1.5))
        f = GROFrameData(step=0, time=0.0, coord_list=split, cell=[BOX] * 3)
        out = make_frames_whole([f], unwrap_steps(_raw(tmp_path)))[0]
        n = (np.asarray(out.coord_list) - np.asarray(split)) / BOX
        assert np.allclose(n, np.round(n), atol=1e-12)

    def test_each_frame_is_fixed_on_its_own(self, tmp_path):
        """前のフレームに頼らない: 箱の半分以上動いたフレームが続いても直る。"""
        steps = unwrap_steps(_raw(tmp_path))
        frames = [GROFrameData(step=k, time=10000.0 * k,
                               coord_list=_wrapped(_whole_coords(shift=s)),
                               cell=[BOX] * 3)
                  for k, s in enumerate([0.0, 1.6, 3.1, 0.7])]  # 1 フレームで半箱以上
        for out in make_frames_whole(frames, steps):
            assert _bond_lengths(out.coord_list, CHN_BONDS) == pytest.approx([0.4] * 3)

    def test_a_frame_without_a_box_is_refused(self, tmp_path):
        f = GROFrameData(step=0, time=0.0, coord_list=_whole_coords(), cell=[0.0] * 3)
        with pytest.raises(ValueError):
            make_frames_whole([f], unwrap_steps(_raw(tmp_path)))

    def test_none_passes_through(self, tmp_path):
        assert make_frames_whole(None, unwrap_steps(_raw(tmp_path))) is None


def _gro(coords):
    names = ["CHN"] * 4 + ["SOL"] * 2 + ["TRI"] * 3
    res = [1, 1, 1, 1, 2, 3, 4, 4, 4]
    lines = ["test", "%5d" % len(coords)]
    for i, (rn, r, xyz) in enumerate(zip(names, res, coords)):
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (r, rn, "X%d" % i, i + 1, *xyz))
    lines.append("%10.5f%10.5f%10.5f" % (BOX, BOX, BOX))
    return "\n".join(lines) + "\n"


def test_cli_needs_no_gmx_and_the_udf_is_whole(tmp_path, capsys):
    """効果を UDF で見る。gmx が無い指定でも通る (nojump を走らせないので)。"""
    pytest.importorskip("UDFManager")
    from UDFManager import UDFManager

    from abmptools.gro2udf.cli import main

    split0 = _wrapped(_whole_coords(shift=1.5))
    split1 = _wrapped(_whole_coords(shift=3.1))
    (tmp_path / "w.top").write_text(TOP)
    (tmp_path / "w.gro").write_text(_gro(split0))
    (tmp_path / "t.gro").write_text(_gro(split0) + _gro(split1))
    out = tmp_path / "w.udf"
    main(["gro2udf", "--from-top", str(tmp_path / "w.top"), str(tmp_path / "w.gro"),
          "--trajectory", str(tmp_path / "t.gro"), "--make-whole",
          "--gmx", "/nonexistent/gmx", "--keep-molecules", "CHN,TRI",
          "--constraints-as-bonds", "--ff", "", "--out", str(out)])
    assert "-pbc nojump is not run" in capsys.readouterr().out

    u = UDFManager(str(out))
    assert u.get("Set_of_Molecules.molecule[].Mol_Name") == ["CHN", "TRI"]
    assert u.totalRecord() == 2
    for rec in range(2):
        u.jump(rec)
        pos = u.get("Structure.Position.mol[].atom[]")
        assert pos is not None
        chn, tri = np.asarray(pos[0]), np.asarray(pos[1])        # A
        assert np.linalg.norm(chn[1:] - chn[:-1], axis=1) == pytest.approx([4.0] * 3, abs=0.02)
        assert max(np.linalg.norm(tri[i] - tri[j]) for i, j in [(0, 1), (0, 2), (1, 2)]) < 3.1
