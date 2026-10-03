# -*- coding: utf-8 -*-
"""
gro2udf ``--keep-molecules``: 指定した分子種だけを UDF に残す。

溶媒和した系はほとんどが溶媒で、粗視化の系では粒子の 3/4 が水ということもある。
全フレーム・全粒子の UDF はビューアで開くのが重い。溶媒を抜くには、これまで
``.top`` を書き換え、``.gro`` / ``.xtc`` もそれに合わせて切り出すしかなく、
他の人は同じ手順を手で再現するしかなかった。

ここでは読み込み時に選ぶ。``.top`` / ``.gro`` / 軌跡は元のまま。
``[ molecules ]`` の並びでどの原子がどのインスタンスかが決まるので、同じ
番号の並びで全フレームを切る。
"""
import pytest

from abmptools.gro2udf.cli import _BUILTIN_TEMPLATE, _from_top_parser, _resolve_options
from abmptools.gro2udf.molecule_select import select_molecules, subset_frames
from abmptools.gro2udf.top_model import GROFrameData
from abmptools.gro2udf.top_parser import TopParser

#: 溶質 2 種 (SOL 以外) が溶媒に挟まれて並ぶ系。並びの途中で種が入れ替わる
#: ので、オフセットの数え違いがあれば座標がずれる。
TOP = """\
[ defaults ]
1  2  no  1.0  1.0

[ atomtypes ]
  PA   45.0 0.000 A 0.43 2.0
  W    72.0 0.000 A 0.47 3.0

[ moleculetype ]
  POL  1

[ atoms ]
  1 PA 1 POL A1 1 0.0 45.0
  2 PA 1 POL A2 2 0.0 45.0
  3 PA 1 POL A3 3 0.0 45.0

[ bonds ]
  1 2 1 0.33 7000
  2 3 1 0.33 7000

[ moleculetype ]
  SOL  1

[ atoms ]
  1 W 1 SOL W 1 0.0 72.0

[ moleculetype ]
  LIG  1

[ atoms ]
  1 PA 1 LIG L1 1 0.0 45.0
  2 PA 1 LIG L2 2 0.0 45.0

[ bonds ]
  1 2 1 0.30 7000

[ system ]
  test

[ molecules ]
  POL  1
  SOL  3
  LIG  1
  SOL  2
  POL  1
"""
# 原子の並び (0 始まり):
#   POL 0-2 | SOL 3,4,5 | LIG 6-7 | SOL 8,9 | POL 10-12   (計 13)
N_FULL = 13


def _raw(tmp_path):
    p = tmp_path / "s.top"
    p.write_text(TOP)
    return TopParser().parse(str(p))


class TestSelect:
    def test_indices_follow_the_molecules_order(self, tmp_path):
        sel, idx = select_molecules(_raw(tmp_path), ["POL", "LIG"])
        assert idx == [0, 1, 2, 6, 7, 10, 11, 12]
        assert sel.mol_instance_list == ["POL", "LIG", "POL"]

    def test_order_of_the_names_does_not_matter(self, tmp_path):
        a = select_molecules(_raw(tmp_path), ["POL", "LIG"])[1]
        b = select_molecules(_raw(tmp_path), ["LIG", "POL", "LIG"])[1]
        assert a == b

    def test_original_raw_is_untouched(self, tmp_path):
        raw = _raw(tmp_path)
        select_molecules(raw, ["LIG"])
        assert raw.mol_instance_list.count("SOL") == 5

    def test_unknown_name_stops_and_lists_what_is_there(self, tmp_path):
        """打ち間違いで種を黙って落とさない。"""
        with pytest.raises(ValueError) as e:
            select_molecules(_raw(tmp_path), ["POL", "LGI"])
        msg = str(e.value)
        assert "LGI" in msg and "POL x2" in msg and "SOL x5" in msg


class TestSubsetFrames:
    def _frame(self, n):
        return GROFrameData(step=0, time=1.5,
                            coord_list=[[float(i), 0.0, 0.0] for i in range(n)],
                            cell=[5.0, 5.0, 5.0])

    def test_every_frame_is_cut_the_same_way(self):
        out = subset_frames([self._frame(N_FULL)] * 3, [0, 6, 12], N_FULL)
        assert [[c[0] for c in f.coord_list] for f in out] == [[0, 6, 12]] * 3
        assert out[0].time == 1.5 and out[0].cell == [5.0, 5.0, 5.0]

    def test_a_frame_of_another_system_is_refused(self):
        with pytest.raises(ValueError):
            subset_frames([self._frame(N_FULL + 1)], [0], N_FULL)

    def test_none_passes_through(self):
        assert subset_frames(None, [0]) is None


class TestCli:
    @staticmethod
    def _keep(argv):
        parser = _from_top_parser()
        args = parser.parse_args(argv)
        return args, _resolve_options(parser, args)[1]

    @pytest.mark.parametrize("argv", [
        ["--keep-molecules", "POL,LIG", "s.top", "s.gro"],
        ["s.top", "s.gro", "--keep-molecules", "POL,LIG"],
        ["--keep-molecules", "POL", "--keep-molecules", "LIG", "s.top", "s.gro"],
        ["s.top", "s.gro", "--keep-molecules", " POL , LIG "],
    ])
    def test_comma_or_repeat_and_any_placement(self, argv):
        """空白区切りで複数取る形だと、後ろの .top / .gro まで分子名として食べる。"""
        args, keep = self._keep(argv)
        assert (args.top_path, args.gro_path) == ("s.top", "s.gro")
        assert keep == ["POL", "LIG"]

    def test_absent_keeps_everything(self):
        assert self._keep(["s.top", "s.gro"])[1] is None

    def test_empty_is_an_error(self, capsys):
        with pytest.raises(SystemExit) as e:
            self._keep(["s.top", "s.gro", "--keep-molecules", ","])
        assert e.value.code == 2


def _gro(n_frames_offset=0.0):
    names = (["POL"] * 3 + ["SOL"] * 3 + ["LIG"] * 2 + ["SOL"] * 2 + ["POL"] * 3)
    res = [1, 1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 8, 8]
    lines = ["test", "%5d" % N_FULL]
    for i, (rn, r) in enumerate(zip(names, res)):
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                     % (r, rn, "X%d" % i, i + 1, 0.1 * i + n_frames_offset, 1.0, 1.0))
    lines.append("%10.5f%10.5f%10.5f" % (5.0, 5.0, 5.0))
    return "\n".join(lines) + "\n"


def test_the_udf_holds_only_the_kept_molecules(tmp_path):
    """効果を UDF で見る: 分子の並び、原子数、座標が元の該当原子と一致する。"""
    pytest.importorskip("UDFManager")
    from UDFManager import UDFManager

    from abmptools.gro2udf.top_exporter import TopExporter

    top = tmp_path / "s.top"
    top.write_text(TOP)
    gro = tmp_path / "s.gro"
    gro.write_text(_gro())
    traj = tmp_path / "t.gro"
    traj.write_text(_gro(0.0) + _gro(0.5))          # 2 フレーム
    out = tmp_path / "s.udf"
    TopExporter().export(str(top), str(gro), _BUILTIN_TEMPLATE, str(out),
                         force_field="", trajectory_path=str(traj),
                         keep_molecules=["POL", "LIG"])
    u = UDFManager(str(out))
    assert u.get("Set_of_Molecules.molecule[].Mol_Name") == ["POL", "LIG", "POL"]
    assert u.totalRecord() == 2
    want = [0, 1, 2, 6, 7, 10, 11, 12]
    for rec, shift in ((0, 0.0), (1, 0.5)):
        u.jump(rec)
        pos = u.get("Structure.Position.mol[].atom[]")
        assert pos is not None                      # パス違いの None 同士一致を防ぐ
        xs = [p[0] for mol in pos for p in mol]     # A
        assert xs == pytest.approx([(0.1 * i + shift) * 10.0 for i in want])
