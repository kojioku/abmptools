# -*- coding: utf-8 -*-
"""
gro2udf ``--constraints-as-bonds``: ``[ constraints ]`` を結合として読む。

``[ constraints ]`` は読まれないので、**拘束だけで結合を定義した分子は UDF の中で
原子が 1 本もつながらない**。Martini では環や小さな剛体部分がこの形で書かれる
(環を拘束だけで組んだ分子で実際に起きた)。
guard はこれを fatal にしていたが、直すには ``.top`` を書き換えるしかなく、
他の人が同じ手順で再現できなかった。

オプションを付けたときだけ、拘束を同じ長さの調和結合 (funct 1) として読む。
つながりと長さは正確に移る。拘束には力の定数が無いので、K は代用の値である。
"""
import pytest

from abmptools.gro2udf.cli import _BUILTIN_TEMPLATE, _from_top_parser
from abmptools.gro2udf.guard import check_top
from abmptools.gro2udf.top_parser import DEFAULT_CONSTRAINT_K, TopParser


#: 拘束だけで組んだ 4 ビーズの分子 (三角形の環 + 環の 2 点に付く置換基) と、
#: 結合と拘束が混ざった 3 ビーズの鎖。comb-rule 2 の普通の型を使い、
#: guard に引っかかるのが [ constraints ] だけになるようにしてある。
TOP = """\
[ defaults ]
1  2  no  1.0  1.0

[ atomtypes ]
  RA   45.0 0.000 A 0.43 2.0
  RB   72.0 0.000 A 0.47 3.0

[ moleculetype ]
; name nrexcl
  RING  1

[ atoms ]
  1 RA 1 RING A1 1 0.0 45.0
  2 RA 1 RING A2 2 0.0 45.0
  3 RA 1 RING A3 3 0.0 45.0
  4 RB 1 RING S4 4 0.0 72.0

[ constraints ]
;  i   j  funct  length
   1   2   1   0.2540   ; ring
   1   3   1   0.2540
   2   3   2   0.2560   ; funct 2 fixes the length too
   2   4   1   0.3686
   3   4   1   0.3686

[ moleculetype ]
  CHAIN  1

[ atoms ]
  1 RA 1 CHAIN C1 1 0.0 45.0
  2 RA 1 CHAIN C2 2 0.0 45.0
  3 RB 1 CHAIN C3 3 0.0 72.0

[ bonds ]
  1 2 1 0.33 7000

[ constraints ]
  2 3 1 0.30

[ system ]
  test

[ molecules ]
  RING   2
  CHAIN  1
"""


def _parse(tmp_path, **kw):
    p = tmp_path / "c.top"
    p.write_text(TOP)
    return TopParser(**kw).parse(str(p))


def _pairs(raw, mol):
    """moleculetype 番号 mol の (atom1, atom2, R0, K) を並び順で返す。"""
    out = []
    for a1, a2, tidx in raw.bondlist[mol]:
        _, _, funct, (b0, k) = raw.bond_types_from_mol[tidx]
        assert funct == 1
        out.append((a1, a2, b0, k))
    return out


class TestOffByDefault:
    """付けなければ今までどおり: 拘束は読まれず、guard が止める。"""

    def test_constraints_are_not_read(self, tmp_path):
        raw = _parse(tmp_path)
        assert raw.bondlist[0] == []                 # RING: 結合 0 本
        assert len(raw.bondlist[1]) == 1             # CHAIN: [ bonds ] の 1 本だけ
        assert raw.n_constraints_as_bonds == 0
        assert raw.sections["constraints"] == 6

    def test_guard_stops_and_names_the_option(self, tmp_path):
        raw = _parse(tmp_path)
        fatal, _ = check_top(raw, raw.sections, raw.dihedral_functs)
        assert len(fatal) == 1
        assert "[ constraints ] has 6 entries" in fatal[0]
        assert "--constraints-as-bonds" in fatal[0]


class TestReadAsBonds:
    def test_every_constraint_becomes_a_bond_of_the_same_length(self, tmp_path):
        raw = _parse(tmp_path, constraints_as_bonds=12345.0)
        assert _pairs(raw, 0) == [
            (1, 2, 0.2540, 12345.0),
            (1, 3, 0.2540, 12345.0),
            (2, 3, 0.2560, 12345.0),     # funct 2 も同じく読む
            (2, 4, 0.3686, 12345.0),
            (3, 4, 0.3686, 12345.0),
        ]
        assert raw.n_constraints_as_bonds == 6

    def test_bonds_and_constraints_in_one_molecule_both_survive(self, tmp_path):
        raw = _parse(tmp_path, constraints_as_bonds=12345.0)
        assert _pairs(raw, 1) == [(1, 2, 0.33, 7000.0), (2, 3, 0.30, 12345.0)]

    def test_the_guard_no_longer_stops(self, tmp_path):
        """拘束は結合として運ばれたので、失われるものは無い。"""
        raw = _parse(tmp_path, constraints_as_bonds=DEFAULT_CONSTRAINT_K)
        assert raw.sections["constraints"] == 0
        fatal, warnings = check_top(raw, raw.sections, raw.dihedral_functs)
        assert fatal == []
        assert not any("constraints" in w for w in warnings)

    def test_other_sections_still_parse(self, tmp_path):
        """[ constraints ] の後ろの [ moleculetype ] / [ molecules ] を取りこぼさない。"""
        raw = _parse(tmp_path, constraints_as_bonds=DEFAULT_CONSTRAINT_K)
        assert raw.mol_types == ["RING", "CHAIN"]
        assert raw.mol_instance_list == ["RING", "RING", "CHAIN"]
        assert [len(a) for a in raw.atomlist] == [4, 3]


class TestCli:
    """値を取らないフラグ + 別の --constraint-k。

    最初は ``--constraints-as-bonds [K]`` (値を省略可) だったが、
    ``--constraints-as-bonds system.top conf.gro`` の並びで .top のパスを K として
    食べ、「float として不正: system.top」で止まった。オプションを先に書く人は
    普通にいて、その文面からは原因にたどり着けない。
    """

    @staticmethod
    def _resolve(argv):
        from abmptools.gro2udf.cli import _resolve_options
        parser = _from_top_parser()
        args = parser.parse_args(argv)
        return args, _resolve_options(parser, args)

    @pytest.mark.parametrize("argv", [
        ["--constraints-as-bonds", "system.top", "conf.gro", "--out", "a.udf"],
        ["system.top", "conf.gro", "--constraints-as-bonds", "--out", "b.udf"],
        ["system.top", "conf.gro", "--out", "c.udf", "--constraints-as-bonds"],
    ])
    def test_any_placement_keeps_the_positionals(self, argv):
        args, (k, _) = self._resolve(argv)
        assert (args.top_path, args.gro_path) == ("system.top", "conf.gro")
        assert k == DEFAULT_CONSTRAINT_K

    def test_constraint_k_sets_the_value(self):
        _, (k, _) = self._resolve(["--constraints-as-bonds", "--constraint-k", "2e5",
                                   "system.top", "conf.gro"])
        assert k == 2e5

    def test_absent_means_off(self):
        _, (k, _) = self._resolve(["system.top", "conf.gro"])
        assert k is None

    def test_constraint_k_alone_is_an_error_not_ignored(self, capsys):
        with pytest.raises(SystemExit) as e:
            self._resolve(["system.top", "conf.gro", "--constraint-k", "2e5"])
        assert e.value.code == 2
        assert "--constraints-as-bonds" in capsys.readouterr().err


def _gro():
    """固定幅 (%5d%-5s%5s%5d%8.3f x3) で書く。"""
    atoms = [(1, "RING", "A1", 1.000, 1.000, 1.000), (1, "RING", "A2", 1.254, 1.000, 1.000),
             (1, "RING", "A3", 1.127, 1.220, 1.000), (1, "RING", "S4", 1.400, 1.300, 1.000),
             (2, "RING", "A1", 2.000, 2.000, 2.000), (2, "RING", "A2", 2.254, 2.000, 2.000),
             (2, "RING", "A3", 2.127, 2.220, 2.000), (2, "RING", "S4", 2.400, 2.300, 2.000),
             (3, "CHAIN", "C1", 3.000, 3.000, 3.000), (3, "CHAIN", "C2", 3.330, 3.000, 3.000),
             (3, "CHAIN", "C3", 3.630, 3.000, 3.000)]
    lines = ["test", "%5d" % len(atoms)]
    for n, (res, rname, aname, x, y, z) in enumerate(atoms, start=1):
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (res, rname, aname, n, x, y, z))
    lines.append("%10.5f%10.5f%10.5f" % (5.0, 5.0, 5.0))
    return "\n".join(lines) + "\n"


def test_the_udf_carries_the_bonds(tmp_path):
    """効果を UDF で見る: 拘束だけの分子にも結合が書かれている。"""
    pytest.importorskip("UDFManager")
    from UDFManager import UDFManager

    from abmptools.gro2udf.top_exporter import TopExporter

    top = tmp_path / "c.top"
    top.write_text(TOP)
    gro = tmp_path / "c.gro"
    gro.write_text(_gro())
    out = tmp_path / "c.udf"
    TopExporter().export(str(top), str(gro), _BUILTIN_TEMPLATE, str(out),
                         force_field="", constraints_as_bonds=DEFAULT_CONSTRAINT_K)
    u = UDFManager(str(out))
    names = u.get("Set_of_Molecules.molecule[].Mol_Name")
    a1 = u.get("Set_of_Molecules.molecule[].bond[].atom1")
    assert names == ["RING", "RING", "CHAIN"]
    assert [len(x) for x in a1] == [5, 5, 2]
    # None 同士の比較で通らないよう、値そのものを見る (パスを間違えると None が返る)
    r0 = u.get("Molecular_Attributes.Bond_Potential[].R0")
    assert r0 is not None
    assert sorted(set(round(x, 4) for x in r0)) == pytest.approx(
        [2.54, 2.56, 3.0, 3.3, 3.686])   # A (UDF 内の単位)
