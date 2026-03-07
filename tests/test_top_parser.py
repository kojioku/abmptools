# -*- coding: utf-8 -*-
"""Tests for abmptools.gro2udf.top_parser."""
import pytest

from abmptools.gro2udf.top_parser import (
    TopParser,
    get_bond_name,
    get_angle_name,
    get_torsion_name,
    is_improper,
    _is_comment,
    _is_blank,
)


# ---------------------------------------------------------------------------
# Pure helper function tests
# ---------------------------------------------------------------------------

class TestGetBondName:
    def test_alphabetical_order(self):
        assert get_bond_name("ca", "ha") == "ca-ha"

    def test_reversed_order(self):
        assert get_bond_name("ha", "ca") == "ca-ha"

    def test_same_names(self):
        assert get_bond_name("ca", "ca") == "ca-ca"


class TestGetAngleName:
    def test_forward(self):
        assert get_angle_name("ca", "cb", "ha") == "ca-cb-ha"

    def test_reversed(self):
        assert get_angle_name("ha", "cb", "ca") == "ca-cb-ha"


class TestGetTorsionName:
    def test_forward(self):
        assert get_torsion_name("a", "b", "c", "d") == "a-b-c-d"

    def test_reversed(self):
        assert get_torsion_name("d", "c", "b", "a") == "a-b-c-d"

    def test_same_outer_tiebreak(self):
        # n1 == n4, so compare n2 vs n3
        assert get_torsion_name("a", "c", "b", "a") == "a-b-c-a"


class TestIsImproper:
    def test_proper_torsion(self):
        # chain: 1-2, 2-3, 3-4 all bonded
        bondlist = [[1, 2, 0], [2, 3, 0], [3, 4, 0]]
        assert is_improper([1, 2, 3, 4], bondlist) is False

    def test_improper_torsion(self):
        # missing bond 2-3
        bondlist = [[1, 2, 0], [3, 4, 0]]
        assert is_improper([1, 2, 3, 4], bondlist) is True


class TestLineClassification:
    def test_comment_semicolon(self):
        assert _is_comment("; this is a comment") is True

    def test_comment_hash(self):
        assert _is_comment("#define FOO") is True

    def test_not_comment(self):
        assert _is_comment("1  ca  1  MOL  ca1  1  -0.1  12.011") is False

    def test_blank(self):
        assert _is_blank("   \n") is True

    def test_not_blank(self):
        assert _is_blank("hello") is False

    def test_empty_string_is_blank_not_comment(self):
        assert _is_blank("") is True
        assert _is_comment("") is False


# ---------------------------------------------------------------------------
# TopParser.parse() integration tests with a minimal .top file
# ---------------------------------------------------------------------------

MINIMAL_TOP = """\
[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1   2   yes   0.5   0.8333

[ atomtypes ]
; name mass charge ptype sigma epsilon
ca   12.011  0.0  A  0.339967  0.359824
ha    1.008  0.0  A  0.259964  0.062760

[ moleculetype ]
; Name nrexcl
MOL  3

[ atoms ]
; nr type resnr residu atom cgnr charge mass
1  ca  1  MOL  ca1  1  -0.1  12.011
2  ha  1  MOL  ha1  2   0.1   1.008

[ bonds ]
; ai aj funct b0 kb
1  2  1  0.14  300000.0

[ molecules ]
MOL  2
"""


class TestTopParserParse:

    @pytest.fixture
    def raw(self, tmp_path):
        top = tmp_path / "test.top"
        top.write_text(MINIMAL_TOP)
        return TopParser().parse(str(top))

    def test_defaults(self, raw):
        assert raw.comb_rule == 2
        assert raw.fudge_lj == pytest.approx(0.5)
        assert raw.fudge_qq == pytest.approx(0.8333)

    def test_atomtypes(self, raw):
        assert len(raw.atomtypes) == 2
        # [name, mass, charge, sigma, epsilon]
        ca = raw.atomtypes[0]
        assert ca[0] == "ca"
        assert ca[1] == pytest.approx(12.011)

    def test_moleculetype_and_atoms(self, raw):
        assert raw.mol_types == ["MOL"]
        assert len(raw.atomlist) == 1
        assert len(raw.atomlist[0]) == 2  # 2 atoms in MOL
        # each atom: [atom_name, type_name, charge, mass]
        assert raw.atomlist[0][0][0] == "ca1"
        assert raw.atomlist[0][0][1] == "ca"

    def test_bonds(self, raw):
        assert len(raw.bondlist) == 1
        assert len(raw.bondlist[0]) == 1  # 1 bond in MOL
        bond = raw.bondlist[0][0]
        assert bond[0] == 1  # atom1
        assert bond[1] == 2  # atom2

    def test_molecules_section(self, raw):
        # "MOL 2" means 2 instances
        assert raw.mol_instance_list == ["MOL", "MOL"]
