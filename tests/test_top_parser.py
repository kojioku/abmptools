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


# ---------------------------------------------------------------------------
# End-of-line comments, in every section
# ---------------------------------------------------------------------------

#: The ``[ defaults ]`` of Martini 3's martini_v3.0.0.itp, verbatim: two
#: fields (the optional gen-pairs / fudgeLJ / fudgeQQ left out) and a comment.
MARTINI3_DEFAULTS = "[ defaults ]\n1 2 ; sigma-epsilon format of LJ parameters\n"

COMMENTED_TOP = MARTINI3_DEFAULTS + """
[ atomtypes ]
; name mass charge ptype sigma epsilon
P5  72.0 0.000 A 0.470 4.99 ; polar bead
C1  72.0 0.000 A 0.470 2.30;no space before the comment

[ bondtypes ]
P5 C1 1 0.47 1250.0 ; a bond type

[ angletypes ]
P5 C1 P5 2 120.0 25.0 ; GROMOS-96 angle

[ moleculetype ]
; molname nrexcl
PEP 1 ; one molecule

[ atoms ]
; nr type resnr residue atom cgnr charge
1 P5 1 ALA BB 1 0.0 ; no mass column, then a comment
2 C1 1 ALA SC1 2 0.0 72.0 ; mass and a comment
3 P5 2 GLY BB 3 0.0

[ bonds ]
1 2 1 0.47 1250.0 ; backbone-side chain
2 3 1 0.47 1250.0

[ angles ]
1 2 3 2 120.0 25.0 ; comment

[ system ]
test

[ molecules ]
PEP 2 ; two copies
"""


class TestEndOfLineComments:

    def test_martini3_defaults_line(self, tmp_path):
        """Martini 3's own [ defaults ] used to stop the parser:
        ValueError: could not convert string to float: 'sigma-epsilon'."""
        top = tmp_path / "m3.top"
        top.write_text(MARTINI3_DEFAULTS)
        raw = TopParser().parse(str(top))
        assert raw.comb_rule == 2
        # GROMACS' defaults for the optional fields
        assert raw.fudge_lj == pytest.approx(1.0)
        assert raw.fudge_qq == pytest.approx(1.0)

    def test_comments_after_data_in_every_section(self, tmp_path):
        top = tmp_path / "commented.top"
        top.write_text(COMMENTED_TOP)
        raw = TopParser().parse(str(top))
        assert [a[0] for a in raw.atomtypes] == ["P5", "C1"]
        assert raw.atomtypes[1][4] == pytest.approx(2.30)      # "2.30;no space"
        assert raw.mol_types == ["PEP"]
        atoms = raw.atomlist[0]
        assert [a[0] for a in atoms] == ["BB", "SC1", "BB"]
        assert atoms[0][3] == 0.0 and atoms[1][3] == pytest.approx(72.0)
        assert len(raw.bondlist[0]) == 2
        assert raw.mol_instance_list == ["PEP", "PEP"]

    def test_the_comment_does_not_change_the_numbers(self, tmp_path):
        """The same topology with and without the comments reads the same."""
        # drop the comments after data; whole-line comments stay (a blank line
        # would end the section in this parser, as in GROMACS' own files)
        plain = "\n".join(line if line.lstrip().startswith(";") else line.split(";", 1)[0].rstrip()
                          for line in COMMENTED_TOP.splitlines()) + "\n"
        a, b = tmp_path / "a.top", tmp_path / "b.top"
        a.write_text(COMMENTED_TOP)
        b.write_text(plain)
        ra, rb = TopParser().parse(str(a)), TopParser().parse(str(b))
        for field in ("comb_rule", "fudge_lj", "fudge_qq", "atomtypes", "mol_types",
                      "atomlist", "bondlist", "mol_instance_list"):
            assert getattr(ra, field) == getattr(rb, field), field
