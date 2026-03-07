import pytest
from abmptools.pdb_io import pdb_io


class TestPdbIoInit:
    """Tests for pdb_io.__init__ defaults."""

    def test_default_solutes(self):
        obj = pdb_io()
        assert obj.solutes == []

    def test_default_getmode(self):
        obj = pdb_io()
        assert obj.getmode == "resnum"

    def test_default_assignresname(self):
        obj = pdb_io()
        assert obj.assignresname is False

    def test_default_refreshatmtype(self):
        obj = pdb_io()
        assert obj.refreshatmtype is False

    def test_default_refreshresid(self):
        obj = pdb_io()
        assert obj.refreshresid is False

    def test_default_cellsize(self):
        obj = pdb_io()
        assert obj.cellsize == 0

    def test_default_empty_lists(self):
        obj = pdb_io()
        for attr in [
            "posRes",
            "atmtypeRes",
            "headRes",
            "labRes",
            "chainRes",
            "resnumRes",
            "totalRes",
            "resnames",
            "codeRes",
            "occRes",
            "tempRes",
            "amarkRes",
            "chargeRes",
            "rescount",
            "gatmlabRes",
        ]:
            assert getattr(obj, attr) == [], f"{attr} should be an empty list"


class TestReadpdb:
    """Tests for pdb_io.readpdb with a minimal PDB file."""

    PDB_CONTENT = (
        "COMPND    test\n"
        "AUTHOR    test\n"
        "HETATM    1  O   UNK     1       0.000   0.000   0.117  1.00  0.00           O\n"
        "HETATM    2  H   UNK     1       0.000   0.757  -0.469  1.00  0.00           H\n"
        "HETATM    3  H   UNK     1       0.000  -0.757  -0.469  1.00  0.00           H\n"
        "END\n"
    )

    @pytest.fixture()
    def pdb_obj(self, tmp_path):
        """Create a pdb_io object with the minimal PDB file read in."""
        pdb_file = tmp_path / "water.pdb"
        pdb_file.write_text(self.PDB_CONTENT)
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        return obj

    def test_totalRes_is_one(self, pdb_obj):
        """All three atoms share resnum 1, so there is exactly one residue."""
        assert pdb_obj.totalRes == 1

    def test_resnames(self, pdb_obj):
        assert pdb_obj.resnames == ["UNK"]

    def test_positions_shape(self, pdb_obj):
        """posRes should be a list of one residue containing three atom positions."""
        assert len(pdb_obj.posRes) == 1
        assert len(pdb_obj.posRes[0]) == 3

    def test_positions_values(self, pdb_obj):
        positions = pdb_obj.posRes[0]
        assert positions[0] == pytest.approx([0.0, 0.0, 0.117], abs=1e-3)
        assert positions[1] == pytest.approx([0.0, 0.757, -0.469], abs=1e-3)
        assert positions[2] == pytest.approx([0.0, -0.757, -0.469], abs=1e-3)

    def test_atom_types(self, pdb_obj):
        """atmtypeRes should contain the raw 4-character atom name fields."""
        atypes = pdb_obj.atmtypeRes[0]
        assert len(atypes) == 3
        # Atom names are parsed from columns 13-16 (0-indexed 12:16)
        assert atypes[0].strip() == "O"
        assert atypes[1].strip() == "H"
        assert atypes[2].strip() == "H"

    def test_head_is_hetatm(self, pdb_obj):
        """All records in the test file are HETATM."""
        for head in pdb_obj.headRes[0]:
            assert head.strip() == "HETATM"

    def test_chain_and_resnum(self, pdb_obj):
        """All atoms belong to the same chain and residue number."""
        chains = pdb_obj.chainRes[0]
        resnums = pdb_obj.resnumRes[0]
        assert len(set(chains)) == 1
        assert len(set(resnums)) == 1

    def test_lab_list_length(self, pdb_obj):
        assert len(pdb_obj.labRes) == 1
        assert len(pdb_obj.labRes[0]) == 3

    def test_rescount(self, pdb_obj):
        """rescount is a Counter over resnames."""
        assert pdb_obj.rescount["UNK"] == 1


class TestReadpdbMultipleResidues:
    """Test readpdb with two residues distinguished by residue number."""

    PDB_CONTENT = (
        "COMPND    test\n"
        "AUTHOR    test\n"
        "HETATM    1  O   UNK     1       0.000   0.000   0.117  1.00  0.00           O\n"
        "HETATM    2  H   UNK     1       0.000   0.757  -0.469  1.00  0.00           H\n"
        "HETATM    3  H   UNK     1       0.000  -0.757  -0.469  1.00  0.00           H\n"
        "HETATM    4  O   UNK     2       3.000   0.000   0.117  1.00  0.00           O\n"
        "HETATM    5  H   UNK     2       3.000   0.757  -0.469  1.00  0.00           H\n"
        "HETATM    6  H   UNK     2       3.000  -0.757  -0.469  1.00  0.00           H\n"
        "END\n"
    )

    @pytest.fixture()
    def pdb_obj(self, tmp_path):
        pdb_file = tmp_path / "two_waters.pdb"
        pdb_file.write_text(self.PDB_CONTENT)
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        return obj

    def test_totalRes_is_two(self, pdb_obj):
        assert pdb_obj.totalRes == 2

    def test_resnames_both_unk(self, pdb_obj):
        assert pdb_obj.resnames == ["UNK", "UNK"]

    def test_positions_per_residue(self, pdb_obj):
        assert len(pdb_obj.posRes) == 2
        assert len(pdb_obj.posRes[0]) == 3
        assert len(pdb_obj.posRes[1]) == 3

    def test_second_residue_x_offset(self, pdb_obj):
        """The second residue's oxygen should be at x=3.0."""
        assert pdb_obj.posRes[1][0][0] == pytest.approx(3.0, abs=1e-3)

    def test_atmtypeRes_per_residue(self, pdb_obj):
        assert len(pdb_obj.atmtypeRes) == 2
        for res_atypes in pdb_obj.atmtypeRes:
            stripped = [a.strip() for a in res_atypes]
            assert stripped == ["O", "H", "H"]

    def test_rescount_two(self, pdb_obj):
        assert pdb_obj.rescount["UNK"] == 2


class TestReadpdbAtomRecord:
    """Test that ATOM records (not just HETATM) are parsed."""

    PDB_CONTENT = (
        "COMPND    test\n"
        "AUTHOR    test\n"
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n"
        "ATOM      2  N   ALA A   1       1.500   2.500   3.500  1.00 10.00           N\n"
        "END\n"
    )

    @pytest.fixture()
    def pdb_obj(self, tmp_path):
        pdb_file = tmp_path / "atom_record.pdb"
        pdb_file.write_text(self.PDB_CONTENT)
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        return obj

    def test_totalRes(self, pdb_obj):
        assert pdb_obj.totalRes == 1

    def test_resnames(self, pdb_obj):
        assert pdb_obj.resnames == ["ALA"]

    def test_head_is_atom(self, pdb_obj):
        for head in pdb_obj.headRes[0]:
            assert head.strip() == "ATOM"

    def test_positions(self, pdb_obj):
        assert pdb_obj.posRes[0][0] == pytest.approx([1.0, 2.0, 3.0], abs=1e-3)
        assert pdb_obj.posRes[0][1] == pytest.approx([1.5, 2.5, 3.5], abs=1e-3)

    def test_chain(self, pdb_obj):
        assert pdb_obj.chainRes[0][0] == "A"
