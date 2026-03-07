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


class TestGetpdbcell:
    """Tests for pdb_io.getpdbcell which reads CRYST1 records."""

    def test_getpdbcell_reads_cryst1(self, tmp_path):
        pdb_file = tmp_path / "cell.pdb"
        pdb_file.write_text(
            "CRYST1   38.376   40.123   42.567  90.00  90.00  90.00               1\n"
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n"
            "END\n"
        )
        obj = pdb_io()
        cell = obj.getpdbcell(str(pdb_file))
        assert cell == pytest.approx([38.376, 40.123, 42.567], abs=1e-3)

    def test_getpdbcell_no_cryst1(self, tmp_path):
        """When no CRYST1 record exists, an empty list is returned."""
        pdb_file = tmp_path / "nocell.pdb"
        pdb_file.write_text(
            "COMPND    test\n"
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n"
            "END\n"
        )
        obj = pdb_io()
        cell = obj.getpdbcell(str(pdb_file))
        assert cell == []

    def test_getpdbcell_cryst_after_atoms(self, tmp_path):
        """CRYST1 can appear anywhere; getpdbcell should still find it."""
        pdb_file = tmp_path / "cell_late.pdb"
        pdb_file.write_text(
            "COMPND    test\n"
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n"
            "CRYST1   10.000   20.000   30.000  90.00  90.00  90.00               1\n"
            "END\n"
        )
        obj = pdb_io()
        cell = obj.getpdbcell(str(pdb_file))
        assert cell == pytest.approx([10.0, 20.0, 30.0], abs=1e-3)


class TestExportardpdbfull:
    """Tests for pdb_io.exportardpdbfull which writes a PDB from internal state."""

    PDB_CONTENT = (
        "COMPND    test\n"
        "AUTHOR    test\n"
        "HETATM    1  O   WAT     1       1.000   2.000   3.000  1.00  0.00           O\n"
        "HETATM    2  H   WAT     1       1.500   2.500   3.500  1.00  0.00           H\n"
        "HETATM    3  O   WAT     2       4.000   5.000   6.000  1.00  0.00           O\n"
        "HETATM    4  H   WAT     2       4.500   5.500   6.500  1.00  0.00           H\n"
        "END\n"
    )

    @pytest.fixture()
    def loaded_obj(self, tmp_path):
        pdb_file = tmp_path / "input.pdb"
        pdb_file.write_text(self.PDB_CONTENT)
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        return obj

    def test_export_creates_file(self, loaded_obj, tmp_path):
        outpath = str(tmp_path / "out")
        loaded_obj.exportardpdbfull(outpath, [0, 1])
        assert (tmp_path / "out.pdb").exists()

    def test_export_contains_end(self, loaded_obj, tmp_path):
        outpath = str(tmp_path / "out")
        loaded_obj.exportardpdbfull(outpath, [0, 1])
        content = (tmp_path / "out.pdb").read_text()
        assert "END" in content

    def test_export_contains_compnd(self, loaded_obj, tmp_path):
        outpath = str(tmp_path / "out")
        loaded_obj.exportardpdbfull(outpath, [0, 1])
        content = (tmp_path / "out.pdb").read_text()
        assert "COMPND" in content

    def test_export_all_atoms_written(self, loaded_obj, tmp_path):
        outpath = str(tmp_path / "out")
        loaded_obj.exportardpdbfull(outpath, [0, 1])
        content = (tmp_path / "out.pdb").read_text()
        hetatm_lines = [l for l in content.splitlines() if l.startswith("HETATM")]
        assert len(hetatm_lines) == 4

    def test_export_subset_of_molecules(self, loaded_obj, tmp_path):
        """Exporting only residue 0 should produce 2 HETATM lines."""
        outpath = str(tmp_path / "subset")
        loaded_obj.exportardpdbfull(outpath, [0])
        content = (tmp_path / "subset.pdb").read_text()
        hetatm_lines = [l for l in content.splitlines() if l.startswith("HETATM")]
        assert len(hetatm_lines) == 2

    def test_export_with_cellsize(self, loaded_obj, tmp_path):
        """When cellsize is set, a CRYST1 record should appear."""
        loaded_obj.cellsize = [10.0, 20.0, 30.0]
        outpath = str(tmp_path / "withcell")
        loaded_obj.exportardpdbfull(outpath, [0, 1])
        content = (tmp_path / "withcell.pdb").read_text()
        assert "CRYST1" in content
        assert "10.000" in content
        assert "20.000" in content
        assert "30.000" in content

    def test_export_no_cryst1_without_cellsize(self, loaded_obj, tmp_path):
        """When cellsize is 0 (default), no CRYST1 record appears."""
        loaded_obj.cellsize = 0
        outpath = str(tmp_path / "nocell")
        loaded_obj.exportardpdbfull(outpath, [0])
        content = (tmp_path / "nocell.pdb").read_text()
        assert "CRYST1" not in content


class TestExportardpdb:
    """Tests for pdb_io.exportardpdb."""

    def test_export_creates_pdb_file(self, tmp_path):
        outpath = str(tmp_path / "result.pdb")
        obj = pdb_io()
        posRes = [[[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]]]
        nameAtom = [["O", "H"]]
        molnames_orig = ["WAT"]
        obj.exportardpdb(outpath, [0], posRes, nameAtom, molnames_orig)
        assert (tmp_path / "result.pdb").exists()

    def test_export_atom_count(self, tmp_path):
        outpath = str(tmp_path / "result.pdb")
        obj = pdb_io()
        posRes = [[[1.0, 2.0, 3.0]], [[4.0, 5.0, 6.0]]]
        nameAtom = [["O"], ["N"]]
        molnames_orig = ["WAT", "NH3"]
        obj.exportardpdb(outpath, [0, 1], posRes, nameAtom, molnames_orig)
        content = (tmp_path / "result.pdb").read_text()
        hetatm_lines = [l for l in content.splitlines() if l.startswith("HETATM")]
        assert len(hetatm_lines) == 2

    def test_export_has_end(self, tmp_path):
        outpath = str(tmp_path / "result.pdb")
        obj = pdb_io()
        posRes = [[[0.0, 0.0, 0.0]]]
        nameAtom = [["C"]]
        molnames_orig = ["MOL"]
        obj.exportardpdb(outpath, [0], posRes, nameAtom, molnames_orig)
        content = (tmp_path / "result.pdb").read_text()
        assert content.strip().endswith("END")


class TestReadpdbEdgeCases:
    """Edge cases and special scenarios for readpdb."""

    def test_mixed_atom_hetatm(self, tmp_path):
        """File with both ATOM and HETATM records in the same residue."""
        pdb_file = tmp_path / "mixed.pdb"
        pdb_file.write_text(
            "COMPND    test\n"
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n"
            "HETATM    2  O   ALA A   1       1.500   2.500   3.500  1.00  0.00           O\n"
            "END\n"
        )
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        assert obj.totalRes == 1
        assert len(obj.posRes[0]) == 2
        assert obj.headRes[0][0].strip() == "ATOM"
        assert obj.headRes[0][1].strip() == "HETATM"

    def test_multiple_different_resnames(self, tmp_path):
        """Multiple residues with different residue names."""
        pdb_file = tmp_path / "multi.pdb"
        pdb_file.write_text(
            "COMPND    test\n"
            "HETATM    1  O   WAT     1       0.000   0.000   0.000  1.00  0.00           O\n"
            "HETATM    2  N   NH3     2       5.000   5.000   5.000  1.00  0.00           N\n"
            "END\n"
        )
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        assert obj.totalRes == 2
        assert obj.resnames == ["WAT", "NH3"]

    def test_cryst1_line_is_skipped_in_readpdb(self, tmp_path):
        """CRYST1 lines should be ignored by readpdb (not treated as atoms)."""
        pdb_file = tmp_path / "cryst.pdb"
        pdb_file.write_text(
            "CRYST1   38.376   38.376   38.376  90.00  90.00  90.00               1\n"
            "HETATM    1  O   WAT     1       1.000   2.000   3.000  1.00  0.00           O\n"
            "END\n"
        )
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        assert obj.totalRes == 1
        assert len(obj.posRes[0]) == 1

    def test_readpdb_occ_and_temp(self, tmp_path):
        """Occupancy and temperature factor fields are parsed."""
        pdb_file = tmp_path / "occ.pdb"
        pdb_file.write_text(
            "COMPND    test\n"
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  0.75 15.30           C\n"
            "END\n"
        )
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        assert obj.occRes[0][0].strip() == "0.75"
        assert obj.tempRes[0][0].strip() == "15.30"

    def test_readpdb_amark_element(self, tmp_path):
        """Element symbol (amark) is parsed from columns 77-78."""
        pdb_file = tmp_path / "elem.pdb"
        pdb_file.write_text(
            "COMPND    test\n"
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C\n"
            "END\n"
        )
        obj = pdb_io()
        obj.readpdb(str(pdb_file))
        assert obj.amarkRes[0][0].strip() == "C"


class TestGetpdbinfowrap:
    """Tests for getpdbinfowrap which combines readpdb and getpdbcell."""

    def test_wrap_with_cryst1(self, tmp_path):
        pdb_file = tmp_path / "wrap.pdb"
        pdb_file.write_text(
            "CRYST1   25.000   30.000   35.000  90.00  90.00  90.00               1\n"
            "HETATM    1  O   WAT     1       1.000   2.000   3.000  1.00  0.00           O\n"
            "END\n"
        )
        obj = pdb_io()
        result = obj.getpdbinfowrap(str(pdb_file))
        assert result.cellsize == pytest.approx([25.0, 30.0, 35.0], abs=1e-3)
        assert result.totalRes == 1

    def test_wrap_without_cryst1(self, tmp_path):
        """Without CRYST1, cellsize should be set to 0."""
        pdb_file = tmp_path / "nowrap.pdb"
        pdb_file.write_text(
            "COMPND    test\n"
            "HETATM    1  O   WAT     1       1.000   2.000   3.000  1.00  0.00           O\n"
            "END\n"
        )
        obj = pdb_io()
        result = obj.getpdbinfowrap(str(pdb_file))
        assert result.cellsize == 0


class TestMovemoltranspdb:
    """Tests for movemoltranspdb (parallel shift)."""

    def test_translate_vector(self):
        import numpy as np
        obj = pdb_io()
        pos = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        trans = np.array([10.0, 20.0, 30.0])
        result = obj.movemoltranspdb(pos, trans)
        expected = np.array([[11.0, 22.0, 33.0], [14.0, 25.0, 36.0]])
        assert result == pytest.approx(expected)


class TestRoundtripReadExport:
    """Test that reading a PDB and exporting it produces a valid PDB."""

    def test_roundtrip_preserves_atom_count(self, tmp_path):
        pdb_content = (
            "COMPND    test\n"
            "HETATM    1  O   WAT     1       1.000   2.000   3.000  1.00  0.00           O\n"
            "HETATM    2  H   WAT     1       1.500   2.500   3.500  1.00  0.00           H\n"
            "HETATM    3  H   WAT     1       0.500   1.500   2.500  1.00  0.00           H\n"
            "END\n"
        )
        infile = tmp_path / "in.pdb"
        infile.write_text(pdb_content)
        obj = pdb_io()
        obj.readpdb(str(infile))
        outpath = str(tmp_path / "roundtrip")
        obj.exportardpdbfull(outpath, [0])
        content = (tmp_path / "roundtrip.pdb").read_text()
        hetatm_lines = [l for l in content.splitlines() if l.startswith("HETATM")]
        assert len(hetatm_lines) == 3

    def test_roundtrip_preserves_positions(self, tmp_path):
        """Positions in the exported file should match the original."""
        pdb_content = (
            "COMPND    test\n"
            "HETATM    1  O   WAT     1      12.345  67.890  -1.234  1.00  0.00           O\n"
            "END\n"
        )
        infile = tmp_path / "in.pdb"
        infile.write_text(pdb_content)
        obj = pdb_io()
        obj.readpdb(str(infile))
        outpath = str(tmp_path / "roundtrip2")
        obj.exportardpdbfull(outpath, [0])
        # Re-read the exported file
        obj2 = pdb_io()
        obj2.readpdb(str(tmp_path / "roundtrip2.pdb"))
        assert obj2.posRes[0][0] == pytest.approx([12.345, 67.890, -1.234], abs=1e-3)
