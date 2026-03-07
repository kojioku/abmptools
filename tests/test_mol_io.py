import os
import pytest
from abmptools.mol_io import mol_io


WATER_XYZ = """\
3
water molecule
O   0.000  0.000  0.117
H   0.000  0.757 -0.469
H   0.000 -0.757 -0.469
"""

MULTI_XYZ = """\
3
water molecule
O   0.000  0.000  0.117
H   0.000  0.757 -0.469
H   0.000 -0.757 -0.469
2
H2 molecule
H   0.000  0.000  0.000
H   0.000  0.000  0.740
"""


@pytest.fixture
def mio():
    return mol_io()


@pytest.fixture
def water_xyz(tmp_path):
    fpath = tmp_path / "water.xyz"
    fpath.write_text(WATER_XYZ)
    return str(fpath)


@pytest.fixture
def multi_xyz(tmp_path):
    fpath = tmp_path / "multi.xyz"
    fpath.write_text(MULTI_XYZ)
    return str(fpath)


# ---------- read_mol_name ----------

class TestReadMolName:

    def test_simple_filename(self, mio):
        assert mio.read_mol_name("water.pdb") == "water"

    def test_xyz_extension(self, mio):
        assert mio.read_mol_name("molecule.xyz") == "molecule"

    def test_with_path(self, mio):
        assert mio.read_mol_name("path/to/mol.xyz") == "path/to/mol"

    def test_no_extension(self, mio):
        assert mio.read_mol_name("noext") == "noext"

    def test_multiple_dots(self, mio):
        assert mio.read_mol_name("my.molecule.pdb") == "my.molecule"


# ---------- read_xyz ----------

class TestReadXyz:

    def test_atom_list(self, mio, water_xyz):
        atom_list, coord_list = mio.read_xyz(water_xyz)
        assert atom_list == ["O", "H", "H"]

    def test_coord_list_length(self, mio, water_xyz):
        atom_list, coord_list = mio.read_xyz(water_xyz)
        assert len(coord_list) == 3

    def test_coord_values(self, mio, water_xyz):
        _, coord_list = mio.read_xyz(water_xyz)
        assert coord_list[0] == pytest.approx([0.0, 0.0, 0.117])
        assert coord_list[1] == pytest.approx([0.0, 0.757, -0.469])
        assert coord_list[2] == pytest.approx([0.0, -0.757, -0.469])

    def test_coords_are_floats(self, mio, water_xyz):
        _, coord_list = mio.read_xyz(water_xyz)
        for coord in coord_list:
            for val in coord:
                assert isinstance(val, float)


# ---------- getatoms (multi-structure) ----------

class TestGetAtoms:

    def test_returns_two_structures(self, mio, multi_xyz):
        atoms, atomnums, poss = mio.getatoms(multi_xyz)
        assert len(atoms) == 2
        assert len(atomnums) == 2
        assert len(poss) == 2

    def test_atomnums(self, mio, multi_xyz):
        _, atomnums, _ = mio.getatoms(multi_xyz)
        assert atomnums[0] == "3"
        assert atomnums[1] == "2"

    def test_first_structure_atoms(self, mio, multi_xyz):
        atoms, _, _ = mio.getatoms(multi_xyz)
        # Each atom is [index, symbol]
        symbols = [a[1] for a in atoms[0]]
        assert symbols == ["O", "H", "H"]

    def test_second_structure_atoms(self, mio, multi_xyz):
        atoms, _, _ = mio.getatoms(multi_xyz)
        symbols = [a[1] for a in atoms[1]]
        assert symbols == ["H", "H"]

    def test_positions_are_strings(self, mio, multi_xyz):
        _, _, poss = mio.getatoms(multi_xyz)
        # getatoms returns positions as strings
        for pos_list in poss:
            for pos in pos_list:
                assert len(pos) == 3
                for val in pos:
                    assert isinstance(val, str)

    def test_single_structure(self, mio, water_xyz):
        atoms, atomnums, poss = mio.getatoms(water_xyz)
        assert len(atoms) == 1
        assert atomnums[0] == "3"


# ---------- getatom (single-structure) ----------

class TestGetAtom:

    def test_atomnum(self, mio, water_xyz):
        _, atomnum, _ = mio.getatom(water_xyz)
        assert atomnum == "3"

    def test_atom_indices_and_symbols(self, mio, water_xyz):
        atom, _, _ = mio.getatom(water_xyz)
        assert len(atom) == 3
        assert atom[0] == [0, "O"]
        assert atom[1] == [1, "H"]
        assert atom[2] == [2, "H"]

    def test_positions_are_strings(self, mio, water_xyz):
        _, _, pos = mio.getatom(water_xyz)
        assert len(pos) == 3
        for p in pos:
            assert len(p) == 3
            for val in p:
                assert isinstance(val, str)

    def test_position_values(self, mio, water_xyz):
        _, _, pos = mio.getatom(water_xyz)
        assert float(pos[0][0]) == pytest.approx(0.0)
        assert float(pos[0][2]) == pytest.approx(0.117)


# ---------- Exportpospdb ----------

class TestExportpospdb:

    def test_writes_pdb_file(self, mio, tmp_path):
        out_file = str(tmp_path / "out.pdb")
        atom = [[0, "O"], [1, "H"], [2, "H"]]
        pos = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]
        mio.Exportpospdb(atom, "3", pos, out_file)

        assert os.path.exists(out_file)

    def test_pdb_contains_hetatm(self, mio, tmp_path):
        out_file = str(tmp_path / "out.pdb")
        atom = [[0, "O"], [1, "H"], [2, "H"]]
        pos = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]
        mio.Exportpospdb(atom, "3", pos, out_file)

        content = open(out_file).read()
        assert "HETATM" in content
        assert content.strip().endswith("END")

    def test_pdb_has_correct_atom_count(self, mio, tmp_path):
        out_file = str(tmp_path / "out.pdb")
        atom = [[0, "O"], [1, "H"], [2, "H"]]
        pos = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]
        mio.Exportpospdb(atom, "3", pos, out_file)

        content = open(out_file).readlines()
        hetatm_lines = [l for l in content if l.startswith("HETATM")]
        assert len(hetatm_lines) == 3

    def test_pdb_contains_element_symbols(self, mio, tmp_path):
        out_file = str(tmp_path / "out.pdb")
        atom = [[0, "O"], [1, "H"], [2, "H"]]
        pos = [[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]]
        mio.Exportpospdb(atom, "3", pos, out_file)

        content = open(out_file).read()
        assert " O" in content
        assert " H" in content

    def test_pdb_with_string_positions(self, mio, tmp_path):
        """Exportpospdb should handle string positions (from getatom)."""
        out_file = str(tmp_path / "out.pdb")
        atom = [[0, "O"]]
        pos = [["1.234", "5.678", "9.012"]]
        mio.Exportpospdb(atom, "1", pos, out_file)

        content = open(out_file).read()
        assert "HETATM" in content
        assert "1.234" in content

    def test_pdb_header(self, mio, tmp_path):
        out_file = str(tmp_path / "out.pdb")
        atom = [[0, "O"]]
        pos = [[0.0, 0.0, 0.0]]
        mio.Exportpospdb(atom, "1", pos, out_file)

        content = open(out_file).readlines()
        assert content[0].startswith("COMPND")
        assert "AUTHOR" in content[1]
        assert "ABMPTools" in content[1]


# ---------- convert_xyz_pdb ----------

class TestConvertXyzPdb:

    def test_creates_pdb_file(self, mio, water_xyz):
        mio.convert_xyz_pdb(water_xyz)
        pdb_path = water_xyz.replace(".xyz", ".pdb")
        assert os.path.exists(pdb_path)

    def test_pdb_content(self, mio, water_xyz):
        mio.convert_xyz_pdb(water_xyz)
        pdb_path = water_xyz.replace(".xyz", ".pdb")
        content = open(pdb_path).read()
        assert "HETATM" in content
        assert "END" in content
        hetatm_lines = [l for l in content.splitlines() if l.startswith("HETATM")]
        assert len(hetatm_lines) == 3
