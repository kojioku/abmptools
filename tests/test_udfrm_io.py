"""Tests for abmptools.udfrm_io module."""

import pytest
from unittest.mock import MagicMock, patch, call
import numpy as np

from abmptools.udfrm_io import udfrm_io


class TestUdfrmIoInit:
    """Tests for udfrm_io.__init__."""

    def test_molflag_default(self):
        obj = udfrm_io()
        assert obj.molflag is False

    def test_cell_default(self):
        obj = udfrm_io()
        assert obj.cell is None


class TestGetmolname:
    """Tests for udfrm_io.getmolname."""

    def test_getmolname_returns_mol_name(self):
        obj = udfrm_io()
        uobj = MagicMock()
        uobj.get.return_value = "Water"

        result = obj.getmolname(0, uobj)

        assert result == "Water"
        uobj.get.assert_called_once_with("Set_of_Molecules.molecule[0].Mol_Name")

    def test_getmolname_with_different_index(self):
        obj = udfrm_io()
        uobj = MagicMock()
        uobj.get.return_value = "Ethanol"

        result = obj.getmolname(5, uobj)

        assert result == "Ethanol"
        uobj.get.assert_called_once_with("Set_of_Molecules.molecule[5].Mol_Name")


class TestMoveintocellRec:
    """Tests for udfrm_io.moveintocell_rec."""

    def test_calls_jump_and_get_cell(self):
        obj = udfrm_io()
        uobj = MagicMock()
        uobj.get.return_value = [10.0, 10.0, 10.0]

        # totalMol=0 means the loop body never runs, so we just test
        # that jump and get are called correctly.
        obj.moveintocell_rec(uobj, 3, 0)

        uobj.jump.assert_called_once_with(3)
        uobj.get.assert_called_once_with("Structure.Unit_Cell.Cell_Size")

    def test_processes_molecules_in_range(self):
        """Test that molecules are iterated and positions are processed."""
        obj = udfrm_io()
        uobj = MagicMock()
        uobj.get.return_value = [10.0, 10.0, 10.0]

        # Mock inherited methods
        obj.getposmol = MagicMock(return_value=np.array([[5.0, 5.0, 5.0]]))
        obj.getCenter = MagicMock(return_value=[5.0, 5.0, 5.0])
        obj.moveMolTrans = MagicMock(return_value=np.array([[5.0, 5.0, 5.0]]))
        obj.putPositionsMol = MagicMock()

        obj.moveintocell_rec(uobj, 0, 2)

        assert obj.getposmol.call_count == 2
        assert obj.getCenter.call_count == 2
        assert obj.moveMolTrans.call_count == 2
        assert obj.putPositionsMol.call_count == 2

    def test_moves_molecule_when_outside_positive(self):
        """Test molecule center > cell is shifted back into cell."""
        obj = udfrm_io()
        uobj = MagicMock()
        uobj.get.return_value = [10.0, 10.0, 10.0]

        # Molecule center at (15, 5, 5) -- x is outside cell
        obj.getposmol = MagicMock(return_value=np.array([[15.0, 5.0, 5.0]]))
        obj.getCenter = MagicMock(return_value=[15.0, 5.0, 5.0])
        obj.moveMolTrans = MagicMock(return_value=np.array([[5.0, 5.0, 5.0]]))
        obj.putPositionsMol = MagicMock()

        obj.moveintocell_rec(uobj, 0, 1)

        # transVec should be [-10, 0, 0] for x-direction
        args = obj.moveMolTrans.call_args
        trans_vec = args[0][1]
        assert trans_vec[0] == pytest.approx(-10.0)
        assert trans_vec[1] == pytest.approx(0.0)
        assert trans_vec[2] == pytest.approx(0.0)

    def test_moves_molecule_when_outside_negative(self):
        """Test molecule center < 0 is shifted back into cell."""
        obj = udfrm_io()
        uobj = MagicMock()
        uobj.get.return_value = [10.0, 10.0, 10.0]

        # Molecule center at (-5, 5, 5) -- x is negative
        obj.getposmol = MagicMock(return_value=np.array([[-5.0, 5.0, 5.0]]))
        obj.getCenter = MagicMock(return_value=[-5.0, 5.0, 5.0])
        obj.moveMolTrans = MagicMock(return_value=np.array([[5.0, 5.0, 5.0]]))
        obj.putPositionsMol = MagicMock()

        obj.moveintocell_rec(uobj, 0, 1)

        args = obj.moveMolTrans.call_args
        trans_vec = args[0][1]
        assert trans_vec[0] == pytest.approx(10.0)
        assert trans_vec[1] == pytest.approx(0.0)
        assert trans_vec[2] == pytest.approx(0.0)


class TestConvertUdfPdb:
    """Tests for udfrm_io.convert_udf_pdb."""

    def _setup_obj_and_uobj(self, molflag=False):
        obj = udfrm_io()
        obj.molflag = molflag
        uobj = MagicMock()
        uobj.get.return_value = [10.0, 10.0, 10.0]

        # Mock inherited methods
        obj.getposmol = MagicMock(return_value=np.array([[1.0, 2.0, 3.0]]))
        obj.getAtomtypename = MagicMock(return_value=["C", "H"])
        obj.getmolname = MagicMock(return_value="MOL")
        obj.Exportpos = MagicMock()
        obj.Exporttgtmolpos = MagicMock()

        return obj, uobj

    def test_calls_jump_and_gets_cell(self):
        obj, uobj = self._setup_obj_and_uobj()
        obj.convert_udf_pdb(0, uobj, 1, "test", writef=False)

        uobj.jump.assert_called_once_with(0)
        uobj.get.assert_called_once_with("Structure.Unit_Cell.Cell_Size")

    def test_no_molflag_iterates_totalMol(self):
        obj, uobj = self._setup_obj_and_uobj(molflag=False)
        result = obj.convert_udf_pdb(0, uobj, 3, "test", writef=False)

        assert obj.getposmol.call_count == 3
        assert obj.getAtomtypename.call_count == 3
        assert obj.getmolname.call_count == 3
        # Returns [typenameMol, posMol, molnamelist]
        assert len(result) == 3
        assert len(result[0]) == 3  # typenameMol has 3 entries
        assert len(result[1]) == 3  # posMol has 3 entries
        assert len(result[2]) == 3  # molnamelist has 3 entries

    def test_molflag_iterates_single_mol(self):
        obj, uobj = self._setup_obj_and_uobj(molflag=True)
        result = obj.convert_udf_pdb(0, uobj, 2, "test", writef=False)

        # With molflag=True, only one molecule (tgtmol=2) is processed
        obj.getposmol.assert_called_once()
        obj.getAtomtypename.assert_called_once()
        obj.getmolname.assert_called_once_with(2, uobj)
        assert len(result[0]) == 1

    def test_no_molflag_calls_exportpos_when_writef(self):
        obj, uobj = self._setup_obj_and_uobj(molflag=False)
        obj.convert_udf_pdb(0, uobj, 2, "outname", writef=True)

        obj.Exportpos.assert_called_once_with(".", 0, 2, uobj, "outname.xyz")

    def test_molflag_calls_exporttgtmolpos_when_writef(self):
        obj, uobj = self._setup_obj_and_uobj(molflag=True)
        obj.convert_udf_pdb(0, uobj, 2, "outname", writef=True)

        obj.Exporttgtmolpos.assert_called_once_with(".", "outname.xyz", 0, [2], uobj)
