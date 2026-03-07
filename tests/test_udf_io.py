import pytest
import numpy as np
from unittest.mock import MagicMock
from abmptools.udf_io import udf_io


class TestUdfIoInit:
    def test_defaults(self):
        obj = udf_io()
        assert obj.verflag is True


class TestGetposatom:
    def test_returns_numpy_array(self):
        obj = udf_io()
        mock_uobj = MagicMock()
        mock_uobj.get.return_value = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        result = obj.getposatom(mock_uobj, 0)
        np.testing.assert_array_equal(result, np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        mock_uobj.get.assert_called_once_with("Structure.Position.mol[].atom[]", [0])

    def test_different_index(self):
        obj = udf_io()
        mock_uobj = MagicMock()
        mock_uobj.get.return_value = [[7.0, 8.0, 9.0]]
        result = obj.getposatom(mock_uobj, 5)
        mock_uobj.get.assert_called_once_with("Structure.Position.mol[].atom[]", [5])


class TestGetposmol:
    def test_returns_numpy_array(self):
        obj = udf_io()
        mock_uobj = MagicMock()
        mock_uobj.get.return_value = [[1.0, 2.0, 3.0]]
        result = obj.getposmol(mock_uobj, 2)
        np.testing.assert_array_equal(result, np.array([[1.0, 2.0, 3.0]]))
        mock_uobj.get.assert_called_once_with("Structure.Position.mol[2].atom[]")


class TestGetposmolrec:
    def test_jumps_to_record(self):
        obj = udf_io()
        mock_uobj = MagicMock()
        mock_uobj.get.return_value = [[0.0, 0.0, 0.0]]
        obj.getposmolrec(mock_uobj, 1, 10)
        mock_uobj.jump.assert_called_once_with(10)


class TestGetnameAtom:
    def test_returns_atom_names(self):
        obj = udf_io()
        mock_uobj = MagicMock()
        mock_uobj.get.return_value = ["C", "H", "O"]
        result = obj.getnameAtom(mock_uobj, 0)
        assert result == ["C", "H", "O"]
        mock_uobj.get.assert_called_once_with(
            "Set_of_Molecules.molecule[].atom[].Atom_Name", [0])


class TestGetAtomtypename:
    def test_returns_type_names(self):
        obj = udf_io()
        mock_uobj = MagicMock()
        mock_uobj.get.return_value = ["CT", "HC", "OH"]
        result = obj.getAtomtypename(mock_uobj, 3)
        assert result == ["CT", "HC", "OH"]


class TestPutPositionsMol:
    def test_puts_all_atoms(self):
        obj = udf_io()
        mock_uobj = MagicMock()
        mock_uobj.size.return_value = 2
        positions = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        obj.putPositionsMol(mock_uobj, 0, positions)
        assert mock_uobj.put.call_count == 6  # 2 atoms * 3 coords
        mock_uobj.size.assert_called_once_with("Structure.Position.mol[].atom[]", [0])
