import pytest
import numpy as np

from abmptools.readcif import getcartesiancellvec, getcartesianmol, intocell, getsymcoord


class TestGetCartesianCellVec:
    """Tests for getcartesiancellvec(angle, length)."""

    def test_cubic_cell(self, capsys):
        angle = [90, 90, 90]
        length = [10, 10, 10]
        a_vec, b_vec, c_vec = getcartesiancellvec(angle, length)
        # Suppress print output check; just verify vectors
        assert a_vec == pytest.approx([10, 0, 0])
        assert b_vec == pytest.approx([0, 10, 0], abs=1e-10)
        assert c_vec == pytest.approx([0, 0, 10], abs=1e-10)

    def test_orthorhombic_cell(self, capsys):
        angle = [90, 90, 90]
        length = [5, 8, 12]
        a_vec, b_vec, c_vec = getcartesiancellvec(angle, length)
        assert a_vec == pytest.approx([5, 0, 0])
        assert b_vec == pytest.approx([0, 8, 0], abs=1e-10)
        assert c_vec == pytest.approx([0, 0, 12], abs=1e-10)

    def test_vector_lengths_match_input(self, capsys):
        angle = [90, 90, 90]
        length = [3.5, 4.2, 6.1]
        a_vec, b_vec, c_vec = getcartesiancellvec(angle, length)
        assert np.linalg.norm(a_vec) == pytest.approx(3.5)
        assert np.linalg.norm(b_vec) == pytest.approx(4.2, abs=1e-10)
        assert np.linalg.norm(c_vec) == pytest.approx(6.1, abs=1e-10)

    def test_nonorthogonal_vector_lengths(self, capsys):
        """For any cell, the magnitudes of a_vec, b_vec, c_vec must equal a, b, c."""
        angle = [90, 90, 120]
        length = [5, 5, 8]
        a_vec, b_vec, c_vec = getcartesiancellvec(angle, length)
        assert np.linalg.norm(a_vec) == pytest.approx(5.0)
        assert np.linalg.norm(b_vec) == pytest.approx(5.0, abs=1e-10)
        assert np.linalg.norm(c_vec) == pytest.approx(8.0, abs=1e-10)


class TestGetCartesianMol:
    """Tests for getcartesianmol(a, b, c, acell, bcell, ccell)."""

    def test_simple_scaling(self):
        acell = [10, 0, 0]
        bcell = [0, 10, 0]
        ccell = [0, 0, 10]
        a_out, b_out, c_out = getcartesianmol(0.5, 0.5, 0.5, acell, bcell, ccell)
        assert a_out == pytest.approx([5, 0, 0])
        assert b_out == pytest.approx([0, 5, 0])
        assert c_out == pytest.approx([0, 0, 5])

    def test_zero_fractional(self):
        acell = [10, 0, 0]
        bcell = [0, 10, 0]
        ccell = [0, 0, 10]
        a_out, b_out, c_out = getcartesianmol(0.0, 0.0, 0.0, acell, bcell, ccell)
        assert a_out == pytest.approx([0, 0, 0])
        assert b_out == pytest.approx([0, 0, 0])
        assert c_out == pytest.approx([0, 0, 0])

    def test_unit_fractional(self):
        acell = [7, 0, 0]
        bcell = [0, 8, 0]
        ccell = [0, 0, 9]
        a_out, b_out, c_out = getcartesianmol(1.0, 1.0, 1.0, acell, bcell, ccell)
        assert a_out == pytest.approx([7, 0, 0])
        assert b_out == pytest.approx([0, 8, 0])
        assert c_out == pytest.approx([0, 0, 9])


class TestIntocell:
    """Tests for intocell(incoords, anum_mol)."""

    def test_already_in_range_single_mol(self):
        incoords = [0.1, 0.2, 0.3]
        result = intocell(incoords, [3])
        assert result == pytest.approx([0.1, 0.2, 0.3])

    def test_already_in_range_two_mol(self):
        incoords = [0.1, 0.2, 0.3, 0.4, 0.5]
        result = intocell(incoords, [2, 3])
        assert result == pytest.approx([0.1, 0.2, 0.3, 0.4, 0.5])

    def test_shift_positive_single_mol(self):
        """Coords with mean > 1 should be shifted down by 1."""
        incoords = [1.1, 1.2, 1.3]
        result = intocell(incoords, [3])
        assert result == pytest.approx([0.1, 0.2, 0.3])

    def test_shift_negative_single_mol(self):
        """Coords with mean < 0 should be shifted up by 1."""
        incoords = [-0.9, -0.8, -0.7]
        result = intocell(incoords, [3])
        assert result == pytest.approx([0.1, 0.2, 0.3])

    def test_shift_two_molecules(self):
        """Each molecule group is shifted independently."""
        # mol1 (2 atoms): mean = 1.5 -> shift -1 -> mean = 0.5
        # mol2 (2 atoms): mean = -0.5 -> shift +1 -> mean = 0.5
        incoords = [1.3, 1.7, -0.3, -0.7]
        result = intocell(incoords, [2, 2])
        assert result == pytest.approx([0.3, 0.7, 0.7, 0.3])

    def test_large_positive_shift(self):
        """Coords far above 1 need multiple shifts."""
        incoords = [2.5, 2.5, 2.5]
        result = intocell(incoords, [3])
        # mean 2.5 -> 1.5 -> 0.5
        assert result == pytest.approx([0.5, 0.5, 0.5])

    def test_large_negative_shift(self):
        """Coords far below 0 need multiple shifts."""
        incoords = [-1.5, -1.5, -1.5]
        result = intocell(incoords, [3])
        # mean -1.5 -> -0.5 -> 0.5
        assert result == pytest.approx([0.5, 0.5, 0.5])


class TestGetSymCoord:
    """Tests for getsymcoord(sympos, incoord, label)."""

    def setup_method(self):
        self.coord = [0.1, 0.2, 0.3]

    def test_x(self):
        assert getsymcoord('x', self.coord, 0) == pytest.approx(0.1)

    def test_y(self):
        assert getsymcoord('y', self.coord, 1) == pytest.approx(0.2)

    def test_z(self):
        assert getsymcoord('z', self.coord, 2) == pytest.approx(0.3)

    def test_neg_x(self):
        assert getsymcoord('-x', self.coord, 0) == pytest.approx(-0.1)

    def test_neg_y(self):
        assert getsymcoord('-y', self.coord, 1) == pytest.approx(-0.2)

    def test_neg_z(self):
        assert getsymcoord('-z', self.coord, 2) == pytest.approx(-0.3)

    def test_half_plus_x(self):
        assert getsymcoord('1/2+x', self.coord, 0) == pytest.approx(0.6)

    def test_half_minus_x(self):
        assert getsymcoord('1/2-x', self.coord, 0) == pytest.approx(0.4)

    def test_half_plus_y(self):
        assert getsymcoord('1/2+y', self.coord, 1) == pytest.approx(0.7)

    def test_half_minus_z(self):
        assert getsymcoord('1/2-z', self.coord, 2) == pytest.approx(0.2)

    def test_x_minus_y(self):
        assert getsymcoord('x-y', self.coord, 0) == pytest.approx(-0.1)

    def test_neg_x_plus_y(self):
        assert getsymcoord('-x+y', self.coord, 0) == pytest.approx(0.1)

    def test_one_third_plus_x(self):
        assert getsymcoord('1/3+x', self.coord, 0) == pytest.approx(1/3 + 0.1)

    def test_two_thirds_minus_y(self):
        assert getsymcoord('2/3-y', self.coord, 1) == pytest.approx(2/3 - 0.2)

    def test_one_quarter_plus_z(self):
        assert getsymcoord('1/4+z', self.coord, 2) == pytest.approx(0.25 + 0.3)

    def test_three_quarters_minus_z(self):
        assert getsymcoord('3/4-z', self.coord, 2) == pytest.approx(0.75 - 0.3)

    def test_neg_half_plus_x(self):
        assert getsymcoord('-1/2+x', self.coord, 0) == pytest.approx(-0.5 + 0.1)

    def test_neg_half_minus_y(self):
        assert getsymcoord('-1/2-y', self.coord, 1) == pytest.approx(-0.5 - 0.2)

    def test_one_sixth_plus_z(self):
        assert getsymcoord('1/6+z', self.coord, 2) == pytest.approx(1/6 + 0.3)

    def test_five_sixths_minus_x(self):
        assert getsymcoord('5/6-x', self.coord, 0) == pytest.approx(5/6 - 0.1)

    def test_one_third_plus_x_minus_y(self):
        assert getsymcoord('1/3+x-y', self.coord, 0) == pytest.approx(1/3 + 0.1 - 0.2)

    def test_two_thirds_minus_x_plus_y(self):
        assert getsymcoord('2/3-x+y', self.coord, 0) == pytest.approx(2/3 - 0.1 + 0.2)

    def test_two_thirds_plus_x_minus_y(self):
        assert getsymcoord('2/3+x-y', self.coord, 0) == pytest.approx(2/3 + 0.1 - 0.2)

    def test_alternate_form_x_plus_half(self):
        assert getsymcoord('x+1/2', self.coord, 0) == pytest.approx(0.6)

    def test_alternate_form_neg_z_plus_half(self):
        assert getsymcoord('-z+1/2', self.coord, 2) == pytest.approx(0.5 - 0.3)
