# -*- coding: utf-8 -*-
import math
import pytest
import numpy as np
from abmptools.molcalc import molcalc


@pytest.fixture
def mc():
    return molcalc()


# ---------------------------------------------------------------------------
# getdist
# ---------------------------------------------------------------------------
class TestGetdist:
    def test_same_point(self, mc):
        p = np.array([1.0, 2.0, 3.0])
        assert mc.getdist(p, p) == pytest.approx(0.0)

    def test_known_distance(self, mc):
        p1 = np.array([0.0, 0.0, 0.0])
        p2 = np.array([3.0, 4.0, 0.0])
        assert mc.getdist(p1, p2) == pytest.approx(5.0)

    def test_3d_distance(self, mc):
        p1 = np.array([1.0, 2.0, 3.0])
        p2 = np.array([4.0, 6.0, 3.0])
        expected = math.sqrt(9 + 16 + 0)
        assert mc.getdist(p1, p2) == pytest.approx(expected)

    def test_negative_coords(self, mc):
        p1 = np.array([-1.0, -2.0, -3.0])
        p2 = np.array([1.0, 2.0, 3.0])
        expected = math.sqrt(4 + 16 + 36)
        assert mc.getdist(p1, p2) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# getdist_list
# ---------------------------------------------------------------------------
class TestGetdistList:
    def test_same_point(self, mc):
        p = [1.0, 2.0, 3.0]
        assert mc.getdist_list(p, p) == pytest.approx(0.0)

    def test_known_distance(self, mc):
        p1 = [0.0, 0.0, 0.0]
        p2 = [3.0, 4.0, 0.0]
        assert mc.getdist_list(p1, p2) == pytest.approx(5.0)

    def test_accepts_plain_lists(self, mc):
        p1 = [1, 2, 3]
        p2 = [4, 6, 3]
        assert mc.getdist_list(p1, p2) == pytest.approx(5.0)

    def test_consistent_with_getdist(self, mc):
        p1 = np.array([2.5, -1.3, 4.7])
        p2 = np.array([-0.5, 3.2, 1.1])
        assert mc.getdist(p1, p2) == pytest.approx(mc.getdist_list(p1, p2))


# ---------------------------------------------------------------------------
# getoriginpos
# ---------------------------------------------------------------------------
class TestGetoriginpos:
    def test_already_centered(self, mc):
        pos = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
        cell = [10.0, 10.0, 10.0]
        result = mc.getoriginpos(pos, cell)
        for atom in result:
            # centroid was (0,0,0) so nothing should change
            pass
        assert result[0][0] == pytest.approx(1.0)
        assert result[1][0] == pytest.approx(-1.0)

    def test_shift_to_origin(self, mc):
        pos = [[2.0, 2.0, 2.0], [4.0, 4.0, 4.0]]
        cell = [10.0, 10.0, 10.0]
        result = mc.getoriginpos(pos, cell)
        # centroid was (3,3,3), so shifted by -3
        assert result[0] == [pytest.approx(-1.0), pytest.approx(-1.0), pytest.approx(-1.0)]
        assert result[1] == [pytest.approx(1.0), pytest.approx(1.0), pytest.approx(1.0)]

    def test_centroid_is_origin_after(self, mc):
        pos = [[1.0, 2.0, 3.0], [5.0, 6.0, 7.0], [3.0, 4.0, 5.0]]
        cell = [20.0, 20.0, 20.0]
        result = mc.getoriginpos(pos, cell)
        cx = sum(r[0] for r in result) / len(result)
        cy = sum(r[1] for r in result) / len(result)
        cz = sum(r[2] for r in result) / len(result)
        assert cx == pytest.approx(0.0, abs=1e-12)
        assert cy == pytest.approx(0.0, abs=1e-12)
        assert cz == pytest.approx(0.0, abs=1e-12)

    def test_modifies_in_place(self, mc):
        pos = [[10.0, 20.0, 30.0], [20.0, 30.0, 40.0]]
        cell = [100.0, 100.0, 100.0]
        returned = mc.getoriginpos(pos, cell)
        assert returned is pos


# ---------------------------------------------------------------------------
# gettranspos
# ---------------------------------------------------------------------------
class TestGettranspos:
    def test_zero_translation(self, mc):
        pos = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        result = mc.gettranspos(pos, [0, 0, 0])
        assert result == pos

    def test_translation(self, mc):
        pos = [[0.0, 0.0, 0.0]]
        result = mc.gettranspos(pos, [1.0, 2.0, 3.0])
        assert result == [[1.0, 2.0, 3.0]]

    def test_multiple_atoms(self, mc):
        pos = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        vec = [10.0, 20.0, 30.0]
        result = mc.gettranspos(pos, vec)
        assert result[0] == [11.0, 22.0, 33.0]
        assert result[1] == [14.0, 25.0, 36.0]


# ---------------------------------------------------------------------------
# getmolradius
# ---------------------------------------------------------------------------
class TestGetmolradius:
    def test_single_atom_at_origin(self, mc):
        pos = [[0.0, 0.0, 0.0]]
        assert mc.getmolradius(pos) == pytest.approx(0.0)

    def test_known_radius(self, mc):
        pos = [[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]]
        assert mc.getmolradius(pos) == pytest.approx(5.0)

    def test_max_is_returned(self, mc):
        pos = [[1.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 2.0]]
        assert mc.getmolradius(pos) == pytest.approx(5.0)


# ---------------------------------------------------------------------------
# getvolume  (same logic as getmolradius)
# ---------------------------------------------------------------------------
class TestGetvolume:
    def test_matches_getmolradius(self, mc):
        pos = [[3.0, 4.0, 0.0], [1.0, 0.0, 0.0]]
        assert mc.getvolume(pos) == pytest.approx(mc.getmolradius(pos))

    def test_single_point(self, mc):
        pos = [[0.0, 0.0, 0.0]]
        assert mc.getvolume(pos) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# calccellsize
# ---------------------------------------------------------------------------
class TestCalccellsize:
    def test_known_values(self, mc):
        totalmass = 100.0  # amu
        density = 1.0      # g/cm^3
        amu_to_grams = 1.66053906660e-24
        cm_to_angstrom = 1e8
        mass_g = totalmass * amu_to_grams
        vol_cm3 = mass_g / density
        expected = math.pow(vol_cm3, 1.0 / 3.0) * cm_to_angstrom
        assert mc.calccellsize(totalmass, density) == pytest.approx(expected)

    def test_higher_density_smaller_cell(self, mc):
        cell1 = mc.calccellsize(1000.0, 1.0)
        cell2 = mc.calccellsize(1000.0, 2.0)
        assert cell2 < cell1

    def test_higher_mass_larger_cell(self, mc):
        cell1 = mc.calccellsize(500.0, 1.0)
        cell2 = mc.calccellsize(1000.0, 1.0)
        assert cell2 > cell1

    def test_water_box_approximate(self, mc):
        # 1000 water molecules: mass ~ 18.015 * 1000 amu, density ~ 1.0 g/cm^3
        totalmass = 18.015 * 1000
        density = 1.0
        cellsize = mc.calccellsize(totalmass, density)
        # Should be around 31 Angstroms
        assert 30.0 < cellsize < 33.0


# ---------------------------------------------------------------------------
# getCenter
# ---------------------------------------------------------------------------
class TestGetCenter:
    def test_single_point(self, mc):
        pos = np.array([[3.0, 4.0, 5.0]])
        center = mc.getCenter(pos)
        np.testing.assert_array_almost_equal(center, [3.0, 4.0, 5.0])

    def test_two_points(self, mc):
        pos = np.array([[0.0, 0.0, 0.0], [4.0, 6.0, 8.0]])
        center = mc.getCenter(pos)
        np.testing.assert_array_almost_equal(center, [2.0, 3.0, 4.0])

    def test_symmetric(self, mc):
        pos = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
        center = mc.getCenter(pos)
        np.testing.assert_array_almost_equal(center, [0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# moveMolTrans
# ---------------------------------------------------------------------------
class TestMoveMolTrans:
    def test_zero_translation(self, mc):
        pos = np.array([[1.0, 2.0, 3.0]])
        result = mc.moveMolTrans(pos, np.array([0.0, 0.0, 0.0]))
        np.testing.assert_array_almost_equal(result, [[1.0, 2.0, 3.0]])

    def test_translation(self, mc):
        pos = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        trans = np.array([10.0, 20.0, 30.0])
        result = mc.moveMolTrans(pos, trans)
        np.testing.assert_array_almost_equal(result, [[11.0, 22.0, 33.0], [14.0, 25.0, 36.0]])


# ---------------------------------------------------------------------------
# getcrossprod
# ---------------------------------------------------------------------------
class TestGetcrossprod:
    def test_orthogonal_unit_vectors(self, mc):
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([0.0, 1.0, 0.0])
        cross = mc.getcrossprod(p1, p2, p3)
        # (1,0,0) x (0,1,0) = (0,0,1)
        np.testing.assert_array_almost_equal(cross, [0.0, 0.0, 1.0])

    def test_parallel_vectors(self, mc):
        p1 = np.array([2.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([3.0, 0.0, 0.0])
        cross = mc.getcrossprod(p1, p2, p3)
        np.testing.assert_array_almost_equal(cross, [0.0, 0.0, 0.0])

    def test_general_case(self, mc):
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([0.0, 0.0, 1.0])
        cross = mc.getcrossprod(p1, p2, p3)
        # (1,0,0) x (0,0,1) = (0*1-0*0, 0*0-1*1, 1*0-0*0) = (0,-1,0)
        np.testing.assert_array_almost_equal(cross, [0.0, -1.0, 0.0])


# ---------------------------------------------------------------------------
# getangle
# ---------------------------------------------------------------------------
class TestGetangle:
    def test_right_angle(self, mc):
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([0.0, 1.0, 0.0])
        length = [1.0, 1.0]
        assert mc.getangle(p1, p2, p3, length) == pytest.approx(90.0)

    def test_straight_angle(self, mc):
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([-1.0, 0.0, 0.0])
        length = [1.0, 1.0]
        assert mc.getangle(p1, p2, p3, length) == pytest.approx(180.0)

    def test_acute_angle(self, mc):
        p1 = np.array([1.0, 1.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([1.0, 0.0, 0.0])
        d1 = math.sqrt(2)
        d2 = 1.0
        angle = mc.getangle(p1, p2, p3, [d1, d2])
        assert angle == pytest.approx(45.0)

    def test_zero_angle_returns_zero(self, mc):
        # When vectors are identical, acos(1.0) can raise ValueError caught by except
        p1 = np.array([1.0, 0.0, 0.0])
        p2 = np.array([0.0, 0.0, 0.0])
        p3 = np.array([2.0, 0.0, 0.0])
        length = [1.0, 2.0]
        assert mc.getangle(p1, p2, p3, length) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# getnppos
# ---------------------------------------------------------------------------
class TestGetnppos:
    def test_converts_to_numpy(self, mc):
        pos = [["1.0", "2.0", "3.0"], ["4.0", "5.0", "6.0"]]
        result = mc.getnppos(pos)
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_almost_equal(result, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    def test_with_floats(self, mc):
        pos = [[1.0, 2.0, 3.0]]
        result = mc.getnppos(pos)
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_almost_equal(result, [[1.0, 2.0, 3.0]])


# ---------------------------------------------------------------------------
# rotate_ardz
# ---------------------------------------------------------------------------
class TestRotateArdz:
    def test_zero_rotation(self, mc):
        pos = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        result = mc.rotate_ardz(0.0, pos)
        np.testing.assert_array_almost_equal(result, pos)

    def test_90_degree_rotation(self, mc):
        pos = np.array([[1.0, 0.0, 0.0]])
        result = mc.rotate_ardz(90.0, pos)
        np.testing.assert_array_almost_equal(result, [[0.0, 1.0, 0.0]], decimal=5)

    def test_180_degree_rotation(self, mc):
        pos = np.array([[1.0, 0.0, 0.0]])
        result = mc.rotate_ardz(180.0, pos)
        np.testing.assert_array_almost_equal(result, [[-1.0, 0.0, 0.0]], decimal=5)

    def test_z_component_unchanged(self, mc):
        pos = np.array([[1.0, 2.0, 5.0]])
        result = mc.rotate_ardz(45.0, pos)
        assert result[0][2] == pytest.approx(5.0)


# ---------------------------------------------------------------------------
# rotate_ardvec
# ---------------------------------------------------------------------------
class TestRotateArdvec:
    def test_rotation_around_z_axis(self, mc):
        pos = np.array([[1.0, 0.0, 0.0]])
        vec = np.array([0.0, 0.0, 1.0])
        result = mc.rotate_ardvec(90.0, vec, pos)
        np.testing.assert_array_almost_equal(result, [[0.0, 1.0, 0.0]], decimal=5)

    def test_zero_rotation(self, mc):
        pos = np.array([[1.0, 2.0, 3.0]])
        vec = np.array([0.0, 0.0, 1.0])
        result = mc.rotate_ardvec(0.0, vec, pos)
        np.testing.assert_array_almost_equal(result, pos, decimal=5)

    def test_rotation_around_x_axis(self, mc):
        pos = np.array([[0.0, 1.0, 0.0]])
        vec = np.array([1.0, 0.0, 0.0])
        result = mc.rotate_ardvec(90.0, vec, pos)
        np.testing.assert_array_almost_equal(result, [[0.0, 0.0, 1.0]], decimal=5)

    def test_360_degree_identity(self, mc):
        pos = np.array([[1.0, 2.0, 3.0]])
        vec = np.array([0.0, 0.0, 1.0])
        result = mc.rotate_ardvec(360.0, vec, pos)
        np.testing.assert_array_almost_equal(result, pos, decimal=5)


# ---------------------------------------------------------------------------
# get_element_name (static)
# ---------------------------------------------------------------------------
class TestGetElementName:
    def test_carbon(self):
        assert molcalc.get_element_name(12.011) == 'C'

    def test_hydrogen(self):
        assert molcalc.get_element_name(1.008) == 'H'

    def test_nitrogen(self):
        assert molcalc.get_element_name(14.007) == 'N'

    def test_oxygen(self):
        assert molcalc.get_element_name(15.999) == 'O'

    def test_fluorine(self):
        assert molcalc.get_element_name(18.998) == 'F'

    def test_unknown_mass(self):
        assert molcalc.get_element_name(200.0) == 'Unknown'


# ---------------------------------------------------------------------------
# get_atomic_radius (static)
# ---------------------------------------------------------------------------
class TestGetAtomicRadius:
    def test_carbon(self):
        assert molcalc.get_atomic_radius('C') == pytest.approx(2 * 1.926)

    def test_hydrogen(self):
        assert molcalc.get_atomic_radius('H') == pytest.approx(2 * 1.443)

    def test_oxygen(self):
        assert molcalc.get_atomic_radius('O') == pytest.approx(2 * 1.750)

    def test_nitrogen(self):
        assert molcalc.get_atomic_radius('N') == pytest.approx(2 * 1.830)

    def test_fluorine(self):
        assert molcalc.get_atomic_radius('F') == pytest.approx(2 * 1.682)

    def test_unknown_element(self):
        assert molcalc.get_atomic_radius('Xx') == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# wrap_to_primary_cell (static)
# ---------------------------------------------------------------------------
class TestWrapToPrimaryCell:
    def test_inside_box(self):
        box = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
        coord = [5.0, 5.0, 5.0]
        result = molcalc.wrap_to_primary_cell(coord, box)
        assert result == [pytest.approx(5.0), pytest.approx(5.0), pytest.approx(5.0)]

    def test_outside_positive(self):
        box = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
        coord = [15.0, 25.0, 35.0]
        result = molcalc.wrap_to_primary_cell(coord, box)
        assert result == [pytest.approx(5.0), pytest.approx(5.0), pytest.approx(5.0)]

    def test_outside_negative(self):
        box = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
        coord = [-3.0, -7.0, -1.0]
        result = molcalc.wrap_to_primary_cell(coord, box)
        assert result == [pytest.approx(7.0), pytest.approx(3.0), pytest.approx(9.0)]

    def test_on_boundary(self):
        box = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
        coord = [10.0, 0.0, 10.0]
        result = molcalc.wrap_to_primary_cell(coord, box)
        assert result == [pytest.approx(0.0), pytest.approx(0.0), pytest.approx(0.0)]

    def test_nonzero_lower_bound(self):
        box = (-5.0, 5.0, -5.0, 5.0, -5.0, 5.0)
        coord = [7.0, -7.0, 12.0]
        result = molcalc.wrap_to_primary_cell(coord, box)
        assert result[0] == pytest.approx(-3.0)
        assert result[1] == pytest.approx(3.0)
        assert result[2] == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# shift_molecule_to_primary_cell (static)
# ---------------------------------------------------------------------------
class TestShiftMoleculeToPrimaryCell:
    def test_molecule_inside(self):
        box = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
        coords = [[4.0, 4.0, 4.0], [6.0, 6.0, 6.0]]
        result = molcalc.shift_molecule_to_primary_cell(coords, box)
        # Centroid is (5,5,5), already inside -> no shift
        assert result[0] == [pytest.approx(4.0), pytest.approx(4.0), pytest.approx(4.0)]
        assert result[1] == [pytest.approx(6.0), pytest.approx(6.0), pytest.approx(6.0)]

    def test_molecule_outside(self):
        box = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
        coords = [[14.0, 14.0, 14.0], [16.0, 16.0, 16.0]]
        result = molcalc.shift_molecule_to_primary_cell(coords, box)
        # Centroid (15,15,15) wraps to (5,5,5), shift = -10
        assert result[0] == [pytest.approx(4.0), pytest.approx(4.0), pytest.approx(4.0)]
        assert result[1] == [pytest.approx(6.0), pytest.approx(6.0), pytest.approx(6.0)]

    def test_preserves_relative_positions(self):
        box = (0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
        coords = [[21.0, 21.0, 21.0], [23.0, 23.0, 23.0]]
        result = molcalc.shift_molecule_to_primary_cell(coords, box)
        dx = result[1][0] - result[0][0]
        dy = result[1][1] - result[0][1]
        dz = result[1][2] - result[0][2]
        assert dx == pytest.approx(2.0)
        assert dy == pytest.approx(2.0)
        assert dz == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# getindex
# ---------------------------------------------------------------------------
class TestGetindex:
    def test_basic(self, mc):
        clistall = [[0, 1, 2], [1, 0, 3], [2, 0], [3, 1]]
        result = mc.getindex(clistall)
        assert result == [0, 1, 2, 3]

    def test_with_duplicates(self, mc):
        clistall = [[0, 1, 1, 2], [1, 0, 0]]
        result = mc.getindex(clistall)
        assert result == [0, 1, 2]

    def test_empty(self, mc):
        clistall = [[], []]
        result = mc.getindex(clistall)
        assert result == []

    def test_single_entry(self, mc):
        clistall = [[5]]
        result = mc.getindex(clistall)
        assert result == [5]

    def test_sorted_output(self, mc):
        clistall = [[3, 1], [2, 0]]
        result = mc.getindex(clistall)
        assert result == [0, 1, 2, 3]


# ---------------------------------------------------------------------------
# getatomisite
# ---------------------------------------------------------------------------
class TestGetatomisite:
    def test_basic_mapping(self, mc):
        isitelist = [["C", "H", "O"], [1.7, 1.2, 1.5]]
        typenameMol = ["H", "C", "O"]
        result = mc.getatomisite(isitelist, typenameMol)
        assert result == [1.2, 1.7, 1.5]

    def test_repeated_types(self, mc):
        isitelist = [["C", "H"], [1.7, 1.2]]
        typenameMol = ["C", "C", "H"]
        result = mc.getatomisite(isitelist, typenameMol)
        assert result == [1.7, 1.7, 1.2]

    def test_empty_input(self, mc):
        isitelist = [[], []]
        typenameMol = []
        result = mc.getatomisite(isitelist, typenameMol)
        assert result == []


# ---------------------------------------------------------------------------
# getmolmass
# ---------------------------------------------------------------------------
class TestGetmolmass:
    def test_basic(self, mc):
        ffname = [["X", "C"], ["Y", "H"], ["Z", "O"]]
        atom_list = [["C", 12.011], ["H", 1.008], ["O", 15.999]]
        result = mc.getmolmass(ffname, atom_list)
        assert result == [12.011, 1.008, 15.999]

    def test_partial_match(self, mc):
        ffname = [["X", "C"], ["Y", "N"]]
        atom_list = [["C", 12.011], ["H", 1.008], ["N", 14.007]]
        result = mc.getmolmass(ffname, atom_list)
        assert result == [12.011, 14.007]

    def test_no_match(self, mc):
        ffname = [["X", "Zz"]]
        atom_list = [["C", 12.011]]
        result = mc.getmolmass(ffname, atom_list)
        assert result == []


# ---------------------------------------------------------------------------
# getrenumindex
# ---------------------------------------------------------------------------
class TestGetrenumindex:
    def test_basic_renumbering(self, mc):
        index = [2, 5, 8]
        clistall = [[2, 5], [5, 2, 8], [8, 5]]
        index_renum, clist_renum = mc.getrenumindex(index, clistall)
        assert index_renum == [0, 1, 2]
        # clistall should now reference renumbered indices
        assert clist_renum[0] == [0, 1]
        assert clist_renum[1] == [1, 0, 2]
        assert clist_renum[2] == [2, 1]

    def test_already_sequential(self, mc):
        index = [0, 1, 2]
        clistall = [[0, 1], [1, 0, 2]]
        index_renum, clist_renum = mc.getrenumindex(index, clistall)
        assert index_renum == [0, 1, 2]
        assert clist_renum[0] == [0, 1]
