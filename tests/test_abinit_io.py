import pytest
from abmptools.abinit_io import abinit_io


@pytest.fixture
def aio():
    """Create a fresh abinit_io instance for each test."""
    return abinit_io()


# --------------------------------------------------------------------------- #
#  __init__ default values
# --------------------------------------------------------------------------- #
class TestInit:

    def test_ajf_method_default(self, aio):
        assert aio.ajf_method == "MP2"

    def test_ajf_basis_set_default(self, aio):
        assert aio.ajf_basis_set == "6-31G*"

    def test_piedaflag_default(self, aio):
        assert aio.piedaflag is True

    def test_cpfflag_default(self, aio):
        assert aio.cpfflag is False

    def test_cpfver_default(self, aio):
        assert aio.cpfver == 23

    def test_cmmflag_default(self, aio):
        assert aio.cmmflag is False

    def test_nofzc_flag_default(self, aio):
        assert aio.nofzc_flag is False

    def test_bsseflag_default(self, aio):
        assert aio.bsseflag is False

    def test_npro_default(self, aio):
        assert aio.npro == 4

    def test_memory_default(self, aio):
        assert aio.memory == "1800"

    def test_solv_flag_default(self, aio):
        assert aio.solv_flag is False

    def test_autofrag_default(self, aio):
        assert aio.autofrag is False

    def test_nbo_default(self, aio):
        assert aio.nbo is True

    def test_resp_default(self, aio):
        assert aio.resp is True

    def test_abinit_ver_default(self, aio):
        assert aio.abinit_ver == "rev22"

    def test_fatomnums_default_empty(self, aio):
        assert aio.fatomnums == []

    def test_fchgs_default_empty(self, aio):
        assert aio.fchgs == []

    def test_fbaas_default_empty(self, aio):
        assert aio.fbaas == []

    def test_fatminfos_default_empty(self, aio):
        assert aio.fatminfos == []

    def test_connects_default_empty(self, aio):
        assert aio.connects == []

    def test_disp_default(self, aio):
        assert aio.disp is False

    def test_is_xyz_default(self, aio):
        assert aio.is_xyz is False

    def test_natom_default(self, aio):
        assert aio.natom == 0

    def test_dgemm_default(self, aio):
        assert aio.dgemm is False

    def test_readgeom_default(self, aio):
        assert aio.readgeom == ""

    def test_writegeom_default(self, aio):
        assert aio.writegeom == ""


# --------------------------------------------------------------------------- #
#  chkdepth
# --------------------------------------------------------------------------- #
class TestChkdepth:

    def test_empty_list_returns_zero(self, aio):
        assert aio.chkdepth([]) == 0

    def test_flat_list_returns_one(self, aio):
        assert aio.chkdepth([1, 2, 3]) == 1

    def test_nested_two_levels(self, aio):
        assert aio.chkdepth([[1, 2], [3, 4]]) == 2

    def test_nested_three_levels(self, aio):
        assert aio.chkdepth([[[1]]]) == 3

    def test_non_list_int_returns_zero(self, aio):
        assert aio.chkdepth(5) == 0

    def test_non_list_string_returns_zero(self, aio):
        assert aio.chkdepth("string") == 0

    def test_none_returns_zero(self, aio):
        assert aio.chkdepth(None) == 0

    def test_mixed_depth_returns_max(self, aio):
        # max depth path is [1, [2, [3]]] -> depth 3
        assert aio.chkdepth([1, [2, [3]]]) == 3

    def test_single_element_flat(self, aio):
        assert aio.chkdepth([42]) == 1

    def test_list_of_empty_lists(self, aio):
        # [[]] -> 1 + max(chkdepth([])) = 1 + 0 = 1
        assert aio.chkdepth([[]]) == 1

    def test_deeply_nested(self, aio):
        assert aio.chkdepth([[[[1]]]]) == 4


# --------------------------------------------------------------------------- #
#  get_fragsection
# --------------------------------------------------------------------------- #
class TestGetFragsection:

    def test_depth1_wraps_in_list(self, aio):
        """When fatomnums is a flat list (depth 1), get_fragsection wraps
        each attribute in another list before processing."""
        aio.fatomnums = [3, 2]
        aio.fchgs = [0, 0]
        aio.fbaas = [0, 0]
        aio.connects = []
        aio.fatminfos = [[1, 2, 3], [4, 5]]

        result = aio.get_fragsection()
        assert isinstance(result, str)
        # frag_atom section: two numbers formatted as %8d
        assert "       3" in result
        assert "       2" in result

    def test_depth2_uses_as_is(self, aio):
        """When fatomnums is a nested list (depth > 1), attributes are used
        directly without wrapping."""
        aio.fatomnums = [[3, 2]]
        aio.fchgs = [[0, 0]]
        aio.fbaas = [[0, 0]]
        aio.connects = [[]]
        aio.fatminfos = [[[1, 2, 3], [4, 5]]]

        result = aio.get_fragsection()
        assert isinstance(result, str)
        assert "       3" in result
        assert "       2" in result

    def test_depth1_and_depth2_produce_same_output(self, aio):
        """Depth-1 data wrapped should produce the same output as the
        equivalent depth-2 data used directly."""
        # Depth 1 instance
        aio1 = abinit_io()
        aio1.fatomnums = [3, 2]
        aio1.fchgs = [0, 0]
        aio1.fbaas = [0, 0]
        aio1.connects = []
        aio1.fatminfos = [[1, 2, 3], [4, 5]]

        # Depth 2 instance (already wrapped)
        aio2 = abinit_io()
        aio2.fatomnums = [[3, 2]]
        aio2.fchgs = [[0, 0]]
        aio2.fbaas = [[0, 0]]
        aio2.connects = [[]]
        aio2.fatminfos = [[[1, 2, 3], [4, 5]]]

        assert aio1.get_fragsection() == aio2.get_fragsection()

    def test_output_contains_formatted_atom_numbers(self, aio):
        """Verify that fragment atom numbers appear with 8-char width
        formatting."""
        aio.fatomnums = [5]
        aio.fchgs = [0]
        aio.fbaas = [1]
        aio.connects = [[0, 1]]
        aio.fatminfos = [[1, 2, 3, 4, 5]]

        result = aio.get_fragsection()
        # %8d formatting: atom index 1 appears as "       1"
        assert "       1" in result

    def test_single_fragment_no_connects(self, aio):
        """A single fragment with no inter-fragment connections."""
        aio.fatomnums = [4]
        aio.fchgs = [0]
        aio.fbaas = [0]
        aio.connects = []
        aio.fatminfos = [[1, 2, 3, 4]]

        result = aio.get_fragsection()
        assert isinstance(result, str)
        # Should contain atom count, charge, baa, and fragment body lines
        lines = [l for l in result.split("\n") if l.strip()]
        assert len(lines) >= 1

    def test_multiple_fragments_depth1(self, aio):
        """Multiple fragments at depth 1 with connect info."""
        aio.fatomnums = [3, 2]
        aio.fchgs = [0, -1]
        aio.fbaas = [1, 1]
        aio.connects = [[2, 3]]
        aio.fatminfos = [[1, 2, 3], [4, 5]]

        result = aio.get_fragsection()
        assert isinstance(result, str)
        # Charges: 0 and -1
        assert "       0" in result
        assert "      -1" in result

    def test_connect_with_three_elements(self, aio):
        """Connections with a third element (e.g., bond order) are included."""
        aio.fatomnums = [2, 2]
        aio.fchgs = [0, 0]
        aio.fbaas = [1, 1]
        aio.connects = [[1, 2, 3]]
        aio.fatminfos = [[1, 2], [3, 4]]

        result = aio.get_fragsection()
        # The third element (3) should appear in the connect section
        assert "       3" in result
