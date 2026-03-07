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


# --------------------------------------------------------------------------- #
#  flatten
# --------------------------------------------------------------------------- #
class TestFlatten:

    def test_flatten_nested_lists(self, aio):
        assert aio.flatten([['1', '2'], ['3', '4']]) == [1, 2, 3, 4]

    def test_flatten_single_inner(self, aio):
        assert aio.flatten([['10']]) == [10]

    def test_flatten_empty(self, aio):
        assert aio.flatten([]) == []


# --------------------------------------------------------------------------- #
#  functor
# --------------------------------------------------------------------------- #
class TestFunctor:

    def test_functor_int_on_nested(self, aio):
        result = aio.functor(int, [['1', '2'], ['3']])
        assert result == [[1, 2], [3]]

    def test_functor_on_scalar(self, aio):
        assert aio.functor(str, 42) == '42'

    def test_functor_deeply_nested(self, aio):
        result = aio.functor(int, [[['5']]])
        assert result == [[[5]]]


# --------------------------------------------------------------------------- #
#  getfraginfo
# --------------------------------------------------------------------------- #
class TestGetfraginfo:

    def test_basic_fraginfo(self, aio):
        seg1 = {'seg_info': [[1, 2, 3], [4, 5]]}
        seg2 = {'seg_info': [[6, 7]]}
        result = aio.getfraginfo(seg1, seg2)
        assert result == [[1, 2], [3]]

    def test_single_frag_each(self, aio):
        seg1 = {'seg_info': [[1]]}
        seg2 = {'seg_info': [[2]]}
        result = aio.getfraginfo(seg1, seg2)
        assert result == [[1], [2]]


# --------------------------------------------------------------------------- #
#  config_read — uses tmp_path to create a temp segment_data.dat
# --------------------------------------------------------------------------- #
class TestConfigRead:

    def test_config_read_finds_matching_segment(self, aio, tmp_path):
        dat_content = """seg_data = [
    {
        'name': 'water',
        'atom': [3],
        'charge': [0],
        'connect_num': [0],
        'connect': [],
        'seg_info': [[1, 2, 3]],
        'nummol_seg': [1],
        'repeat': [1],
        'pair_file': [],
        'multi_xyz': 'none'
    },
]
"""
        (tmp_path / "segment_data.dat").write_text(dat_content)
        aio.mainpath = str(tmp_path)
        result = aio.config_read('water')
        assert result['name'] == 'water'
        assert result['atom'] == [3]
        assert result['charge'] == [0]
        assert result['seg_info'] == [[1, 2, 3]]

    def test_config_read_returns_default_when_not_found(self, aio, tmp_path):
        dat_content = "seg_data = []\n"
        (tmp_path / "segment_data.dat").write_text(dat_content)
        aio.mainpath = str(tmp_path)
        result = aio.config_read('unknown', seg_atom=5)
        assert result['name'] == 'unknown'
        assert result['atom'] == [5]
        assert result['seg_info'] == [[1, 2, 3, 4, 5]]
        assert result['nummol_seg'] == [1]

    def test_config_read_overrides_defaults_with_seg_data(self, aio, tmp_path):
        dat_content = """seg_data = [
    {
        'name': 'methane',
        'atom': [5],
        'charge': [-1],
    },
]
"""
        (tmp_path / "segment_data.dat").write_text(dat_content)
        aio.mainpath = str(tmp_path)
        result = aio.config_read('methane', seg_atom=0)
        assert result['charge'] == [-1]
        # Keys not in seg_data retain defaults
        assert result['connect'] == []
        assert result['nterm'] == 2


# --------------------------------------------------------------------------- #
#  gen_ajf_body — generates AJF file content from internal state
# --------------------------------------------------------------------------- #
class TestGenAjfBody:

    def _make_aio(self):
        obj = abinit_io()
        obj.para_job = 1
        return obj

    def test_contains_method(self):
        aio = self._make_aio()
        aio.ajf_method = "HF"
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "Method='HF'" in body

    def test_contains_basis_set(self):
        aio = self._make_aio()
        aio.ajf_basis_set = "cc-pVDZ"
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "BasisSet='cc-pVDZ'" in body

    def test_single_fragment_uses_fmo_off(self):
        aio = self._make_aio()
        param = [0, "", 1, 0]
        body = aio.gen_ajf_body(param)
        assert "FMO='OFF'" in body

    def test_multiple_fragments_uses_fmo_on(self):
        aio = self._make_aio()
        param = [0, "frag_data_here", 3, 0]
        body = aio.gen_ajf_body(param)
        assert "FMO='ON'" in body
        assert "NF= 3" in body

    def test_charge_appears_in_cntrl(self):
        aio = self._make_aio()
        param = [-2, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "Charge=-2" in body

    def test_pieda_yes_when_flag_true(self):
        aio = self._make_aio()
        aio.piedaflag = True
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "PIEDA='YES'" in body

    def test_pieda_no_when_flag_false(self):
        aio = self._make_aio()
        aio.piedaflag = False
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "PIEDA='NO'" in body

    def test_bsse_on_when_flag_true(self):
        aio = self._make_aio()
        aio.bsseflag = True
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "CP='ON'" in body

    def test_bsse_absent_when_flag_false(self):
        aio = self._make_aio()
        aio.bsseflag = False
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "CP='ON'" not in body

    def test_dgemm_sections_present(self):
        aio = self._make_aio()
        aio.dgemm = True
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "MOD1ST='GEMM'" in body
        assert "MOD4TH='GEMM'" in body

    def test_dgemm_absent_by_default(self):
        aio = self._make_aio()
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "MOD1ST='GEMM'" not in body

    def test_solv_on(self):
        aio = self._make_aio()
        aio.solv_flag = True
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "EFFECT='ON'" in body

    def test_solv_off(self):
        aio = self._make_aio()
        aio.solv_flag = False
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "EFFECT='OFF'" in body

    def test_disp_on_in_lrd(self):
        aio = self._make_aio()
        aio.disp = True
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "DISP='ON'" in body

    def test_autofrag_on(self):
        aio = self._make_aio()
        aio.autofrag = True
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "AutoFrag='ON'" in body

    def test_autofrag_off(self):
        aio = self._make_aio()
        aio.autofrag = False
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "AutoFrag='OFF'" in body

    def test_fragment_section_included(self):
        aio = self._make_aio()
        param = [0, "MY_FRAG_DATA", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "&FRAGMENT" in body
        assert "MY_FRAG_DATA" in body

    def test_resp_and_nbo(self):
        aio = self._make_aio()
        aio.resp = True
        aio.nbo = True
        param = [0, "", 2, 0]
        body = aio.gen_ajf_body(param)
        assert "ESPFIT='ON'" in body
        assert "NBOANL='ON'" in body

    def test_mome_flag(self):
        aio = self._make_aio()
        param = [0, "", 3, 1]
        body = aio.gen_ajf_body(param)
        assert "LMOTYP='MOME'" in body


# --------------------------------------------------------------------------- #
#  read_ifie — parses IFIE data from log file
# --------------------------------------------------------------------------- #
class TestReadIfie:

    IFIE_LOG = """\
 some header line
 ## MP2-IFIE

           IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE
                      / A      APPROX.   / Hartree  / Hartree
        ----------------------------------------------------------------------------------------------------
    2    1    0.001000   F      -0.050000  -0.020000  -0.030000  -0.010000  -0.005000  -0.002000
    3    1    0.002000   F      -0.100000  -0.040000  -0.060000  -0.020000  -0.010000  -0.004000
 ## Mulliken
 rest of file
"""

    def test_read_ifie_basic(self, aio, tmp_path):
        f = tmp_path / "test.out"
        f.write_text(self.IFIE_LOG)
        result = aio.read_ifie(str(f))
        assert len(result) == 2
        # columns: [0]=i, [1]=j, [2]=dist, [3]=F, [4]=dimer_es, [5]=hf_ifie, [6]=mp2_ifie
        assert result[0][0] == '2'
        assert result[0][1] == '1'
        assert float(result[0][4]) == pytest.approx(-0.05)
        assert float(result[1][4]) == pytest.approx(-0.1)

    def test_read_ifie_missing_file(self, aio, tmp_path):
        result = aio.read_ifie(str(tmp_path / "nonexistent.out"))
        assert result == []

    def test_read_ifie_no_ifie_section(self, aio, tmp_path):
        f = tmp_path / "empty.out"
        f.write_text("no ifie data here\njust some text\n")
        result = aio.read_ifie(str(f))
        assert result == []

    def test_read_ifie_clamps_large_negative(self, aio, tmp_path):
        """Values < -2 in columns 4 or 5 are zeroed out."""
        content = """\
 header
 ## MP2-IFIE

           IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE
                      / A      APPROX.   / Hartree  / Hartree
        ----------------------------------------------------------------------------------------------------
    2    1    0.001000   F      -3.000000  -0.020000  -0.030000  -0.010000  -0.005000  -0.002000
 ## Mulliken
"""
        f = tmp_path / "clamp.out"
        f.write_text(content)
        result = aio.read_ifie(str(f))
        assert len(result) == 1
        # columns [4]=dimer_es, [5]=hf_ifie clamped to 0.0 because [4] < -2
        assert result[0][4] == 0.0
        assert result[0][5] == 0.0

    def test_read_ifie_stops_at_pieda(self, aio, tmp_path):
        content = """\
 header
 ## MP2-IFIE

           IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE
                      / A      APPROX.   / Hartree  / Hartree
        ----------------------------------------------------------------------------------------------------
    2    1    0.001000   F      -0.050000  -0.020000  -0.030000  -0.010000  -0.005000  -0.002000
 ## PIEDA
    2    1    0.001000   F      -0.999000  -0.999000  -0.999000  -0.999000  -0.999000  -0.999000
 ## Mulliken
"""
        f = tmp_path / "pieda.out"
        f.write_text(content)
        result = aio.read_ifie(str(f))
        assert len(result) == 1
        assert float(result[0][4]) == pytest.approx(-0.05)


# --------------------------------------------------------------------------- #
#  getmo_or_fmo — detects FMO vs MO calculation
# --------------------------------------------------------------------------- #
class TestGetmoOrFmo:

    def test_fmo_detected(self, aio, tmp_path):
        f = tmp_path / "fmo.out"
        f.write_text("some line\n FMO ON\nanother line\n")
        assert aio.getmo_or_fmo(str(f)) is True

    def test_mo_detected(self, aio, tmp_path):
        f = tmp_path / "mo.out"
        f.write_text("some line\nno fmo here\nanother line\n")
        assert aio.getmo_or_fmo(str(f)) is False


# --------------------------------------------------------------------------- #
#  getifiesum — pure computation on energy list
# --------------------------------------------------------------------------- #
class TestGetifiesum:

    def test_empty_energy_returns_zeros(self, aio):
        result = aio.getifiesum([], [[1], [2]])
        assert result == [0, 0]
