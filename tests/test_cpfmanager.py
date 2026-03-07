"""Tests for abmptools.cpfmanager.CPFManager."""

import os
import pytest

from abmptools.cpfmanager import CPFManager


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_CPF_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir,
    "sample", "convertcpf", "gly5",
)

SAMPLE_CPF_23 = os.path.join(SAMPLE_CPF_DIR, "gly5-23.cpf")
SAMPLE_CPF_10 = os.path.join(SAMPLE_CPF_DIR, "gly5-10.cpf")


@pytest.fixture
def cpfm():
    """Return a fresh CPFManager instance."""
    return CPFManager()


# ---------------------------------------------------------------------------
# __init__ defaults
# ---------------------------------------------------------------------------

class TestInit:
    def test_default_cpfver(self, cpfm):
        assert cpfm.cpfver == 23

    def test_default_is_gz(self, cpfm):
        assert cpfm.is_gz is False

    def test_default_is_difie(self, cpfm):
        assert cpfm.is_difie is False

    def test_default_tgtfrag(self, cpfm):
        assert cpfm.tgtfrag == 0

    def test_default_atominfo(self, cpfm):
        assert cpfm.atominfo == {}

    def test_default_fraginfo(self, cpfm):
        assert cpfm.fraginfo == []

    def test_default_condition(self, cpfm):
        assert cpfm.condition == []

    def test_default_static_data(self, cpfm):
        assert cpfm.static_data == {}

    def test_default_mominfo(self, cpfm):
        assert cpfm.mominfo == {}

    def test_default_diminfo(self, cpfm):
        assert cpfm.diminfo == {}

    def test_default_labels(self, cpfm):
        assert cpfm.labels == {}


# ---------------------------------------------------------------------------
# flatten
# ---------------------------------------------------------------------------

class TestFlatten:
    def test_single_inner_list(self):
        assert CPFManager.flatten([["1", "2", "3"]]) == [1, 2, 3]

    def test_multiple_inner_lists(self):
        assert CPFManager.flatten([["1", "2"], ["3", "4"]]) == [1, 2, 3, 4]

    def test_empty_inner_lists(self):
        assert CPFManager.flatten([[], []]) == []

    def test_empty_outer_list(self):
        assert CPFManager.flatten([]) == []

    def test_single_element_lists(self):
        assert CPFManager.flatten([["10"], ["20"]]) == [10, 20]

    def test_already_integer_strings(self):
        result = CPFManager.flatten([["100", "200", "300"]])
        assert all(isinstance(x, int) for x in result)
        assert result == [100, 200, 300]

    def test_mixed_size_inner_lists(self):
        assert CPFManager.flatten([["1"], ["2", "3", "4"]]) == [1, 2, 3, 4]


# ---------------------------------------------------------------------------
# functor
# ---------------------------------------------------------------------------

class TestFunctor:
    def test_apply_int_to_flat_list(self):
        result = CPFManager.functor(int, ["1", "2", "3"])
        assert result == [1, 2, 3]

    def test_apply_int_to_nested_list(self):
        result = CPFManager.functor(int, [["1", "2"], ["3", "4"]])
        assert result == [[1, 2], [3, 4]]

    def test_apply_str_to_flat_list(self):
        result = CPFManager.functor(str, [1, 2, 3])
        assert result == ["1", "2", "3"]

    def test_apply_to_scalar(self):
        result = CPFManager.functor(int, "42")
        assert result == 42

    def test_apply_to_empty_list(self):
        result = CPFManager.functor(int, [])
        assert result == []

    def test_deeply_nested(self):
        result = CPFManager.functor(int, [[["1"]], [["2"]]])
        assert result == [[[1]], [[2]]]

    def test_apply_float(self):
        result = CPFManager.functor(float, ["1.5", "2.5"])
        assert result == [1.5, 2.5]


# ---------------------------------------------------------------------------
# selectfrag
# ---------------------------------------------------------------------------

class TestSelectfrag:
    def test_range_string(self):
        result = CPFManager.selectfrag("1-5")
        assert result == [1, 2, 3, 4, 5]

    def test_range_string_single(self):
        result = CPFManager.selectfrag("3-3")
        assert result == [3]

    def test_comma_separated_string(self):
        result = CPFManager.selectfrag("1,3,5")
        assert result == [1, 3, 5]

    def test_single_number_string(self):
        result = CPFManager.selectfrag("7")
        assert result == [7]

    def test_integer_input(self):
        result = CPFManager.selectfrag(42)
        assert result == 42

    def test_range_inclusive(self):
        result = CPFManager.selectfrag("2-4")
        assert result == [2, 3, 4]


# ---------------------------------------------------------------------------
# read_header (static method, needs a file-like object)
# ---------------------------------------------------------------------------

class TestReadHeader:
    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_read_header_rev23(self):
        with open(SAMPLE_CPF_23, "rt") as f:
            header, natom, nfrag, labels, cpfver = CPFManager.read_header(f)
        assert cpfver == 23
        assert natom == 38
        assert nfrag == 5
        assert "CPF Open1.0 rev23" in header
        assert "charge" in labels
        assert "DPM" in labels
        assert "monomer" in labels
        assert "dimer" in labels
        assert "MUL-HF" in labels["charge"]

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_10),
        reason="Sample CPF rev10 file not found",
    )
    def test_read_header_rev10(self):
        with open(SAMPLE_CPF_10, "rt") as f:
            header, natom, nfrag, labels, cpfver = CPFManager.read_header(f)
        assert cpfver == 10
        assert "charge" in labels
        assert "dimer" in labels


# ---------------------------------------------------------------------------
# parse (integration test using sample CPF file)
# ---------------------------------------------------------------------------

class TestParse:
    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_basic(self, cpfm):
        result = cpfm.parse(SAMPLE_CPF_23)
        # parse returns self
        assert result is cpfm
        assert cpfm.cpfver == 23

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_atominfo(self, cpfm):
        cpfm.parse(SAMPLE_CPF_23)
        import pandas as pd
        assert isinstance(cpfm.atominfo, pd.DataFrame)
        assert len(cpfm.atominfo) == 38
        assert "alabels" in cpfm.atominfo.columns
        assert "elems" in cpfm.atominfo.columns
        assert "xcoords" in cpfm.atominfo.columns

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_fraginfo(self, cpfm):
        cpfm.parse(SAMPLE_CPF_23)
        assert "natoms" in cpfm.fraginfo
        assert "baas" in cpfm.fraginfo
        assert "connects" in cpfm.fraginfo
        assert len(cpfm.fraginfo["natoms"]) == 5

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_diminfo(self, cpfm):
        cpfm.parse(SAMPLE_CPF_23)
        import pandas as pd
        assert isinstance(cpfm.diminfo, pd.DataFrame)
        # 5 fragments -> 5*4/2 = 10 dimer pairs
        assert len(cpfm.diminfo) == 10
        assert "fragi" in cpfm.diminfo.columns
        assert "fragj" in cpfm.diminfo.columns
        assert "min-dist" in cpfm.diminfo.columns

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_mominfo(self, cpfm):
        cpfm.parse(SAMPLE_CPF_23)
        import pandas as pd
        assert isinstance(cpfm.mominfo, pd.DataFrame)
        assert len(cpfm.mominfo) == 5

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_condition(self, cpfm):
        cpfm.parse(SAMPLE_CPF_23)
        assert cpfm.condition["basis_set"] == "STO-3G"
        assert cpfm.condition["electronic_state"] == "S1"
        assert cpfm.condition["calculation_method"] == "HF"

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_static_data(self, cpfm):
        cpfm.parse(SAMPLE_CPF_23)
        assert cpfm.static_data["natom"] == 38
        assert cpfm.static_data["nfrag"] == 5
        assert cpfm.static_data["ndimer"] == 10
        assert pytest.approx(cpfm.static_data["nuclear_repulsion_energy"],
                             rel=1e-6) == 1889.49334720085

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_parse_rev23_labels(self, cpfm):
        cpfm.parse(SAMPLE_CPF_23)
        assert "charge" in cpfm.labels
        assert "DPM" in cpfm.labels
        assert "monomer" in cpfm.labels
        assert "dimer" in cpfm.labels

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_10),
        reason="Sample CPF rev10 file not found",
    )
    def test_parse_rev10(self, cpfm):
        cpfm.parse(SAMPLE_CPF_10)
        assert cpfm.cpfver == 10


# ---------------------------------------------------------------------------
# write (round-trip test)
# ---------------------------------------------------------------------------

class TestWrite:
    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_write_creates_file(self, cpfm, tmp_path):
        cpfm.parse(SAMPLE_CPF_23)
        outfile = str(tmp_path / "out.cpf")
        cpfm.write("round-trip-test", outfile, cpfver=23)
        assert os.path.isfile(outfile)

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_write_contains_end(self, cpfm, tmp_path):
        cpfm.parse(SAMPLE_CPF_23)
        outfile = str(tmp_path / "out.cpf")
        cpfm.write("round-trip-test", outfile, cpfver=23)
        with open(outfile, "r") as f:
            content = f.read()
        assert content.rstrip().endswith("END")

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_CPF_23),
        reason="Sample CPF rev23 file not found",
    )
    def test_write_roundtrip_header(self, cpfm, tmp_path):
        cpfm.parse(SAMPLE_CPF_23)
        outfile = str(tmp_path / "out.cpf")
        cpfm.write("round-trip-test", outfile, cpfver=23)
        # Re-parse the written file
        cpfm2 = CPFManager()
        cpfm2.parse(outfile)
        assert cpfm2.cpfver == 23
        assert cpfm2.static_data["natom"] == cpfm.static_data["natom"]
        assert cpfm2.static_data["nfrag"] == cpfm.static_data["nfrag"]


# ---------------------------------------------------------------------------
# setupfragstr
# ---------------------------------------------------------------------------

class TestSetupfragstr:
    def test_basic_fragstr(self):
        fraginfo = {
            "natoms": [16, 30, 30, 30, 54],
            "baas": [0, 1, 1, 1, 1],
            "connects": [[2, 3], [10, 11], [17, 18], [24, 25]],
        }
        result = CPFManager.setupfragstr(fraginfo)
        # natoms and baas are formatted 8-wide
        assert "16" in result
        assert "30" in result
        assert "54" in result
        # connects formatted 10-wide pairs
        assert "2" in result
        assert "3" in result

    def test_fragstr_newline_per_10(self):
        fraginfo = {
            "natoms": list(range(1, 12)),  # 11 items -> 2 lines
            "baas": [0] * 11,
            "connects": [],
        }
        result = CPFManager.setupfragstr(fraginfo)
        lines = result.strip().split("\n")
        # natoms: 2 lines (10 + 1), baas: 2 lines (10 + 1)
        assert len(lines) == 4
