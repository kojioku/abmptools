"""Tests for abmptools.logmanager.LOGManager."""

import io
import os
import pytest

from abmptools.logmanager import LOGManager


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, "sample", "log2cpf"
)

SAMPLE_HF_LOG = os.path.join(SAMPLE_DIR, "gly5-HF-STO-3G-resp.log")
SAMPLE_MP2_LOG = os.path.join(SAMPLE_DIR, "gly5-MP2-STO-3G-nbo.log")


@pytest.fixture
def manager():
    """Return a fresh LOGManager instance."""
    return LOGManager()


# ---------------------------------------------------------------------------
# __init__ defaults
# ---------------------------------------------------------------------------

class TestInit:

    def test_cpfver_default(self, manager):
        assert manager.cpfver == 23

    def test_is_gz_default(self, manager):
        assert manager.is_gz is False

    def test_atominfo_default(self, manager):
        assert manager.atominfo == {}

    def test_fraginfo_default(self, manager):
        assert manager.fraginfo == []

    def test_condition_default(self, manager):
        assert manager.condition == []

    def test_static_data_default(self, manager):
        assert manager.static_data == {}

    def test_mominfo_default(self, manager):
        assert manager.mominfo == {}

    def test_diminfo_default(self, manager):
        assert manager.diminfo == {}

    def test_labels_default(self, manager):
        assert manager.labels == {}

    def test_tgtfrag_default(self, manager):
        assert manager.tgtfrag == 0


# ---------------------------------------------------------------------------
# getversion with StringIO
# ---------------------------------------------------------------------------

class TestGetVersion:

    def test_version1(self):
        """Version 1 log header (Open Ver. 1 Rev. 23)."""
        text = (
            "\n"
            "   ----------------------------------------------\n"
            "     ABINIT-MP - Open Ver. 1 Rev. 23 (BINDS Ver. 1) / 20230712\n"
            "   ----------------------------------------------\n"
        )
        f = io.StringIO(text)
        version = LOGManager.getversion(f)
        assert version == 1

    def test_version2(self):
        """Version 2 log header (Open Ver. 2 Rev. 8)."""
        text = (
            "\n"
            "   ----------------------------------------------\n"
            "     ABINIT-MP - Open Ver. 2 Rev. 8 / 20240101\n"
            "   ----------------------------------------------\n"
        )
        f = io.StringIO(text)
        version = LOGManager.getversion(f)
        assert version == 2


# ---------------------------------------------------------------------------
# getcondition with StringIO
# ---------------------------------------------------------------------------

class TestGetCondition:

    @staticmethod
    def _make_condition_block(
        method="HF",
        elecstate="S1",
        basisset="STO-3G",
        laoc="0.000",
        lptc="2.000",
        ldimer="2.000",
        readgeom="test.pdb",
        autofrag="ON",
        nprint="3",
        nboanl="OFF",
        esptyp="NONE",
    ):
        """Build a minimal condition section that getcondition can parse."""
        lines = [
            f"        Method                  = {method}",
            f"        Nprint                  = {nprint:>5s}",
            f"        ElecState               = {elecstate}",
            f"        BasisSet                = {basisset}",
            f"        Laoc                    = {laoc:>9s}",
            f"        Lptc                    = {lptc:>9s}",
            f"        ReadGeom                = {readgeom}",
            f"        AutoFrag                = {autofrag}",
            f"        Ldimer                  = {ldimer:>9s}",
            f"        NBOANL                  = {nboanl}",
            f"        ESPTYP                  = {esptyp}",
            "     ## CHECK AVAILABLE MEMORY",
        ]
        return io.StringIO("\n".join(lines) + "\n")

    def test_hf_defaults(self):
        f = self._make_condition_block()
        (Method, ElecState, BasisSet, Laoc, Lptc, Ldimer,
         ReadGeom, fragmode, is_npa, is_resp, Nprint) = LOGManager.getcondition(f)
        assert Method == "HF"
        assert ElecState == "S1"
        assert BasisSet == "STO-3G"
        assert Laoc == 0.0
        assert Lptc == 2.0
        assert Ldimer == 2.0
        assert ReadGeom == "test.pdb"
        assert fragmode == "auto"
        assert is_npa is False
        assert is_resp is False
        assert Nprint == 3

    def test_mp2_with_npa_resp(self):
        f = self._make_condition_block(
            method="MP2", nboanl="ON", esptyp="RESP", nprint="1"
        )
        (Method, _, _, _, _, _, _, _, is_npa, is_resp, Nprint) = \
            LOGManager.getcondition(f)
        assert Method == "MP2"
        assert is_npa is True
        assert is_resp is True
        assert Nprint == 1

    def test_manual_fragmode(self):
        f = self._make_condition_block(autofrag="OFF")
        (_, _, _, _, _, _, _, fragmode, _, _, _) = LOGManager.getcondition(f)
        assert fragmode == "manual"


# ---------------------------------------------------------------------------
# readelemslog with StringIO
# ---------------------------------------------------------------------------

class TestReadElemsLog:

    def test_basic_elements(self):
        """Parse a small element block."""
        lines = [
            "  ==================================",
            "     ## READ MOLECULAR STRUCTURE FROM PDB-FILE",
            "  ==================================",
            "       1   N   nitrogen",
            "       2   C   carbon",
            "       3   H   hydrogen",
            "     ## Molecular formula",
        ]
        f = io.StringIO("\n".join(lines) + "\n")
        elems = LOGManager.readelemslog(f)
        assert elems == ["N", "C", "H"]


# ---------------------------------------------------------------------------
# Integration tests with sample log files
# ---------------------------------------------------------------------------

class TestParseIntegration:

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_HF_LOG),
        reason="Sample HF log file not available",
    )
    def test_parse_hf_log(self, monkeypatch):
        monkeypatch.chdir(SAMPLE_DIR)
        mgr = LOGManager()
        mgr.parse(SAMPLE_HF_LOG)

        # condition
        assert mgr.condition["calculation_method"] == "HF"
        assert mgr.condition["basis_set"] == "STO-3G"
        assert mgr.condition["electronic_state"] == "S1"

        # static_data
        assert mgr.static_data["nfrag"] == 5
        assert mgr.static_data["natom"] == 38
        assert mgr.static_data["ndimer"] == 10

        # labels
        assert "HF" in mgr.labels["monomer"]
        assert "ES" in mgr.labels["dimer"]

        # atominfo should be a DataFrame with correct row count
        assert len(mgr.atominfo) == 38

        # mominfo should have one row per fragment
        assert len(mgr.mominfo) == 5

        # diminfo should have one row per dimer pair
        assert len(mgr.diminfo) == 10

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_MP2_LOG),
        reason="Sample MP2 log file not available",
    )
    def test_parse_mp2_log(self, monkeypatch):
        monkeypatch.chdir(SAMPLE_DIR)
        mgr = LOGManager()
        mgr.parse(SAMPLE_MP2_LOG)

        assert mgr.condition["calculation_method"] == "MP2"
        assert mgr.static_data["nfrag"] == 5
        assert "MP2" in mgr.labels["monomer"]
        assert "MP2" in mgr.labels["dimer"]

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_HF_LOG),
        reason="Sample HF log file not available",
    )
    def test_getversion_from_real_file(self):
        with open(SAMPLE_HF_LOG, "r") as f:
            version = LOGManager.getversion(f)
        assert version == 1

    @pytest.mark.skipif(
        not os.path.isfile(SAMPLE_HF_LOG),
        reason="Sample HF log file not available",
    )
    def test_getcondition_from_real_file(self):
        with open(SAMPLE_HF_LOG, "r") as f:
            # Skip past the version header first (getversion consumes lines)
            _ = LOGManager.getversion(f)
            result = LOGManager.getcondition(f)
        Method, ElecState, BasisSet, Laoc, Lptc, Ldimer, \
            ReadGeom, fragmode, is_npa, is_resp, Nprint = result
        assert Method == "HF"
        assert ElecState == "S1"
        assert BasisSet == "STO-3G"
        assert fragmode == "auto"
        assert is_resp is True
        assert Nprint == 3
