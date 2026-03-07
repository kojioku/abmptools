"""Tests for abmptools.setfmo.setfmo class."""
import pytest
from abmptools.setfmo import setfmo


@pytest.fixture
def fmo():
    """Create a setfmo instance for testing."""
    return setfmo()


# ---------- __init__ defaults ----------

class TestInit:
    def test_ajf_method_default(self, fmo):
        assert fmo.ajf_method == 'HF'

    def test_ajf_basis_set_default(self, fmo):
        assert fmo.ajf_basis_set == '6-31Gdag'

    def test_cpfflag_default(self, fmo):
        assert fmo.cpfflag is True

    def test_solv_flag_default(self, fmo):
        assert fmo.solv_flag is False

    def test_piedaflag_default(self, fmo):
        assert fmo.piedaflag is True

    def test_cutmode_default(self, fmo):
        assert fmo.cutmode == 'sphere'

    def test_abinit_ver_default(self, fmo):
        assert fmo.abinit_ver == 'rev22'

    def test_memory_default(self, fmo):
        assert fmo.memory == 3000

    def test_npro_default(self, fmo):
        assert fmo.npro == 8

    def test_para_job_default(self, fmo):
        assert fmo.para_job == 1

    def test_mainpath_default(self, fmo):
        assert fmo.mainpath == '.'

    def test_cmmflag_default(self, fmo):
        assert fmo.cmmflag is False

    def test_molname_default(self, fmo):
        assert fmo.molname == []

    def test_criteria_default(self, fmo):
        assert fmo.criteria == []

    def test_tgtpos_default(self, fmo):
        assert fmo.tgtpos == []


# ---------- setrfmoparam ----------

class TestSetrfmoparam:
    def test_empty_dict_gives_all_defaults(self, fmo):
        fmo.setrfmoparam({})
        assert fmo.cutmode == 'sphere'
        assert fmo.getmode == 'resnum'
        assert fmo.ajf_method == 'HF'
        assert fmo.ajf_basis_set == '6-31Gdag'
        assert fmo.cpfflag is True
        assert fmo.solv_flag is False
        assert fmo.piedaflag is True
        assert fmo.cmmflag is False
        assert fmo.abinit_ver == 'rev22'
        assert fmo.molname == []
        assert fmo.criteria == []
        assert fmo.tgtpos == []
        assert fmo.solutes == []
        assert fmo.solutes_charge == 0
        assert fmo.ionname == []
        assert fmo.ionmode == 'remain'
        assert fmo.memory == 3000
        assert fmo.npro == 8

    def test_all_keys_provided(self, fmo):
        params = {
            'cutmode': 'cube',
            'getmode': 'atomnum',
            'ajf_method': 'MP2',
            'ajf_basis_set': 'cc-pVDZ',
            'cpfflag': False,
            'solv_flag': True,
            'piedaflag': False,
            'cmmflag': True,
            'abinit_ver': 'rev23',
            'molname': ['WAT', 'NA'],
            'criteria': [5.0],
            'tgtpos': [1.0, 2.0, 3.0],
            'solutes': [0, 1],
            'solutes_charge': 2,
            'ionname': ['NA', 'CL'],
            'ionmode': 'remove',
            'memory': 8000,
            'npro': 16,
        }
        fmo.setrfmoparam(params)

        assert fmo.cutmode == 'cube'
        assert fmo.getmode == 'atomnum'
        assert fmo.ajf_method == 'MP2'
        assert fmo.ajf_basis_set == 'cc-pVDZ'
        assert fmo.cpfflag is False
        assert fmo.solv_flag is True
        assert fmo.piedaflag is False
        assert fmo.cmmflag is True
        assert fmo.abinit_ver == 'rev23'
        assert fmo.molname == ['WAT', 'NA']
        assert fmo.criteria == [5.0]
        assert fmo.tgtpos == [1.0, 2.0, 3.0]
        assert fmo.solutes == [0, 1]
        assert fmo.solutes_charge == 2
        assert fmo.ionname == ['NA', 'CL']
        assert fmo.ionmode == 'remove'
        assert fmo.memory == 8000
        assert fmo.npro == 16

    def test_partial_keys_uses_defaults_for_missing(self, fmo):
        params = {
            'ajf_method': 'DFT',
            'memory': 4000,
            'solutes': [0],
        }
        fmo.setrfmoparam(params)

        # Provided values
        assert fmo.ajf_method == 'DFT'
        assert fmo.memory == 4000
        assert fmo.solutes == [0]

        # Defaults for missing keys
        assert fmo.cutmode == 'sphere'
        assert fmo.getmode == 'resnum'
        assert fmo.ajf_basis_set == '6-31Gdag'
        assert fmo.cpfflag is True
        assert fmo.solv_flag is False
        assert fmo.piedaflag is True
        assert fmo.cmmflag is False
        assert fmo.abinit_ver == 'rev22'
        assert fmo.molname == []
        assert fmo.criteria == []
        assert fmo.tgtpos == []
        assert fmo.solutes_charge == 0
        assert fmo.ionname == []
        assert fmo.ionmode == 'remain'
        assert fmo.npro == 8

    def test_setrfmoparam_overwrites_init_values(self, fmo):
        """Calling setrfmoparam should overwrite values set by __init__."""
        assert fmo.ajf_method == 'HF'
        fmo.setrfmoparam({'ajf_method': 'MP2'})
        assert fmo.ajf_method == 'MP2'

    def test_setrfmoparam_called_twice_second_overrides(self, fmo):
        """A second call to setrfmoparam should fully reset to its own values."""
        fmo.setrfmoparam({'ajf_method': 'MP2', 'npro': 32})
        assert fmo.ajf_method == 'MP2'
        assert fmo.npro == 32

        fmo.setrfmoparam({})
        assert fmo.ajf_method == 'HF'
        assert fmo.npro == 8
