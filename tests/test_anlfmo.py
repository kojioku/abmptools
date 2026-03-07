"""Tests for abmptools.anlfmo module."""

import pytest
from unittest.mock import patch, MagicMock

from abmptools.anlfmo import anlfmo


class TestAnlfmoInit:
    """Tests for anlfmo.__init__ default values."""

    def test_anlmode_default(self):
        obj = anlfmo()
        assert obj.anlmode == 'frag'

    def test_fragmode_default(self):
        obj = anlfmo()
        assert obj.fragmode == 'auto'

    def test_cpfflag_default(self):
        obj = anlfmo()
        assert obj.cpfflag is True

    def test_solv_flag_default(self):
        obj = anlfmo()
        assert obj.solv_flag is False

    def test_memory_default(self):
        obj = anlfmo()
        assert obj.memory == 3000

    def test_npro_default(self):
        obj = anlfmo()
        assert obj.npro == 8

    def test_para_job_default(self):
        obj = anlfmo()
        assert obj.para_job == 1

    def test_cutmode_default(self):
        obj = anlfmo()
        assert obj.cutmode == 'sphere'

    def test_abinit_ver_default(self):
        obj = anlfmo()
        assert obj.abinit_ver == 'rev11'

    def test_piedaflag_default(self):
        obj = anlfmo()
        assert obj.piedaflag is True

    def test_dist_default(self):
        obj = anlfmo()
        assert obj.dist == 1000.0

    def test_tgt1frag_default(self):
        obj = anlfmo()
        assert obj.tgt1frag is None

    def test_rpdbflag_default(self):
        obj = anlfmo()
        assert obj.rpdbflag is False

    def test_pdbname_default(self):
        obj = anlfmo()
        assert obj.pdbname is None

    def test_is_disp_default(self):
        obj = anlfmo()
        assert obj.is_disp is False

    def test_tgt2type_default(self):
        obj = anlfmo()
        assert obj.tgt2type == 'frag'

    def test_selecttype_default(self):
        obj = anlfmo()
        assert obj.selecttype == 'molid'

    def test_tgtmolid_default(self):
        obj = anlfmo()
        assert obj.tgtmolid is None

    def test_logMethod_default(self):
        obj = anlfmo()
        assert obj.logMethod == 'MP2'

    def test_molname_default_empty_list(self):
        obj = anlfmo()
        assert obj.molname == []

    def test_criteria_default_empty_list(self):
        obj = anlfmo()
        assert obj.criteria == []

    def test_icolumn_has_expected_entries(self):
        obj = anlfmo()
        assert 'I' in obj.icolumn
        assert 'J' in obj.icolumn
        assert 'HF-IFIE' in obj.icolumn
        assert 'MP2-IFIE' in obj.icolumn

    def test_pcolumn_has_expected_entries(self):
        obj = anlfmo()
        assert 'ES' in obj.pcolumn
        assert 'EX' in obj.pcolumn
        assert 'CT-mix' in obj.pcolumn
        assert 'DI(MP2)' in obj.pcolumn

    def test_pynp_default(self):
        obj = anlfmo()
        assert obj.pynp == 3

    def test_multi_file_settings_default_none(self):
        obj = anlfmo()
        assert obj.ilog_head is None
        assert obj.ilog_tail is None
        assert obj.pdb_head is None
        assert obj.pdb_tail is None
        assert obj.start is None
        assert obj.end is None
        assert obj.interval is None

    def test_inherits_pdb_io_attributes(self):
        """Verify that pdb_io parent attributes are present via super().__init__."""
        obj = anlfmo()
        # pdb_io sets these in its __init__
        assert hasattr(obj, 'solutes')
        assert hasattr(obj, 'getmode')
        assert obj.getmode == 'resnum'
