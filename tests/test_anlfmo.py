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


# ---------------------------------------------------------------------------
# Tests for depth()
# ---------------------------------------------------------------------------
class TestDepth:
    """Tests for anlfmo.depth — recursion depth helper."""

    def setup_method(self):
        self.obj = anlfmo()

    def test_depth_none(self):
        assert self.obj.depth(None) == 0

    def test_depth_empty_list(self):
        assert self.obj.depth([]) == 0

    def test_depth_flat_list(self):
        assert self.obj.depth([1, 2, 3]) == 1

    def test_depth_nested_list(self):
        assert self.obj.depth([[1, 2], [3, 4]]) == 2

    def test_depth_deeply_nested(self):
        assert self.obj.depth([[[1]]]) == 3

    def test_depth_scalar(self):
        assert self.obj.depth(42) == 0

    def test_depth_mixed_nesting(self):
        # max depth wins: [1] is depth 1, [[2]] is depth 2
        assert self.obj.depth([[2], [1]]) == 2


# ---------------------------------------------------------------------------
# Tests for getisdisp()
# ---------------------------------------------------------------------------
class TestGetIsDisp:
    """Tests for anlfmo.getisdisp — checks dispersion flag in log."""

    def setup_method(self):
        self.obj = anlfmo()

    def test_disp_on(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Some header\nDisp = ON\nMore stuff\n")
        assert self.obj.getisdisp(str(log)) is True

    def test_disp_absent(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Some header\nMethod = MP2\nMore stuff\n")
        assert self.obj.getisdisp(str(log)) is False

    def test_disp_off_not_matched(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Disp = OFF\n")
        assert self.obj.getisdisp(str(log)) is False


# ---------------------------------------------------------------------------
# Tests for getlogmethod()
# ---------------------------------------------------------------------------
class TestGetLogMethod:
    """Tests for anlfmo.getlogmethod — gets method from log."""

    def setup_method(self):
        self.obj = anlfmo()

    def test_method_mp2(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Header line\nMethod = MP2\nAnother line\n")
        assert self.obj.getlogmethod(str(log)) == 'MP2'

    def test_method_hf(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Method = HF\n")
        assert self.obj.getlogmethod(str(log)) == 'HF'

    def test_method_missing_returns_none(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("No method line here\n")
        assert self.obj.getlogmethod(str(log)) is None


# ---------------------------------------------------------------------------
# Tests for getpbflag()
# ---------------------------------------------------------------------------
class TestGetPbFlag:
    """Tests for anlfmo.getpbflag — checks PB flag in log."""

    def setup_method(self):
        self.obj = anlfmo()

    def test_pb_on(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Header\nEFFECT = ON\n## CHECK AVAILABLE\nTrailing\n")
        self.obj.getpbflag(str(log))
        assert self.obj.pbflag is True

    def test_pb_off(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Header\nEFFECT = OFF\n## CHECK AVAILABLE\n")
        self.obj.getpbflag(str(log))
        assert self.obj.pbflag is False

    def test_pb_absent(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("Nothing relevant\n## CHECK AVAILABLE\n")
        self.obj.getpbflag(str(log))
        assert self.obj.pbflag is False


# ---------------------------------------------------------------------------
# DataFrame helper tests — need pandas
# ---------------------------------------------------------------------------
import pandas as pd


class TestGetIfiedf:
    """Tests for anlfmo.getifiedf — creates DataFrame from IFIE data."""

    def setup_method(self):
        self.obj = anlfmo()
        # icolumn: ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE',
        #           'PR-TYPE1', 'GRIMME', 'JUNG', 'HILL']
        self.sample_ifie = [
            ['1', '2', '3.5', 'T', '0.001', '0.002', '0.003', '0.004', '0.005', '0.006'],
            ['1', '3', '5.0', 'F', '0.010', '0.020', '0.030', '0.040', '0.050', '0.060'],
        ]

    def test_returns_dataframe(self):
        df = self.obj.getifiedf(self.sample_ifie)
        assert isinstance(df, pd.DataFrame)

    def test_column_names(self):
        df = self.obj.getifiedf(self.sample_ifie)
        assert list(df.columns) == self.obj.icolumn

    def test_row_count(self):
        df = self.obj.getifiedf(self.sample_ifie)
        assert len(df) == 2

    def test_int_columns(self):
        df = self.obj.getifiedf(self.sample_ifie)
        assert df['I'].dtype in (int, 'int64')
        assert df['J'].dtype in (int, 'int64')

    def test_hf_ifie_converted_to_kcal(self):
        df = self.obj.getifiedf(self.sample_ifie)
        # 0.001 * 627.5095
        assert df.iloc[0]['HF-IFIE'] == pytest.approx(0.001 * 627.5095)

    def test_mp2_ifie_converted_to_kcal(self):
        df = self.obj.getifiedf(self.sample_ifie)
        assert df.iloc[0]['MP2-IFIE'] == pytest.approx(0.002 * 627.5095)

    def test_solv_merge(self):
        solv = [['1', '2', '0.123'], ['1', '3', '0.456']]
        df = self.obj.getifiedf(self.sample_ifie, solv=solv)
        assert 'Solv(ES)' in df.columns
        assert df.iloc[0]['Solv(ES)'] == pytest.approx(0.123)


class TestGetPiedadf:
    """Tests for anlfmo.getpiedadf — creates DataFrame from PIEDA data."""

    def setup_method(self):
        self.obj = anlfmo()
        # pcolumn: ['I', 'J', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']
        self.sample_pieda = [
            ['1', '2', '-1.5', '0.3', '-0.2', '-0.8', '0.01'],
            ['1', '3', '-2.0', '0.5', '-0.1', '-1.0', '0.02'],
        ]

    def test_returns_dataframe(self):
        df = self.obj.getpiedadf(self.sample_pieda)
        assert isinstance(df, pd.DataFrame)

    def test_column_types(self):
        df = self.obj.getpiedadf(self.sample_pieda)
        assert df['ES'].dtype == float
        assert df['EX'].dtype == float
        assert df['DI(MP2)'].dtype == float

    def test_values(self):
        df = self.obj.getpiedadf(self.sample_pieda)
        assert df.iloc[0]['ES'] == pytest.approx(-1.5)
        assert df.iloc[1]['q(I=>J)'] == pytest.approx(0.02)


class TestGetMomenedf:
    """Tests for anlfmo.getmomenedf — creates DataFrame from monomer energy."""

    def setup_method(self):
        self.obj = anlfmo()

    def test_basic(self):
        data = [['1', '-75.123', '-75.456'], ['2', '-80.111', '-80.222']]
        df = self.obj.getmomenedf(data)
        assert list(df.columns) == ['Frag.', 'HF', 'MP2']
        assert len(df) == 2
        assert df.iloc[0]['Frag.'] == 1
        assert df.iloc[0]['HF'] == pytest.approx(-75.123)
        assert df.iloc[1]['MP2'] == pytest.approx(-80.222)


class TestGetDimenedf:
    """Tests for anlfmo.getdimenedf — creates DataFrame from dimer energy."""

    def setup_method(self):
        self.obj = anlfmo()

    def test_basic(self):
        data = [['1', '2', '-150.0', '-150.5'], ['1', '3', '-155.0', '-155.5']]
        df = self.obj.getdimenedf(data)
        assert list(df.columns) == ['I', 'J', 'DIMER-HF', 'DIMER-MP2']
        assert len(df) == 2
        assert df.iloc[0]['I'] == 1
        assert df.iloc[1]['DIMER-MP2'] == pytest.approx(-155.5)


# ---------------------------------------------------------------------------
# Tests for gettgtdf_ff and gettgtdf_ffs
# ---------------------------------------------------------------------------
class TestGettgtdfFf:
    """Tests for anlfmo.gettgtdf_ff — filter DataFrame by fragment pair."""

    def setup_method(self):
        self.obj = anlfmo()
        self.df = pd.DataFrame({
            'I': [1, 1, 2, 3],
            'J': [2, 3, 3, 4],
            'HF-IFIE': [0.1, 0.2, 0.3, 0.4],
        })

    def test_finds_pair_ij(self):
        result = self.obj.gettgtdf_ff(self.df, 1, 2)
        assert len(result) == 1
        assert result.iloc[0]['HF-IFIE'] == pytest.approx(0.1)

    def test_finds_pair_ji(self):
        # frag1=2, frag2=1 should still find the (1,2) row
        result = self.obj.gettgtdf_ff(self.df, 2, 1)
        assert len(result) == 1

    def test_no_match(self):
        result = self.obj.gettgtdf_ff(self.df, 1, 4)
        assert len(result) == 0

    def test_string_args_converted(self):
        result = self.obj.gettgtdf_ff(self.df, '1', '3')
        assert len(result) == 1


class TestGettgtdfFfs:
    """Tests for anlfmo.gettgtdf_ffs — filter by frag1 vs list of frag2."""

    def setup_method(self):
        self.obj = anlfmo()
        self.df = pd.DataFrame({
            'I': [1, 1, 2, 3],
            'J': [2, 3, 3, 4],
            'HF-IFIE': [0.1, 0.2, 0.3, 0.4],
        })

    def test_multiple_targets(self):
        result = self.obj.gettgtdf_ffs(self.df, 1, [2, 3])
        assert len(result) == 2

    def test_single_target_list(self):
        result = self.obj.gettgtdf_ffs(self.df, 1, [2])
        assert len(result) == 1

    def test_symmetric_lookup(self):
        # frag2 contains 1, frag1=2 should find the (1,2) row via swapped condition
        result = self.obj.gettgtdf_ffs(self.df, 2, [1])
        assert len(result) == 1
