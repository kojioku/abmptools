# -*- coding: utf-8 -*-
"""Tests for abmptools.gro2udf.mdp_parser."""
import pytest

from abmptools.gro2udf.mdp_parser import parse_mdp, MdpParams, load_mdp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

MINIMAL_MDP = """\
; GROMACS mdp file
integrator  = md
ref-t       = 300.0
tau-t       = 0.5
coulombtype = PME
rcoulomb    = 1.2
constraints = h-bonds
"""


# ---------------------------------------------------------------------------
# parse_mdp tests
# ---------------------------------------------------------------------------

class TestParseMdp:

    def test_basic_parsing(self, tmp_path):
        mdp = tmp_path / "run.mdp"
        mdp.write_text(MINIMAL_MDP)
        d = parse_mdp(str(mdp))
        assert d["integrator"] == "md"
        assert d["rcoulomb"] == "1.2"

    def test_comment_lines_skipped(self, tmp_path):
        mdp = tmp_path / "run.mdp"
        mdp.write_text("; only comments\n; nothing else\n")
        d = parse_mdp(str(mdp))
        assert d == {}

    def test_hyphen_normalization(self, tmp_path):
        mdp = tmp_path / "run.mdp"
        mdp.write_text("ref-t = 400.0\ntau-t = 0.2\n")
        d = parse_mdp(str(mdp))
        assert "ref_t" in d
        assert "tau_t" in d

    def test_inline_comment_stripped(self, tmp_path):
        mdp = tmp_path / "run.mdp"
        mdp.write_text("ref_t = 300.0 ; temperature\n")
        d = parse_mdp(str(mdp))
        assert d["ref_t"] == "300.0"


# ---------------------------------------------------------------------------
# MdpParams tests
# ---------------------------------------------------------------------------

class TestMdpParams:

    def test_defaults_when_keys_missing(self):
        p = MdpParams({})
        assert p.ref_t == pytest.approx(300.0)
        assert p.tau_t == pytest.approx(0.1)
        assert p.coulombtype == "PME"
        assert p.rcoulomb == pytest.approx(0.9)
        assert p.constraints == "none"
        assert p.integrator == "md"

    def test_space_separated_ref_t(self):
        """When ref_t contains multiple values, the first is taken."""
        p = MdpParams({"ref_t": "300.0 300.0"})
        assert p.ref_t == pytest.approx(300.0)


# ---------------------------------------------------------------------------
# load_mdp tests
# ---------------------------------------------------------------------------

class TestLoadMdp:

    def test_none_returns_none(self):
        assert load_mdp(None) is None

    def test_valid_path_returns_mdp_params(self, tmp_path):
        mdp = tmp_path / "run.mdp"
        mdp.write_text(MINIMAL_MDP)
        result = load_mdp(str(mdp))
        assert isinstance(result, MdpParams)
        assert result.ref_t == pytest.approx(300.0)
        assert result.constraints == "h-bonds"


# ---------------------------------------------------------------------------
# MdpParams: time-related accessors (dt / nsteps / nstxout-compressed / nstenergy)
# ---------------------------------------------------------------------------

def test_mdp_params_dt_default():
    """MdpParams.dt は dt が未設定なら 0.001 ps default。"""
    from abmptools.gro2udf.mdp_parser import MdpParams
    assert MdpParams({}).dt == 0.001


def test_mdp_params_dt_parsed(tmp_path):
    from abmptools.gro2udf.mdp_parser import MdpParams
    assert MdpParams({"dt": "0.002"}).dt == 0.002


def test_mdp_params_nsteps_default():
    from abmptools.gro2udf.mdp_parser import MdpParams
    assert MdpParams({}).nsteps == 0


def test_mdp_params_nsteps_parsed():
    from abmptools.gro2udf.mdp_parser import MdpParams
    assert MdpParams({"nsteps": "500000"}).nsteps == 500000


def test_mdp_params_nstxout_compressed_prefers_dashed_key():
    """`nstxout-compressed = 5000` を優先して読む。"""
    from abmptools.gro2udf.mdp_parser import MdpParams
    p = MdpParams({"nstxout-compressed": "5000",
                   "nstxout_compressed": "0",
                   "nstxout": "0"})
    assert p.nstxout_compressed == 5000


def test_mdp_params_nstxout_compressed_falls_back_to_nstxout():
    """nstxout-compressed が 0 / 未設定なら nstxout に fallback。"""
    from abmptools.gro2udf.mdp_parser import MdpParams
    p = MdpParams({"nstxout-compressed": "0", "nstxout": "1000"})
    assert p.nstxout_compressed == 1000


def test_mdp_params_nstenergy_parsed():
    from abmptools.gro2udf.mdp_parser import MdpParams
    assert MdpParams({"nstenergy": "1000"}).nstenergy == 1000


def test_mdp_params_nstenergy_default():
    from abmptools.gro2udf.mdp_parser import MdpParams
    assert MdpParams({}).nstenergy == 0
