# -*- coding: utf-8 -*-
"""Tests for abmptools.udf2gro.gromacs.writers.mdp_writer."""
import pytest

from abmptools.core.system_model import (
    CellGeometry,
    SimulationParams,
    SystemModel,
)
from abmptools.udf2gro.gromacs.writers.mdp_writer import MdpWriter


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _minimal_model(**sp_overrides) -> SystemModel:
    """Build a SystemModel with SimulationParams for MDP testing."""
    sp_defaults = dict(
        title="test mdp",
        algorithm="NVT_Nose_Hoover",
        nsteps=10000,
        dt=0.001,
        outputinterval=1000,
        outputinterval2=100,
        integrator="md-vv",
    )
    sp_defaults.update(sp_overrides)
    sp = SimulationParams(**sp_defaults)

    return SystemModel(
        title="test mdp",
        udf_path="/tmp/test.udf",
        comb_rule=2,
        flags14=0,
        fudgeLJ=0.5,
        fudgeQQ=0.8333,
        calcQQ=0,
        sim_params=sp,
        cell=CellGeometry(a=3.0, b=3.0, c=3.0),
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestMdpWriter:

    def test_write_creates_file(self, tmp_path):
        model = _minimal_model()
        filepath = str(tmp_path / "out.mdp")
        MdpWriter().write(model, filepath)

        assert (tmp_path / "out.mdp").exists()

    def test_output_contains_title_and_integrator(self, tmp_path):
        model = _minimal_model()
        filepath = str(tmp_path / "out.mdp")
        MdpWriter().write(model, filepath)

        content = open(filepath).read()
        assert "title" in content
        assert "integrator" in content
        assert "md-vv" in content

    def test_nvt_settings(self, tmp_path):
        """When t_coupl is not 'no', tc_grps / tau_t / ref_t must appear."""
        model = _minimal_model(
            t_coupl="nose-hoover",
            tau_t=0.5,
            ref_t=350.0,
        )
        filepath = str(tmp_path / "out.mdp")
        MdpWriter().write(model, filepath)

        content = open(filepath).read()
        assert "tc_grps" in content
        assert "tau_t" in content
        assert "ref_t" in content
        assert "350.0" in content

    def test_npt_settings(self, tmp_path):
        """When p_coupl is not 'no', pcoupl / tau_p / ref_p must appear."""
        model = _minimal_model(
            t_coupl="nose-hoover",
            p_coupl="MTTK",
            tau_p=2.0,
            ref_p=1.0,
        )
        filepath = str(tmp_path / "out.mdp")
        MdpWriter().write(model, filepath)

        content = open(filepath).read()
        assert "pcoupl" in content
        assert "tau_p" in content
        assert "ref_p" in content

    def test_vel_gen_true(self, tmp_path):
        """vel_gen=True should produce 'gen_vel = yes' in output."""
        model = _minimal_model(vel_gen=True, gen_temp=300.0)
        filepath = str(tmp_path / "out.mdp")
        MdpWriter().write(model, filepath)

        content = open(filepath).read()
        assert "gen_vel" in content
        assert "yes" in content
        assert "gen_temp" in content

    def test_constraints_rattle_bond(self, tmp_path):
        """rattle_bond=True should produce 'all-bonds' constraint."""
        model = _minimal_model(rattle_bond=True)
        filepath = str(tmp_path / "out.mdp")
        MdpWriter().write(model, filepath)

        content = open(filepath).read()
        assert "all-bonds" in content
