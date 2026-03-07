# -*- coding: utf-8 -*-
"""Tests for abmptools.core.system_model dataclasses."""

import pytest

from abmptools.core.system_model import (
    AtomType,
    AtomRecord,
    BondRecord,
    PairRecord,
    AngleRecord,
    DihedralRecord,
    MoleculeTopology,
    AtomPosition,
    CellGeometry,
    NdxData,
    SimulationParams,
    AnnealProtocol,
    SystemModel,
)


# ---------------------------------------------------------------------------
# CellGeometry
# ---------------------------------------------------------------------------

class TestCellGeometry:

    def test_defaults_alpha_beta_gamma(self):
        cell = CellGeometry(a=1.0, b=2.0, c=3.0)
        assert cell.alpha == pytest.approx(90.0)
        assert cell.beta == pytest.approx(90.0)
        assert cell.gamma == pytest.approx(90.0)

    def test_is_rectangular_with_defaults(self):
        cell = CellGeometry(a=1.0, b=2.0, c=3.0)
        assert cell.is_rectangular() is True

    def test_is_rectangular_false_for_triclinic(self):
        cell = CellGeometry(a=1.0, b=2.0, c=3.0, alpha=80.0, beta=90.0, gamma=90.0)
        assert cell.is_rectangular() is False

    def test_is_rectangular_within_threshold(self):
        cell = CellGeometry(a=1.0, b=1.0, c=1.0,
                            alpha=90.0 + 1e-7, beta=90.0 - 1e-7, gamma=90.0)
        assert cell.is_rectangular() is True

    def test_is_rectangular_custom_threshold(self):
        cell = CellGeometry(a=1.0, b=1.0, c=1.0, alpha=90.5)
        assert cell.is_rectangular(thres=1.0) is True
        assert cell.is_rectangular(thres=0.1) is False


# ---------------------------------------------------------------------------
# SimulationParams defaults
# ---------------------------------------------------------------------------

class TestSimulationParams:

    def test_defaults(self):
        sp = SimulationParams(
            title="test", algorithm="NVT", nsteps=1000, dt=0.001,
            outputinterval=100, outputinterval2=10, integrator="md-vv",
        )
        assert sp.vel_gen is False
        assert sp.calcQQ == 0
        assert sp.t_coupl == "no"
        assert sp.ref_t == pytest.approx(300.0)
        assert sp.tau_t == pytest.approx(0.1)
        assert sp.p_coupl == "no"
        assert sp.pcoupltype == "isotropic"
        assert sp.tau_p == pytest.approx(2.0)
        assert sp.ref_p == pytest.approx(1.0)
        assert sp.compressibility == pytest.approx(4.5e-5)
        assert sp.tail_correction == 0
        assert sp.pbc == "xyz"
        assert sp.periodic_mol is False
        assert sp.lj_cutoff == pytest.approx(0.0)
        assert sp.coulomb_cutoff == pytest.approx(0.0)
        assert sp.ld_seed is None
        assert sp.gen_temp is None
        assert sp.deform_vel is None
        assert sp.freeze_grps is None


# ---------------------------------------------------------------------------
# AnnealProtocol defaults
# ---------------------------------------------------------------------------

class TestAnnealProtocol:

    def test_defaults(self):
        ap = AnnealProtocol()
        assert ap.T_high == pytest.approx(600.0)
        assert ap.T_low == pytest.approx(300.0)
        assert ap.P_ref == pytest.approx(1.0)
        assert ap.em_steps == 50000
        assert ap.em_tol == pytest.approx(1000.0)
        assert ap.nvt_high_nsteps == 100000
        assert ap.npt_high_nsteps == 200000
        assert ap.anneal_nsteps == 500000
        assert ap.npt_low_nsteps == 500000
        assert ap.dt == pytest.approx(0.001)
        assert ap.tau_t == pytest.approx(0.1)
        assert ap.tau_p == pytest.approx(2.0)
        assert ap.nstxout_compressed == 5000
        assert ap.nstenergy == 1000
        assert ap.gen_seed is None


# ---------------------------------------------------------------------------
# SystemModel
# ---------------------------------------------------------------------------

class TestSystemModel:

    def test_list_fields_default_to_empty(self):
        sm = SystemModel(
            title="test.udf", udf_path="/tmp/test.udf",
            comb_rule=2, flags14=1, fudgeLJ=0.5, fudgeQQ=0.8333, calcQQ=1,
        )
        assert sm.atom_types == []
        assert sm.mol_topologies == []
        assert sm.mol_sequence == []
        assert sm.atom_positions == []
        assert sm.cell is None
        assert sm.sim_params is None
        assert sm.ndx_data is None

    def test_list_fields_are_independent_instances(self):
        """Each SystemModel instance must get its own list objects."""
        sm1 = SystemModel(
            title="a", udf_path="a", comb_rule=2, flags14=0,
            fudgeLJ=1.0, fudgeQQ=1.0, calcQQ=0,
        )
        sm2 = SystemModel(
            title="b", udf_path="b", comb_rule=2, flags14=0,
            fudgeLJ=1.0, fudgeQQ=1.0, calcQQ=0,
        )
        sm1.atom_types.append(AtomType(name="ca", mass=12.0, sigma=0.3, epsilon=0.5))
        assert len(sm2.atom_types) == 0


# ---------------------------------------------------------------------------
# MoleculeTopology
# ---------------------------------------------------------------------------

class TestMoleculeTopology:

    def test_list_defaults(self):
        mt = MoleculeTopology(udf_name="Water", gro_name="M0000", nrexcl=3)
        assert mt.atoms == []
        assert mt.bonds == []
        assert mt.pairs == []
        assert mt.angles == []
        assert mt.dihedrals == []
