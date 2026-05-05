# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.models dataclasses."""
from __future__ import annotations

import pytest

from abmptools.genesis.grest.models import (
    EquilibrationStage,
    GrestBuildConfig,
    GrestStage,
    MinimizationStage,
    ReplicaTemperatureSpec,
    RESTSelectionSpec,
    VALID_PARAM_TYPES,
)


# ---------------------------------------------------------------------------
# RESTSelectionSpec
# ---------------------------------------------------------------------------

class TestRESTSelectionSpec:
    def test_defaults_explicit_empty_raises(self):
        # Default mode is 'explicit' but residues is empty -> error.
        with pytest.raises(ValueError, match="non-empty residues"):
            RESTSelectionSpec()

    def test_explicit_valid(self):
        spec = RESTSelectionSpec(mode="explicit", residues=["1-138"])
        assert spec.mode == "explicit"
        assert spec.residues == ["1-138"]
        assert spec.param_types == ["C", "L"]

    def test_explicit_integer_residues_normalised_to_str(self):
        spec = RESTSelectionSpec(mode="explicit", residues=[21, 96, 274])
        assert spec.residues == ["21", "96", "274"]

    def test_around_valid(self):
        spec = RESTSelectionSpec(mode="around", center="rno:96", radius_A=4.5)
        assert spec.mode == "around"
        assert spec.center == "rno:96"
        assert spec.radius_A == 4.5

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError, match="mode must be"):
            RESTSelectionSpec(mode="bogus", residues=["1"])

    def test_around_empty_center_raises(self):
        with pytest.raises(ValueError, match="non-empty center"):
            RESTSelectionSpec(mode="around", center="", radius_A=5.0)

    def test_around_zero_radius_raises(self):
        with pytest.raises(ValueError, match="radius_A must be > 0"):
            RESTSelectionSpec(mode="around", center="rno:96", radius_A=0)

    def test_around_negative_radius_raises(self):
        with pytest.raises(ValueError, match="radius_A must be > 0"):
            RESTSelectionSpec(mode="around", center="rno:96", radius_A=-1)

    def test_unknown_param_type_raises(self):
        with pytest.raises(ValueError, match="Unknown param_type"):
            RESTSelectionSpec(
                mode="explicit", residues=["1"], param_types=["D", "Q"]
            )

    def test_all_valid_param_types(self):
        # No exception when all VALID_PARAM_TYPES are listed.
        spec = RESTSelectionSpec(
            mode="explicit",
            residues=["1"],
            param_types=sorted(VALID_PARAM_TYPES),
        )
        assert set(spec.param_types) == VALID_PARAM_TYPES


# ---------------------------------------------------------------------------
# ReplicaTemperatureSpec
# ---------------------------------------------------------------------------

class TestReplicaTemperatureSpec:
    def test_defaults_auto_geometric(self):
        rt = ReplicaTemperatureSpec()
        assert rt.mode == "auto"
        assert rt.method == "geometric"
        assert rt.T_min == 300.0
        assert rt.T_max == 357.10
        assert rt.n_replicas == 4

    def test_manual_valid(self):
        rt = ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.0, 318.11, 337.11, 357.10],
        )
        assert rt.mode == "manual"
        assert rt.temperatures[0] == 300.0
        assert rt.temperatures[-1] == 357.10

    def test_manual_non_monotonic_raises(self):
        with pytest.raises(ValueError, match="monotonically increasing"):
            ReplicaTemperatureSpec(
                mode="manual", temperatures=[300.0, 320.0, 310.0]
            )

    def test_manual_too_few_raises(self):
        with pytest.raises(ValueError, match=">= 2 temperatures"):
            ReplicaTemperatureSpec(mode="manual", temperatures=[300.0])

    def test_manual_negative_temp_raises(self):
        with pytest.raises(ValueError, match="must be > 0"):
            ReplicaTemperatureSpec(
                mode="manual", temperatures=[-1.0, 300.0]
            )

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError, match="mode must be"):
            ReplicaTemperatureSpec(mode="weird")

    def test_auto_n_too_small_raises(self):
        with pytest.raises(ValueError, match="n_replicas must be >= 2"):
            ReplicaTemperatureSpec(mode="auto", n_replicas=1)

    def test_auto_T_max_le_T_min_raises(self):
        with pytest.raises(ValueError, match="T_max"):
            ReplicaTemperatureSpec(
                mode="auto", T_min=350.0, T_max=300.0, n_replicas=4
            )

    def test_auto_T_min_zero_raises(self):
        with pytest.raises(ValueError, match="T_min must be > 0"):
            ReplicaTemperatureSpec(
                mode="auto", T_min=0.0, T_max=357.0, n_replicas=4
            )

    def test_auto_unknown_method_raises(self):
        with pytest.raises(ValueError, match="method must be"):
            ReplicaTemperatureSpec(mode="auto", method="quadratic")


# ---------------------------------------------------------------------------
# MinimizationStage
# ---------------------------------------------------------------------------

class TestMinimizationStage:
    def test_defaults(self):
        mini = MinimizationStage()
        assert mini.method == "SD"
        assert mini.nsteps == 10_000
        assert mini.cutoffdist_A == 12.0

    def test_invalid_method_raises(self):
        with pytest.raises(ValueError, match="method must be"):
            MinimizationStage(method="CG")

    def test_negative_nsteps_raises(self):
        with pytest.raises(ValueError, match="nsteps must be >= 0"):
            MinimizationStage(nsteps=-1)

    def test_negative_posres_raises(self):
        with pytest.raises(ValueError, match="posres_force"):
            MinimizationStage(posres_force_kcal_per_mol=-0.5)

    def test_pairlist_lt_cutoff_raises(self):
        with pytest.raises(ValueError, match="pairlistdist_A"):
            MinimizationStage(cutoffdist_A=12.0, pairlistdist_A=10.0)


# ---------------------------------------------------------------------------
# EquilibrationStage
# ---------------------------------------------------------------------------

class TestEquilibrationStage:
    def test_defaults(self):
        eq = EquilibrationStage()
        assert eq.integrator == "VVER"
        assert eq.timestep_ps == 0.002
        assert eq.ensemble == "NPT"
        assert eq.hydrogen_mr is True

    def test_invalid_integrator_raises(self):
        with pytest.raises(ValueError, match="integrator"):
            EquilibrationStage(integrator="DUMMY")

    def test_invalid_ensemble_raises(self):
        with pytest.raises(ValueError, match="ensemble"):
            EquilibrationStage(ensemble="NVE2")

    def test_zero_timestep_raises(self):
        with pytest.raises(ValueError, match="timestep_ps"):
            EquilibrationStage(timestep_ps=0.0)

    def test_temperature_zero_raises(self):
        with pytest.raises(ValueError, match="temperature_K"):
            EquilibrationStage(temperature_K=0.0)

    def test_pairlist_lt_cutoff_raises(self):
        with pytest.raises(ValueError, match="pairlistdist_A"):
            EquilibrationStage(cutoffdist_A=8.0, pairlistdist_A=7.0)


# ---------------------------------------------------------------------------
# GrestStage
# ---------------------------------------------------------------------------

class TestGrestStage:
    def test_defaults(self):
        grest = GrestStage()
        assert grest.integrator == "VRES"
        assert grest.timestep_ps == 0.0035
        assert grest.exchange_period == 3000
        assert grest.nsteps == 3_000_000
        assert grest.analysis_grest is True

    def test_nsteps_not_multiple_of_2_exchange_raises(self):
        # 3_000_000 / (2 * 3000) = 500, OK. Try 3_000_001.
        with pytest.raises(ValueError, match="2\\*exchange_period"):
            GrestStage(nsteps=3_000_001, exchange_period=3000)

    def test_nsteps_not_multiple_of_rstout_raises(self):
        # 3_000_000 % 30_001 != 0
        with pytest.raises(ValueError, match="rstout_period"):
            GrestStage(nsteps=3_000_000, exchange_period=3000, rstout_period=30_001)

    def test_zero_exchange_period_raises(self):
        with pytest.raises(ValueError, match="exchange_period"):
            GrestStage(exchange_period=0)

    def test_zero_rstout_raises(self):
        with pytest.raises(ValueError, match="rstout_period"):
            GrestStage(rstout_period=0)

    def test_zero_nsteps_raises(self):
        with pytest.raises(ValueError, match="nsteps"):
            GrestStage(nsteps=0)

    def test_invalid_integrator_raises(self):
        with pytest.raises(ValueError, match="integrator"):
            GrestStage(integrator="JUNK")

    def test_invalid_ensemble_raises(self):
        with pytest.raises(ValueError, match="ensemble"):
            GrestStage(ensemble="NVT2")

    def test_pairlist_lt_cutoff_raises(self):
        with pytest.raises(ValueError, match="pairlistdist_A"):
            GrestStage(cutoffdist_A=8.0, pairlistdist_A=7.0)


# ---------------------------------------------------------------------------
# GrestBuildConfig
# ---------------------------------------------------------------------------

def _valid_cfg_kwargs(**overrides):
    """Build a minimal valid kwargs dict for GrestBuildConfig."""
    base = dict(
        input_pdb="/tmp/protein.pdb",
        rest_selection=RESTSelectionSpec(
            mode="explicit", residues=["1-3"]
        ),
        replica_temperatures=ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.0, 318.11, 337.11, 357.10],
        ),
    )
    base.update(overrides)
    return base


class TestGrestBuildConfig:
    def test_defaults_with_required(self):
        cfg = GrestBuildConfig(**_valid_cfg_kwargs())
        assert cfg.input_pdb == "/tmp/protein.pdb"
        assert cfg.project_name == "grest_run"
        assert cfg.ff_protein == "leaprc.protein.ff19SB"
        assert cfg.ff_water == "leaprc.water.tip3p"
        assert cfg.solvatebox_padding_A == 10.0
        assert cfg.mpi_processes_per_replica == 2

    def test_missing_input_pdb_raises(self):
        with pytest.raises(ValueError, match="input_pdb is required"):
            GrestBuildConfig(**_valid_cfg_kwargs(input_pdb=""))

    def test_empty_project_name_raises(self):
        with pytest.raises(ValueError, match="project_name"):
            GrestBuildConfig(**_valid_cfg_kwargs(project_name=""))

    def test_zero_padding_raises(self):
        with pytest.raises(ValueError, match="solvatebox_padding_A"):
            GrestBuildConfig(**_valid_cfg_kwargs(solvatebox_padding_A=0))

    def test_negative_salt_raises(self):
        with pytest.raises(ValueError, match="salt_concentration_M"):
            GrestBuildConfig(**_valid_cfg_kwargs(salt_concentration_M=-0.1))

    def test_zero_mpi_raises(self):
        with pytest.raises(ValueError, match="mpi_processes_per_replica"):
            GrestBuildConfig(**_valid_cfg_kwargs(mpi_processes_per_replica=0))

    def test_zero_omp_raises(self):
        with pytest.raises(ValueError, match="omp_num_threads"):
            GrestBuildConfig(**_valid_cfg_kwargs(omp_num_threads=0))

    def test_lowest_T_mismatch_raises(self):
        # grest stage default temperature_K == 300.0; ladder begins at 310.0.
        with pytest.raises(ValueError, match="match the lowest replica"):
            GrestBuildConfig(
                **_valid_cfg_kwargs(
                    replica_temperatures=ReplicaTemperatureSpec(
                        mode="manual",
                        temperatures=[310.0, 320.0],
                    ),
                )
            )

    def test_lowest_T_match_auto_mode(self):
        cfg = GrestBuildConfig(
            **_valid_cfg_kwargs(
                replica_temperatures=ReplicaTemperatureSpec(
                    mode="auto", T_min=300.0, T_max=357.10, n_replicas=4
                ),
            )
        )
        assert cfg.replica_temperatures.mode == "auto"

    def test_to_json_roundtrip(self, tmp_path):
        cfg = GrestBuildConfig(**_valid_cfg_kwargs())
        out = tmp_path / "cfg.json"
        cfg.to_json(str(out))
        loaded = GrestBuildConfig.from_json(str(out))
        assert loaded.input_pdb == cfg.input_pdb
        assert loaded.project_name == cfg.project_name
        assert loaded.rest_selection.mode == "explicit"
        assert loaded.rest_selection.residues == ["1-3"]
        assert loaded.replica_temperatures.mode == "manual"
        assert loaded.replica_temperatures.temperatures == [
            300.0, 318.11, 337.11, 357.10,
        ]
        assert loaded.minimize.nsteps == cfg.minimize.nsteps
        assert loaded.equilibrate.timestep_ps == cfg.equilibrate.timestep_ps
        assert loaded.grest.exchange_period == cfg.grest.exchange_period

    def test_to_json_roundtrip_around_mode(self, tmp_path):
        cfg = GrestBuildConfig(
            **_valid_cfg_kwargs(
                rest_selection=RESTSelectionSpec(
                    mode="around", center="rno:96", radius_A=5.0
                ),
            )
        )
        out = tmp_path / "cfg.json"
        cfg.to_json(str(out))
        loaded = GrestBuildConfig.from_json(str(out))
        assert loaded.rest_selection.mode == "around"
        assert loaded.rest_selection.center == "rno:96"
        assert loaded.rest_selection.radius_A == 5.0
