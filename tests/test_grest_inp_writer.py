# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.inp_writer."""
from __future__ import annotations

import pytest

from abmptools.genesis.grest.inp_writer import (
    SystemMetadata,
    render_equilibrate_inp,
    render_grest_inp,
    render_minimize_inp,
    render_remd_convert_inp,
)
from abmptools.genesis.grest.models import (
    GrestBuildConfig,
    ReplicaTemperatureSpec,
    RESTSelectionSpec,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def minimal_cfg() -> GrestBuildConfig:
    return GrestBuildConfig(
        input_pdb="/tmp/protein.pdb",
        rest_selection=RESTSelectionSpec(
            mode="explicit", residues=["21,96,274-275"]
        ),
        replica_temperatures=ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.0, 318.11, 337.11, 357.10],
        ),
    )


@pytest.fixture
def meta() -> SystemMetadata:
    return SystemMetadata(
        box_size_A=(81.056806, 83.517861, 92.905634),
        n_protein_residues=350,
        rest_selection_string="rno:21,96,274-275",
    )


@pytest.fixture
def poc_ladder() -> list:
    return [300.0, 318.11, 337.11, 357.10]


# ---------------------------------------------------------------------------
# render_minimize_inp
# ---------------------------------------------------------------------------

class TestRenderMinimizeInp:
    def test_basic_structure(self, minimal_cfg, meta):
        text = render_minimize_inp(minimal_cfg, meta)
        # All required GENESIS sections present.
        for section in [
            "[INPUT]",
            "[OUTPUT]",
            "[ENERGY]",
            "[MINIMIZE]",
            "[BOUNDARY]",
            "[SELECTION]",
            "[RESTRAINTS]",
        ]:
            assert section in text

    def test_amber_force_field(self, minimal_cfg, meta):
        text = render_minimize_inp(minimal_cfg, meta)
        assert "forcefield       = AMBER" in text
        # AMBER convention: no vdw_force_switch.
        assert "vdw_force_switch = NO" in text

    def test_box_dimensions_inserted(self, minimal_cfg, meta):
        text = render_minimize_inp(minimal_cfg, meta)
        assert "box_size_x = 81.056806" in text
        assert "box_size_y = 83.517861" in text
        assert "box_size_z = 92.905634" in text

    def test_posres_uses_n_protein_residues(self, minimal_cfg, meta):
        text = render_minimize_inp(minimal_cfg, meta)
        assert "group1 = rno:1-350 and heavy" in text

    def test_minimize_method_and_nsteps(self, minimal_cfg, meta):
        text = render_minimize_inp(minimal_cfg, meta)
        assert "method        = SD" in text
        assert "nsteps        = 10000" in text

    def test_custom_filenames(self, minimal_cfg, meta):
        text = render_minimize_inp(
            minimal_cfg, meta,
            prmtop_name="custom.prmtop",
            coor_name="init.coor",
            ref_coor_name="ref.coor",
        )
        assert "prmtopfile = custom.prmtop" in text
        assert "ambcrdfile = init.coor" in text
        assert "ambreffile = ref.coor" in text


# ---------------------------------------------------------------------------
# render_equilibrate_inp
# ---------------------------------------------------------------------------

class TestRenderEquilibrateInp:
    def test_required_sections(self, minimal_cfg, meta):
        text = render_equilibrate_inp(minimal_cfg, meta)
        for section in [
            "[INPUT]",
            "[OUTPUT]",
            "[ENERGY]",
            "[DYNAMICS]",
            "[CONSTRAINTS]",
            "[ENSEMBLE]",
            "[BOUNDARY]",
            "[SELECTION]",
            "[RESTRAINTS]",
        ]:
            assert section in text

    def test_npt_ensemble(self, minimal_cfg, meta):
        text = render_equilibrate_inp(minimal_cfg, meta)
        assert "ensemble    = NPT" in text
        assert "tpcontrol   = BUSSI" in text

    def test_hmr_yes_default(self, minimal_cfg, meta):
        text = render_equilibrate_inp(minimal_cfg, meta)
        assert "hydrogen_mr      = YES" in text
        assert "hmr_ratio        = 3.00" in text

    def test_integrator_vver(self, minimal_cfg, meta):
        text = render_equilibrate_inp(minimal_cfg, meta)
        assert "integrator       = VVER" in text

    def test_temperature_default_300(self, minimal_cfg, meta):
        text = render_equilibrate_inp(minimal_cfg, meta)
        assert "temperature = 300.000" in text


# ---------------------------------------------------------------------------
# render_grest_inp
# ---------------------------------------------------------------------------

class TestRenderGrestInp:
    def test_remd_block_present(self, minimal_cfg, meta, poc_ladder):
        text = render_grest_inp(minimal_cfg, meta, poc_ladder)
        assert "[REMD]" in text
        assert "type1           = REST" in text
        assert "nreplica1       = 4" in text
        assert "parameters1     = 300.000 318.110 337.110 357.100" in text

    def test_default_param_type_C_L(self, minimal_cfg, meta, poc_ladder):
        text = render_grest_inp(minimal_cfg, meta, poc_ladder)
        # POC SSCR setting: CHARGE + LJ.
        assert "param_type1     = C L" in text

    def test_custom_param_type(self, minimal_cfg, meta, poc_ladder):
        minimal_cfg.rest_selection = RESTSelectionSpec(
            mode="explicit",
            residues=["21,96"],
            param_types=["D", "L"],
        )
        text = render_grest_inp(minimal_cfg, meta, poc_ladder)
        assert "param_type1     = D L" in text

    def test_selection_uses_rest_string(self, minimal_cfg, meta, poc_ladder):
        text = render_grest_inp(minimal_cfg, meta, poc_ladder)
        assert "group1 = rno:21,96,274-275" in text

    def test_replica_count_from_ladder_length(self, minimal_cfg, meta):
        ladder_8 = [300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0]
        # Update cross-field constraint: lowest-T must match grest.temperature_K.
        text = render_grest_inp(minimal_cfg, meta, ladder_8)
        assert "nreplica1       = 8" in text

    def test_exchange_period_default_3000(self, minimal_cfg, meta, poc_ladder):
        text = render_grest_inp(minimal_cfg, meta, poc_ladder)
        assert "exchange_period = 3000" in text

    def test_integrator_vres(self, minimal_cfg, meta, poc_ladder):
        text = render_grest_inp(minimal_cfg, meta, poc_ladder)
        assert "integrator         = VRES" in text

    def test_replica_per_file_template(self, minimal_cfg, meta, poc_ladder):
        text = render_grest_inp(minimal_cfg, meta, poc_ladder)
        # GENESIS uses {} as the per-replica placeholder.
        assert "step3_rep{}.dcd" in text
        assert "step3_rep{}.rem" in text
        assert "step3_rep{}.ene" in text


# ---------------------------------------------------------------------------
# render_remd_convert_inp
# ---------------------------------------------------------------------------

class TestRenderRemdConvertInp:
    def test_basic_structure(self, minimal_cfg, meta, poc_ladder):
        text = render_remd_convert_inp(minimal_cfg, meta, poc_ladder)
        for section in [
            "[INPUT]",
            "[OUTPUT]",
            "[SELECTION]",
            "[FITTING]",
            "[OPTION]",
        ]:
            assert section in text

    def test_default_convert_ids_empty(self, minimal_cfg, meta, poc_ladder):
        # POC convention: empty -> lowest-T extraction only.
        text = render_remd_convert_inp(minimal_cfg, meta, poc_ladder)
        assert "convert_ids     = " in text

    def test_explicit_convert_ids(self, minimal_cfg, meta, poc_ladder):
        text = render_remd_convert_inp(
            minimal_cfg, meta, poc_ladder, convert_ids=[1, 2]
        )
        assert "convert_ids     = 1 2" in text

    def test_replica_count(self, minimal_cfg, meta, poc_ladder):
        text = render_remd_convert_inp(minimal_cfg, meta, poc_ladder)
        assert "num_replicas    = 4" in text

    def test_per_replica_template(self, minimal_cfg, meta, poc_ladder):
        text = render_remd_convert_inp(minimal_cfg, meta, poc_ladder)
        assert "dcdfile    = step3_rep{}.dcd" in text
        assert "trjfile = param{}.dcd" in text

    def test_convert_type_parameter(self, minimal_cfg, meta, poc_ladder):
        text = render_remd_convert_inp(minimal_cfg, meta, poc_ladder)
        assert "convert_type    = PARAMETER" in text


# ---------------------------------------------------------------------------
# Snapshot test (full POC reproduction)
# ---------------------------------------------------------------------------

def test_full_minimize_snapshot(minimal_cfg, meta):
    """Pin the exact rendering against a captured POC-derived snapshot.

    If the template format changes intentionally, regenerate by running:

        python -c "from abmptools.genesis.grest.inp_writer import \
            render_minimize_inp, SystemMetadata; ..."
    """
    text = render_minimize_inp(minimal_cfg, meta)
    # Pin the load-bearing tokens (sections + key=value pairs) without
    # over-constraining whitespace inside comments.
    expected_substrings = [
        "[INPUT]",
        "prmtopfile = system.prmtop",
        "ambcrdfile = system.coor",
        "ambreffile = system.coor",
        "[OUTPUT]",
        "dcdfile = step1.dcd",
        "rstfile = step1.rst",
        "[ENERGY]",
        "forcefield       = AMBER",
        "[MINIMIZE]",
        "method        = SD",
        "nsteps        = 10000",
        "[BOUNDARY]",
        "type       = PBC",
        "box_size_x = 81.056806",
        "box_size_y = 83.517861",
        "box_size_z = 92.905634",
        "[SELECTION]",
        "group1 = rno:1-350 and heavy",
        "[RESTRAINTS]",
        "function1     = POSI",
        "constant1     = 1.0000",
        "select_index1 = 1",
    ]
    for sub in expected_substrings:
        assert sub in text, f"missing token: {sub!r}"
