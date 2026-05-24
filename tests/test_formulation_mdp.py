# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.mdp_templates."""
from __future__ import annotations

from pathlib import Path

import pytest

from abmptools.formulation.mdp_templates import (
    DEFAULT_TC_GROUPS,
    render_em_mdp,
    render_npt_eq_mdp,
    render_nvt_mdp,
    render_production_mdp,
    write_equilibration_mdps,
)
from abmptools.formulation.models import (
    EquilibrationProtocol,
    ProductionProtocol,
)


@pytest.fixture
def equil():
    return EquilibrationProtocol()


@pytest.fixture
def prod():
    return ProductionProtocol()


# ---------------------------------------------------------------------------
# em.mdp
# ---------------------------------------------------------------------------


def test_em_mdp_uses_steep(equil):
    text = render_em_mdp(equilibration=equil)
    assert "integrator               = steep" in text
    assert "nsteps                   = 10000" in text


def test_em_mdp_no_thermostat(equil):
    text = render_em_mdp(equilibration=equil)
    assert "tcoupl" not in text
    assert "pcoupl" not in text


def test_em_mdp_force_switch_vdw(equil):
    text = render_em_mdp(equilibration=equil)
    assert "vdw-modifier             = Force-switch" in text
    assert "rvdw-switch              = 1.0" in text
    assert "rvdw                     = 1.2" in text


# ---------------------------------------------------------------------------
# nvt.mdp
# ---------------------------------------------------------------------------


def test_nvt_mdp_two_group_thermostat_default(equil):
    text = render_nvt_mdp(equilibration=equil)
    g1, g2 = DEFAULT_TC_GROUPS
    assert f"tc-grps                  = {g1} {g2}" in text


def test_nvt_mdp_custom_tc_groups(equil):
    text = render_nvt_mdp(
        equilibration=equil, tc_groups=("Protein", "Water_and_ions")
    )
    assert "Protein Water_and_ions" in text


def test_nvt_mdp_temperature_in_ref_t(equil):
    e = EquilibrationProtocol(temperature_K=323.15)
    text = render_nvt_mdp(equilibration=e)
    assert "323.15 323.15" in text


def test_nvt_mdp_no_pressure_coupling(equil):
    text = render_nvt_mdp(equilibration=equil)
    assert "pcoupl                   = no" in text


def test_nvt_mdp_constraints_h_bonds(equil):
    text = render_nvt_mdp(equilibration=equil)
    assert "constraints              = h-bonds" in text


def test_nvt_mdp_gen_seed_default_minus_one(equil):
    text = render_nvt_mdp(equilibration=equil)
    assert "gen_seed                 = -1" in text


def test_nvt_mdp_gen_seed_explicit():
    e = EquilibrationProtocol(gen_seed=42)
    text = render_nvt_mdp(equilibration=e)
    assert "gen_seed                 = 42" in text


# ---------------------------------------------------------------------------
# npt eq.mdp
# ---------------------------------------------------------------------------


def test_npt_eq_uses_c_rescale_isotropic(equil):
    text = render_npt_eq_mdp(equilibration=equil)
    assert "pcoupl                   = c-rescale" in text
    assert "pcoupltype               = isotropic" in text


def test_npt_eq_ref_p_default(equil):
    text = render_npt_eq_mdp(equilibration=equil)
    assert "ref_p                    = 1.0" in text


# ---------------------------------------------------------------------------
# production.mdp
# ---------------------------------------------------------------------------


def test_production_uses_parrinello_rahman(prod):
    text = render_production_mdp(production=prod)
    assert "pcoupl                   = Parrinello-Rahman" in text
    assert "pcoupltype               = isotropic" in text


def test_production_default_500000_steps(prod):
    text = render_production_mdp(production=prod)
    assert "nsteps                   = 500000" in text


def test_production_dt_is_2fs(prod):
    text = render_production_mdp(production=prod)
    assert "dt                       = 0.002" in text


def test_production_pme_1_2nm_cutoff(prod):
    text = render_production_mdp(production=prod)
    assert "coulombtype              = PME" in text
    assert "rcoulomb                 = 1.2" in text


# ---------------------------------------------------------------------------
# write_equilibration_mdps
# ---------------------------------------------------------------------------


def test_write_equilibration_mdps_creates_all_files(tmp_path, equil, prod):
    paths = write_equilibration_mdps(
        equilibration=equil, production=prod, out_dir=str(tmp_path),
    )
    assert set(paths.keys()) == {"em", "nvt", "npt", "prod"}
    for p in paths.values():
        assert Path(p).is_file()


def test_write_equilibration_mdps_idempotent(tmp_path, equil, prod):
    out = str(tmp_path / "md")
    write_equilibration_mdps(
        equilibration=equil, production=prod, out_dir=out,
    )
    # Calling twice should not error (overwrite)
    paths = write_equilibration_mdps(
        equilibration=equil, production=prod, out_dir=out,
    )
    assert Path(paths["em"]).is_file()
