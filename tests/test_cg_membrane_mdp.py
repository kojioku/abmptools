# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.mdp_templates."""
from __future__ import annotations

from abmptools.cg.membrane import mdp_templates
from abmptools.cg.membrane.models import (
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
)


def _cfg():
    return MembraneCGBuildConfig(
        lipids=[LipidMix()],
        peptide=PeptideMembraneSpec(),
        martini_itp_dir="./ff",
    )


# ---------------------------------------------------------------------------
# Equilibration MDPs
# ---------------------------------------------------------------------------

class TestEmMdp:
    def test_em_uses_steep(self):
        text = mdp_templates.render_em_mdp(config=_cfg())
        assert "integrator" in text and "steep" in text

    def test_em_has_reaction_field(self):
        text = mdp_templates.render_em_mdp(config=_cfg())
        assert "reaction-field" in text
        assert "rcoulomb                 = 1.1" in text


class TestNvtMdp:
    def test_nvt_two_group_thermostat(self):
        text = mdp_templates.render_nvt_mdp(config=_cfg())
        assert "tc-grps" in text and "Bilayer Non_Bilayer" in text

    def test_nvt_no_pcoupl(self):
        text = mdp_templates.render_nvt_mdp(config=_cfg())
        assert "pcoupl                   = no" in text

    def test_nvt_dt_20fs(self):
        text = mdp_templates.render_nvt_mdp(config=_cfg())
        assert "dt                       = 0.0200" in text

    def test_nvt_no_constraints(self):
        text = mdp_templates.render_nvt_mdp(config=_cfg())
        assert "constraints              = none" in text


class TestNptMdp:
    def test_npt_semiisotropic(self):
        text = mdp_templates.render_npt_mdp(config=_cfg())
        assert "pcoupltype               = semiisotropic" in text

    def test_npt_compressibility_3em4(self):
        text = mdp_templates.render_npt_mdp(config=_cfg())
        assert "compressibility          = 3e-4 3e-4" in text

    def test_npt_tau_p_12(self):
        text = mdp_templates.render_npt_mdp(config=_cfg())
        assert "tau_p                    = 12.00" in text

    def test_npt_two_group_thermostat(self):
        text = mdp_templates.render_npt_mdp(config=_cfg())
        assert "Bilayer Non_Bilayer" in text


# ---------------------------------------------------------------------------
# Pull / window MDPs
# ---------------------------------------------------------------------------

class TestPullMdp:
    def test_pull_uses_direction_periodic(self):
        text = mdp_templates.render_pull_mdp(
            config=_cfg(), pull_init_nm=3.0,
        )
        assert "direction-periodic" in text

    def test_pull_no_pcoupl(self):
        """Pulling requires no pressure coupling (NVT chassis)."""
        text = mdp_templates.render_pull_mdp(
            config=_cfg(), pull_init_nm=3.0,
        )
        assert "pcoupl                   = no" in text

    def test_pull_rate_nonzero(self):
        text = mdp_templates.render_pull_mdp(
            config=_cfg(), pull_init_nm=3.0,
        )
        assert "pull-coord1-rate           = 0.001000" in text

    def test_pull_init_value(self):
        text = mdp_templates.render_pull_mdp(
            config=_cfg(), pull_init_nm=2.5,
        )
        assert "pull-coord1-init           = 2.5000" in text

    def test_pull_pbc_atom_inserted(self):
        text = mdp_templates.render_pull_mdp(
            config=_cfg(), pull_init_nm=3.0,
            pbc_atom_g1=1234,
        )
        assert "pull-group1-pbcatom        = 1234" in text

    def test_pull_groups_default(self):
        text = mdp_templates.render_pull_mdp(
            config=_cfg(), pull_init_nm=3.0,
        )
        assert "pull-group1-name           = Bilayer" in text
        assert "pull-group2-name           = Peptide" in text

    def test_pull_uses_nsteps_from_pulling(self):
        cfg = _cfg()
        cfg.pulling.nsteps = 12345
        text = mdp_templates.render_pull_mdp(
            config=cfg, pull_init_nm=3.0,
        )
        assert "nsteps                   = 12345" in text


class TestWindowMdp:
    def test_window_uses_direction_no_periodic(self):
        """Window MDP uses direction (not direction-periodic) for dynamic box."""
        text = mdp_templates.render_window_mdp(
            config=_cfg(), window_z_nm=0.0,
        )
        assert "direction" in text
        # Window allows semiisotropic Pcoupl, so direction-periodic is wrong.
        assert "direction-periodic" not in text

    def test_window_has_pcoupl(self):
        text = mdp_templates.render_window_mdp(
            config=_cfg(), window_z_nm=0.0,
        )
        assert "pcoupltype               = semiisotropic" in text

    def test_window_static_pull(self):
        text = mdp_templates.render_window_mdp(
            config=_cfg(), window_z_nm=-1.0,
        )
        assert "pull-coord1-rate           = 0.000000" in text

    def test_window_init_at_target_z(self):
        text = mdp_templates.render_window_mdp(
            config=_cfg(), window_z_nm=-0.75,
        )
        assert "pull-coord1-init           = -0.7500" in text

    def test_window_uses_force_constant_from_config(self):
        cfg = _cfg()
        cfg.umbrella.force_constant_kj_mol_nm2 = 2500.0
        text = mdp_templates.render_window_mdp(
            config=cfg, window_z_nm=0.0,
        )
        assert "pull-coord1-k              = 2500.00" in text

    def test_window_uses_window_nsteps(self):
        cfg = _cfg()
        cfg.umbrella.window_nsteps = 99999
        text = mdp_templates.render_window_mdp(
            config=cfg, window_z_nm=0.0,
        )
        assert "nsteps                   = 99999" in text


# ---------------------------------------------------------------------------
# Bulk writer
# ---------------------------------------------------------------------------

def test_write_equilibration_mdps_creates_3_files(tmp_path):
    paths = mdp_templates.write_equilibration_mdps(
        config=_cfg(), equil_dir=tmp_path,
    )
    assert "em" in paths and paths["em"].exists()
    assert "nvt" in paths and paths["nvt"].exists()
    assert "npt" in paths and paths["npt"].exists()


def test_write_equilibration_mdps_content_distinct(tmp_path):
    paths = mdp_templates.write_equilibration_mdps(
        config=_cfg(), equil_dir=tmp_path,
    )
    em = paths["em"].read_text()
    nvt = paths["nvt"].read_text()
    npt = paths["npt"].read_text()
    assert "steep" in em
    assert "v-rescale" in nvt and "pcoupl                   = no" in nvt
    assert "semiisotropic" in npt
