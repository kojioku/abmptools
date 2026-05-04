# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.mdp_templates module."""
from __future__ import annotations

import pytest

from abmptools.cg.peptide.mdp_templates import (
    generate_mdp_files,
    render_mdp,
)
from abmptools.cg.peptide.models import PeptideBuildConfig


class TestRenderMdp:
    def test_em(self):
        cfg = PeptideBuildConfig(em_steps=12345)
        text = render_mdp("em", cfg)
        assert "integrator  = steep" in text
        assert "nsteps      = 12345" in text
        assert "rcoulomb        = 1.1" in text  # M3 cutoff

    def test_nvt_uses_dt_and_temp(self):
        cfg = PeptideBuildConfig(
            dt_fs=20.0, nvt_nsteps=100000, temperature=300.0,
        )
        text = render_mdp("nvt", cfg)
        assert "dt          = 0.0200" in text
        assert "nsteps      = 100000" in text
        assert "ref_t       = 300.0" in text
        assert "tcoupl      = v-rescale" in text
        assert "pcoupl      = no" in text

    def test_npt_has_pressure_coupling(self):
        cfg = PeptideBuildConfig(npt_nsteps=200000, temperature=310.0)
        text = render_mdp("npt", cfg)
        assert "pcoupl          = c-rescale" in text
        assert "ref_t       = 310.0" in text
        assert "nsteps      = 200000" in text

    def test_md_default(self):
        cfg = PeptideBuildConfig(md_nsteps=5_000_000)
        text = render_mdp("md", cfg)
        assert "nsteps      = 5000000" in text
        assert "Production MD" in text

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError, match="Unknown"):
            render_mdp("bogus", PeptideBuildConfig())


class TestGenerateMdpFiles:
    def test_all_four_default(self, tmp_path):
        cfg = PeptideBuildConfig()
        out = generate_mdp_files(cfg, tmp_path)
        assert set(out.keys()) == {"em", "nvt", "npt", "md"}
        for p in out.values():
            assert p.exists()

    def test_disable_md(self, tmp_path):
        cfg = PeptideBuildConfig(mdp_md=False)
        out = generate_mdp_files(cfg, tmp_path)
        assert "md" not in out
        assert "em" in out
        assert "nvt" in out
        assert "npt" in out

    def test_disable_all(self, tmp_path):
        cfg = PeptideBuildConfig(
            mdp_em=False, mdp_nvt=False, mdp_npt=False, mdp_md=False,
        )
        out = generate_mdp_files(cfg, tmp_path)
        assert out == {}
