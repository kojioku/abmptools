# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.umbrella_release."""
from __future__ import annotations

from pathlib import Path

import pytest

from abmptools.formulation.models import (
    BileSaltSpec,
    EnhancerSpec,
    FormulationBuildConfig,
    PeptideSpec,
    SystemSpec,
    USProtocol,
)
from abmptools.formulation.umbrella_release import (
    render_pull_mdp,
    render_window_mdp,
    run_release_us,
    write_release_us_layout,
)


def _make_cfg(output_dir: Path) -> FormulationBuildConfig:
    return FormulationBuildConfig(
        system=SystemSpec(
            peptides=[PeptideSpec(name="p", sequence="AAA", n_copies=2)],
            enhancers=[EnhancerSpec(
                name="caprate", resname="CPR",
                smiles_neutral="C", n_neutral=2,
            )],
            bile_salts=[BileSaltSpec(
                name="tch", resname="TCH",
                smiles="O", n_copies=1,
            )],
        ),
        release_us=USProtocol(n_windows=3, window_spacing_nm=0.5,
                              force_constant_kj_mol_nm2=500.0,
                              window_nsteps=1000),
        output_dir=str(output_dir),
    )


def test_render_pull_mdp_contains_umbrella_section():
    us = USProtocol(pull_nsteps=10000, force_constant_kj_mol_nm2=1000.0)
    text = render_pull_mdp(us)
    assert "pull                     = yes" in text
    assert "pull-coord1-type         = umbrella" in text
    assert "pull-group1-name         = Enhancer" in text
    assert "pull-group2-name         = Peptide" in text
    assert "pull-coord1-k            = 1000.0" in text


def test_render_window_mdp_uses_target_distance():
    us = USProtocol(window_nsteps=500000, force_constant_kj_mol_nm2=1000.0)
    text = render_window_mdp(us, target_dist_nm=1.5)
    assert "pull-coord1-init         = 1.5" in text
    assert "pull-coord1-k            = 1000.0" in text


def test_write_release_us_layout_emits_pull_and_windows(tmp_path):
    cfg = _make_cfg(tmp_path)
    art = write_release_us_layout(config=cfg, out_dir=str(tmp_path),
                                  initial_distance_nm=0.0)
    assert Path(art.pull_mdp).is_file()
    assert len(art.window_mdps) == 3
    for w in art.window_mdps:
        assert Path(w).is_file()
    # targets: 0.0, 0.5, 1.0
    assert art.window_targets_nm == [0.0, 0.5, 1.0]


def test_write_release_us_layout_requires_us_protocol(tmp_path):
    cfg = _make_cfg(tmp_path)
    cfg.release_us = None
    with pytest.raises(ValueError, match="release_us is None"):
        write_release_us_layout(config=cfg, out_dir=str(tmp_path))


def test_run_release_us_populates_default_us(tmp_path):
    cfg = _make_cfg(tmp_path)
    cfg.release_us = None
    art = run_release_us(config=cfg, target_peptide_idx=1)
    assert cfg.release_us is not None
    assert cfg.release_us.target_peptide_idx == 1
    assert art.n_windows == cfg.release_us.n_windows
