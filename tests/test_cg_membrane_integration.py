# -*- coding: utf-8 -*-
"""End-to-end integration tests for abmptools.cg.membrane.

Skipped unless ``insane`` / ``martinize2`` / ``gmx`` are on PATH AND the 4
required Martini 3 ITPs are present in ``./ff``. tleap is optional
(extended-backbone fallback through cg.peptide).

Marked ``slow`` -- run explicitly with::

    pytest tests/test_cg_membrane_integration.py -m slow -v
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from abmptools.cg.membrane.builder import MembraneCGBuilder
from abmptools.cg.membrane.forcefield_check import REQUIRED_MARTINI_FILES
from abmptools.cg.membrane.models import (
    EquilibrationCGProtocol,
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
    PullingCGProtocol,
    UmbrellaCGProtocol,
)


def _ff_dir_or_skip() -> Path:
    candidate = Path("./ff").resolve()
    missing = [
        f for f in REQUIRED_MARTINI_FILES if not (candidate / f).exists()
    ]
    if missing:
        pytest.skip(
            f"Martini 3 force field files missing in {candidate}: {missing}. "
            "Download martini_v300.zip from https://cgmartini.nl/."
        )
    return candidate


def _tools_or_skip() -> None:
    for tool in ("insane", "martinize2", "gmx"):
        if not shutil.which(tool):
            pytest.skip(
                f"{tool} not on PATH (install via pip / mamba or check env)."
            )


def _short_protocol_cfg(tmp_path, ff_dir):
    """A compact config: 100-step MD durations, KGG x1, POPC small bilayer."""
    return MembraneCGBuildConfig(
        lipids=[LipidMix(resname="POPC", n_per_leaflet=32)],
        peptide=PeptideMembraneSpec(name="kgg", sequence="KGG"),
        insane_d_nm=6.0,
        box_z_nm=12.0,
        martini_itp_dir=str(ff_dir),
        output_dir=str(tmp_path),
        equilibration=EquilibrationCGProtocol(
            em_steps=100, nvt_nsteps=100, npt_nsteps=100,
        ),
        pulling=PullingCGProtocol(nsteps=100),
        umbrella=UmbrellaCGProtocol(window_nsteps=100),
        grompp_maxwarn=10,
    )


@pytest.mark.slow
def test_kgg_popc_build_grompp_em(tmp_path):
    """KGG x1 + POPC -> MembraneCGBuilder.build() -> gmx grompp em.mdp PASS."""
    _tools_or_skip()
    ff_dir = _ff_dir_or_skip()

    cfg = _short_protocol_cfg(tmp_path, ff_dir)
    result = MembraneCGBuilder(cfg).build()

    # Output structure exists.
    assert result["gro"].exists()
    assert result["top"].exists()
    assert result["ndx"].exists()
    em_mdp = tmp_path / "mdp" / "em.mdp"
    assert em_mdp.exists()

    proc = subprocess.run(
        [
            "gmx", "grompp",
            "-f", str(em_mdp),
            "-c", str(result["gro"]),
            "-p", str(result["top"]),
            "-n", str(result["ndx"]),
            "-o", str(tmp_path / "em.tpr"),
            "-maxwarn", "10",
        ],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert proc.returncode == 0, (
        f"gmx grompp em.mdp failed:\nstdout:{proc.stdout}\nstderr:{proc.stderr}"
    )
    assert (tmp_path / "em.tpr").exists()


@pytest.mark.slow
def test_kgg_popc_build_grompp_pull(tmp_path):
    """Pull MDP grompps (NVT-chassis + direction-periodic + pull block)."""
    _tools_or_skip()
    ff_dir = _ff_dir_or_skip()

    cfg = _short_protocol_cfg(tmp_path, ff_dir)
    result = MembraneCGBuilder(cfg).build()

    pull_mdp = tmp_path / "pull" / "pull.mdp"
    assert pull_mdp.exists()

    proc = subprocess.run(
        [
            "gmx", "grompp",
            "-f", str(pull_mdp),
            "-c", str(result["gro"]),
            "-p", str(result["top"]),
            "-n", str(result["ndx"]),
            "-o", str(tmp_path / "pull.tpr"),
            "-maxwarn", "10",
        ],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert proc.returncode == 0, (
        f"gmx grompp pull.mdp failed:\nstdout:{proc.stdout}\nstderr:{proc.stderr}"
    )
    assert (tmp_path / "pull.tpr").exists()


@pytest.mark.slow
def test_kgg_popc_build_grompp_central_window(tmp_path):
    """Window 006 (z=0) MDP grompps (NPT-chassis + direction + static pull)."""
    _tools_or_skip()
    ff_dir = _ff_dir_or_skip()

    cfg = _short_protocol_cfg(tmp_path, ff_dir)
    result = MembraneCGBuilder(cfg).build()

    win6 = tmp_path / "windows" / "win_006" / "window.mdp"
    assert win6.exists()

    proc = subprocess.run(
        [
            "gmx", "grompp",
            "-f", str(win6),
            "-c", str(result["gro"]),
            "-p", str(result["top"]),
            "-n", str(result["ndx"]),
            "-o", str(tmp_path / "windows" / "win_006" / "window.tpr"),
            "-maxwarn", "10",
        ],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert proc.returncode == 0, (
        f"gmx grompp window.mdp failed:\nstdout:{proc.stdout}\nstderr:{proc.stderr}"
    )
    assert (tmp_path / "windows" / "win_006" / "window.tpr").exists()


@pytest.mark.slow
def test_kgg_popc_build_run_script_executable(tmp_path):
    """run.sh exists, is executable, and contains expected stages."""
    _tools_or_skip()
    ff_dir = _ff_dir_or_skip()

    cfg = _short_protocol_cfg(tmp_path, ff_dir)
    result = MembraneCGBuilder(cfg).build()

    run_sh = tmp_path / "run.sh"
    assert run_sh.exists()
    assert run_sh.stat().st_mode & 0o100, "run.sh not executable"
    text = run_sh.read_text()
    for marker in (
        "Stage 1: energy minimisation",
        "Stage 2: NVT",
        "Stage 3: NPT",
        "Stage 4: pulling",
        "Stage 5: extract per-window",
        "Stage 6: per-window MD",
        "Stage 7: PMF",
        "python -m abmptools.cg.membrane make-windows",
        "python -m abmptools.cg.membrane wham",
    ):
        assert marker in text, f"missing in run.sh: {marker!r}"
