# -*- coding: utf-8 -*-
"""End-to-end integration test for abmptools.cg.peptide.

Skipped unless ``martinize2`` and ``gmx`` are on PATH **and** the four
required Martini 3 force field files are present in ``./ff``. tleap is
optional (extended-backbone fallback is exercised when missing).

Marked ``slow`` -- run explicitly with::

    pytest tests/test_cg_peptide_integration.py -m slow -v
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from abmptools.cg.peptide.builder import PeptideCGBuilder
from abmptools.cg.peptide.forcefield_check import REQUIRED_MARTINI_FILES
from abmptools.cg.peptide.models import PeptideBuildConfig, PeptideSpec


def _ff_dir_or_skip() -> Path:
    candidate = Path("./ff").resolve()
    missing = [
        f for f in REQUIRED_MARTINI_FILES if not (candidate / f).exists()
    ]
    if missing:
        pytest.skip(
            f"Martini 3 force field files missing in {candidate}: "
            f"{missing}. Download from https://cgmartini.nl/."
        )
    return candidate


@pytest.mark.slow
def test_kgg_build_grompp_em(tmp_path):
    """KGG (1 copy, 4 nm cubic) -> gmx grompp -f em.mdp PASS."""
    if not shutil.which("martinize2"):
        pytest.skip("martinize2 not on PATH (vermouth not installed)")
    if not shutil.which("gmx"):
        pytest.skip("gmx not on PATH (GROMACS not installed)")

    ff_dir = _ff_dir_or_skip()

    # solvent_enabled=False keeps the test self-contained:
    # cgmartini.nl does not distribute a Martini 3 water.gro directly,
    # so end-to-end gmx solvate would require a user-generated water box.
    # The build still exercises atomistic + martinize2 + gmx insert-molecules
    # + topology + MDP + grompp.
    cfg = PeptideBuildConfig(
        peptides=[PeptideSpec(name="kgg", sequence="KGG", count=1)],
        box_size_nm=4.0,
        martini_itp_dir=str(ff_dir),
        output_dir=str(tmp_path),
        solvent_enabled=False,
        em_steps=100,
        nvt_nsteps=100,
        npt_nsteps=100,
        md_nsteps=100,
    )
    result = PeptideCGBuilder(cfg).build()

    # Build outputs exist
    assert result["gro"].exists()
    assert result["top"].exists()
    em_mdp = tmp_path / "mdp" / "em.mdp"
    assert em_mdp.exists()

    # gmx grompp must succeed on the built system
    proc = subprocess.run(
        [
            "gmx", "grompp",
            "-f", str(em_mdp),
            "-c", str(result["gro"]),
            "-p", str(result["top"]),
            "-o", str(tmp_path / "em.tpr"),
            "-maxwarn", "5",
        ],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert proc.returncode == 0, (
        f"gmx grompp failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    assert (tmp_path / "em.tpr").exists()


@pytest.mark.slow
def test_kgg_build_with_solvate_and_water_autogen(tmp_path):
    """Full pipeline including ``gmx solvate`` + ``gmx genion``.

    cgmartini.nl does not distribute ``martini_v3.0.0_water.gro``; the
    builder's water_box helper should auto-generate one via
    ``gmx insert-molecules``. This test verifies the full 6-stage flow
    ends with a working ``system_ions.gro`` that ``gmx grompp`` accepts.
    """
    if not shutil.which("martinize2"):
        pytest.skip("martinize2 not on PATH (vermouth not installed)")
    if not shutil.which("gmx"):
        pytest.skip("gmx not on PATH (GROMACS not installed)")

    ff_dir = _ff_dir_or_skip()

    cfg = PeptideBuildConfig(
        peptides=[PeptideSpec(name="kgg", sequence="KGG", count=1)],
        box_size_nm=5.0,
        martini_itp_dir=str(ff_dir),
        output_dir=str(tmp_path),
        solvent_enabled=True,
        neutralize=True,
        nacl_molar=0.15,
        em_steps=100,
        nvt_nsteps=100,
        npt_nsteps=100,
        md_nsteps=100,
    )
    result = PeptideCGBuilder(cfg).build()

    # Final coordinates are post-genion (system_ions.gro)
    assert result["gro"].name == "system_ions.gro"
    assert result["gro"].exists()
    # Auto-generated water box was placed in output_dir
    auto_water = tmp_path / "_auto_martini_w_box.gro"
    assert auto_water.exists()

    # gmx grompp on em.mdp must succeed against the solvated system
    proc = subprocess.run(
        [
            "gmx", "grompp",
            "-f", str(tmp_path / "mdp" / "em.mdp"),
            "-c", str(result["gro"]),
            "-p", str(result["top"]),
            "-o", str(tmp_path / "em.tpr"),
            "-maxwarn", "5",
        ],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert proc.returncode == 0, (
        f"gmx grompp failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    assert (tmp_path / "em.tpr").exists()
