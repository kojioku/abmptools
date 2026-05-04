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

    cfg = PeptideBuildConfig(
        peptides=[PeptideSpec(name="kgg", sequence="KGG", count=1)],
        box_size_nm=4.0,
        martini_itp_dir=str(ff_dir),
        output_dir=str(tmp_path),
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
