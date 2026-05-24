# -*- coding: utf-8 -*-
"""End-to-end smoke test for abmptools.formulation.

Requires tleap + acpype + packmol + gmx on PATH. Skipped otherwise.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

REQUIRED_TOOLS = ("tleap", "acpype", "packmol", "gmx")
MISSING = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        bool(MISSING),
        reason=f"missing external tools: {MISSING}",
    ),
]


SAMPLE_CONFIG = (
    Path(__file__).resolve().parents[2]
    / "sample" / "formulation" / "kggggg_smoke" / "config.json"
)


def test_kggggg_smoke_end_to_end(tmp_path):
    """Build kggggg_smoke and confirm gmx grompp -f em.mdp passes."""
    out = tmp_path / "kggggg"
    result = subprocess.run(
        [
            "python", "-m", "abmptools.formulation", "build",
            "--config", str(SAMPLE_CONFIG),
            "--output-dir", str(out),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"build failed: {result.stderr}"
    assert (out / "system.top").is_file()
    assert (out / "system.gro").is_file()
    assert (out / "system.ndx").is_file()
    assert (out / "run.sh").is_file()
    for mdp in ("em", "nvt", "npt", "prod"):
        assert (out / "md" / f"{mdp}.mdp").is_file()

    # grompp em.mdp
    em_tpr = tmp_path / "em.tpr"
    result = subprocess.run(
        [
            "gmx", "grompp",
            "-f", str(out / "md" / "em.mdp"),
            "-p", str(out / "system.top"),
            "-c", str(out / "system.gro"),
            "-n", str(out / "system.ndx"),
            "-o", str(em_tpr), "-maxwarn", "5",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"grompp failed: {result.stderr}"
    assert em_tpr.is_file()


def test_kggggg_smoke_atom_count_reasonable(tmp_path):
    """Build kggggg_smoke and assert the atom count is in the expected
    range for the 6 nm box (peptides + 16 caprate + 2 taurocholate +
    TIP3P + ions ≈ 15-20k atoms)."""
    out = tmp_path / "kggggg"
    subprocess.run(
        [
            "python", "-m", "abmptools.formulation", "build",
            "--config", str(SAMPLE_CONFIG),
            "--output-dir", str(out),
        ],
        check=True, capture_output=True, text=True,
    )
    gro = (out / "system.gro").read_text().splitlines()
    n_atoms = int(gro[1].strip().split()[0])
    assert 5000 < n_atoms < 25000, f"unexpected atom count {n_atoms}"
