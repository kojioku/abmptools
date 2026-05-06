# -*- coding: utf-8 -*-
"""Integration tests for abmptools.genesis.mmgbsa (slow, gated).

Auto-skip when ``atdyn`` / ``tleap`` / ``acpype`` / ``mpirun`` are
missing on the PATH. Run explicitly with::

    pytest -m slow tests/test_mmgbsa_integration.py -v
"""
from __future__ import annotations

import shutil

import pytest

REQUIRED_TOOLS = ["atdyn", "tleap", "acpype", "mpirun"]


def _missing_tools() -> list:
    return [t for t in REQUIRED_TOOLS if shutil.which(t) is None]


pytestmark = pytest.mark.slow


@pytest.mark.skipif(
    bool(_missing_tools()),
    reason=f"Missing required tools: {_missing_tools()}",
)
class TestMMGBSASmokeIntegration:
    def test_pipeline_runs(self, tmp_path):
        """Single-target trivial smoke: tleap + acpype + atdyn must run.

        Uses a synthetic 1-residue ligand. ΔG_bind is checked only for
        plausibility, not exactness, since the fake ligand is not
        physically meaningful.
        """
        from abmptools.genesis.mmgbsa.builder import MMGBSAOrchestrator
        from abmptools.genesis.mmgbsa.models import (
            MMGBSABuildConfig,
            TargetSpec,
        )

        # Synthetic protein-ligand complex (not chemically meaningful,
        # only verifies that the toolchain plumbing works end-to-end).
        pdb_text = """\
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C
ATOM      3  C   GLY A   1       2.000   1.400   0.000  1.00 20.00           C
ATOM      4  O   GLY A   1       1.300   2.400   0.000  1.00 20.00           O
HETATM    5  C1  LIG A 201       8.000   8.000   0.000  1.00 20.00           C
HETATM    6  C2  LIG A 201       9.000   8.500   0.000  1.00 20.00           C
TER
END
"""
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        (input_dir / "smoke.pdb").write_text(pdb_text)

        cfg = MMGBSABuildConfig(
            targets=[TargetSpec(pdb="smoke.pdb", ligand_resno=201)],
            input_dir=str(input_dir),
            output_dir=str(tmp_path / "out"),
            project_name="smoke_integration",
        )
        orch = MMGBSAOrchestrator(cfg)
        result = orch.run()

        # Pipeline finished and emitted analysis outputs.
        assert result["n_succeeded"] >= 1
        assert result["csv_path"] is not None
        assert result["png_path"] is not None
        # Plausible range check (real LIG could be anything; we
        # accept a wide window).
        for t in result["targets"]:
            assert -10000.0 < t["dg_bind"] < 10000.0
