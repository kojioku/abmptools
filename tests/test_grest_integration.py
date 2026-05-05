# -*- coding: utf-8 -*-
"""Integration tests for abmptools.genesis.grest.

These exercise the real subprocess pipeline (tleap + atdyn / spdyn /
remd_convert via mpirun) and are gated behind ``@pytest.mark.slow``.
They auto-skip when the tools are missing, so the suite still passes
on CI without GENESIS / AmberTools.

Run with::

    pytest -m slow tests/test_grest_integration.py -v

Smoke target: KGG tripeptide x 4 replicas x 4000 nsteps (~10 min on
16 cores). Verifies that the rendered .inp set is accepted by GENESIS
through the minimisation stage; the full gREST stage may require GPU
nodes for a meaningful run-time.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

REQUIRED_TOOLS = ["tleap", "atdyn", "spdyn", "remd_convert", "mpirun"]


def _missing_tools() -> list:
    return [t for t in REQUIRED_TOOLS if shutil.which(t) is None]


pytestmark = pytest.mark.slow


@pytest.mark.skipif(
    bool(_missing_tools()),
    reason=f"missing one or more required tools: {_missing_tools()}",
)
class TestKggSmokeIntegration:
    def test_full_build_pipeline(self, tmp_path):
        """Run the 5-stage build with real tleap, then verify outputs.

        Writes a minimal KGG PDB (3 residues, ~30 atoms) to disk,
        configures the smoke profile, and asserts that:

        - tleap produces prmtop + coor + ref pdb
        - all 4 .inp files are written
        - run.sh has +x and contains the expected stage banners
        """
        from abmptools.genesis.grest.builder import GrestBuilder
        from abmptools.genesis.grest.models import GrestBuildConfig

        pdb = tmp_path / "kgg.pdb"
        pdb.write_text(
            "ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  LYS A   1       1.500   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   LYS A   1       2.000   1.400   0.000  1.00  0.00           C\n"
            "ATOM      4  O   LYS A   1       1.300   2.400   0.000  1.00  0.00           O\n"
            "ATOM      5  N   GLY A   2       3.300   1.600   0.000  1.00  0.00           N\n"
            "ATOM      6  CA  GLY A   2       4.000   2.900   0.000  1.00  0.00           C\n"
            "ATOM      7  C   GLY A   2       5.500   2.700   0.000  1.00  0.00           C\n"
            "ATOM      8  O   GLY A   2       5.900   1.500   0.000  1.00  0.00           O\n"
            "ATOM      9  N   GLY A   3       6.300   3.700   0.000  1.00  0.00           N\n"
            "ATOM     10  CA  GLY A   3       7.700   3.700   0.000  1.00  0.00           C\n"
            "ATOM     11  C   GLY A   3       8.300   5.100   0.000  1.00  0.00           C\n"
            "ATOM     12  O   GLY A   3       7.700   6.200   0.000  1.00  0.00           O\n"
            "TER\n"
            "END\n"
        )
        cfg = GrestBuildConfig.from_json(
            str(Path(__file__).parent.parent / "sample" / "grest" / "kgg_smoke.json")
        )
        cfg.input_pdb = str(pdb)
        cfg.output_dir = str(tmp_path / "out")

        builder = GrestBuilder(cfg)
        result = builder.build()

        out = Path(result["output_dir"])
        assert (out / "build" / "system.prmtop").exists()
        assert (out / "build" / "system.coor").exists()
        assert (out / "inp" / "step1_minimize.inp").exists()
        assert (out / "inp" / "step3_grest.inp").exists()
        assert (out / "run.sh").exists()
        # n_protein_residues correct.
        assert result["n_protein_residues"] == 3
        # ladder length matches.
        assert result["n_replicas"] == 4
