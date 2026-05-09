# -*- coding: utf-8 -*-
"""CLI tests for ``abmptools.crystal.cli``.

Covers argparse wiring, ``example`` / ``validate`` / ``nearest``
plus a YAML-config-driven ``pipeline`` smoke run on csp7 R00001 that
re-uses the byte-equivalence guarantee from
:mod:`tests.test_crystal_builder`.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys

import pytest

TESTS_DIR = os.path.dirname(__file__)
REF_MAIN = os.path.join(TESTS_DIR, "regression", "reference", "main")
CRYSTAL_FIXTURE = os.path.join(REF_MAIN, "crystal_csp7")
SAMPLE_CSP7_SMOKE = os.path.abspath(
    os.path.join(TESTS_DIR, os.pardir, "sample", "crystal", "csp7_smoke"),
)
UNK_AJF_TEMPLATE = os.path.join(SAMPLE_CSP7_SMOKE, "UNK.ajf")

_HAS_FIXTURE = os.path.isdir(os.path.join(CRYSTAL_FIXTURE, "R00001"))
_HAS_TEMPLATE = os.path.isfile(UNK_AJF_TEMPLATE)


def _run_cli(*args, cwd=None):
    return subprocess.run(
        [sys.executable, "-m", "abmptools.crystal", *args],
        cwd=cwd,
        capture_output=True, text=True, timeout=600,
    )


# ---------------------------------------------------------------------------
# example / validate / nearest (lightweight)
# ---------------------------------------------------------------------------

def test_example_yaml_is_self_loadable(tmp_path):
    """`example` output must be parseable as a CrystalBuildConfig."""
    pytest.importorskip("yaml")
    out = _run_cli("example")
    assert out.returncode == 0, out.stderr
    yaml_text = out.stdout
    assert "project_name:" in yaml_text
    cfg_path = tmp_path / "from_example.yaml"
    cfg_path.write_text(yaml_text)
    from abmptools.crystal.models import CrystalBuildConfig
    cfg = CrystalBuildConfig.from_yaml(str(cfg_path))
    assert cfg.fmo.is_xyz is True  # crystal default
    assert cfg.cif_engine.engine == "legacy"


def test_validate_no_config():
    """`validate` without --config still runs (Phase A behaviour)."""
    out = _run_cli("validate")
    assert out.returncode == 0
    assert "abmptools.crystal" in out.stdout


def test_help_shows_eight_subcommands():
    out = _run_cli("--help")
    assert out.returncode == 0
    for name in ("example", "validate", "expand", "fragment",
                 "jobs", "pipeline", "postproc", "nearest"):
        assert name in out.stdout, f"missing subcommand in help: {name}"


def test_nearest_subcommand(tmp_path):
    """``nearest`` standalone on a synthetic PDB."""
    pdb = tmp_path / "synth.pdb"
    pdb.write_text(
        "HETATM    1  C1  UNK     1      -0.500   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  UNK     1       0.500   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  N1  UNK     2       0.000   2.000   0.000  1.00  0.00           N\n"
        "HETATM    4  O1  UNK     3       0.000   0.000   5.000  1.00  0.00           O\n"
        "END\n"
    )
    out = _run_cli("nearest", "--pdb", str(pdb), "--center", "1", "-n", "2")
    assert out.returncode == 0, out.stderr
    # Closest non-self: N at 2.0, O at 5.0.
    assert "N" in out.stdout and "2.0000" in out.stdout
    assert "O" in out.stdout and "5.0000" in out.stdout


# ---------------------------------------------------------------------------
# YAML-driven pipeline (slow, byte-equivalence)
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.skipif(
    not (_HAS_FIXTURE and _HAS_TEMPLATE),
    reason="crystal_csp7 fixture or csp7_smoke UNK.ajf missing",
)
def test_cli_pipeline_csp7_r00001_matches_fixture(tmp_path):
    """`abmp-crystal pipeline --config <yaml>` on csp7 R00001 must match
    the Phase B fixture byte-for-byte (legacy engine, -xyz mode)."""
    pytest.importorskip("yaml")
    from .test_regression import _compare_output_dir

    fixture_dir = os.path.join(CRYSTAL_FIXTURE, "R00001")
    cif_src = os.path.join(fixture_dir, "XXXI-MMFF-R00001.cif")

    work = tmp_path / "csp7"
    work.mkdir()
    shutil.copy(cif_src, work / "XXXI-MMFF-R00001.cif")
    shutil.copy(UNK_AJF_TEMPLATE, work / "UNK.ajf")

    cfg_yaml = work / "crystal.yaml"
    cfg_yaml.write_text(f"""\
project_name: csp7_r00001_cli_smoke
output_dir: ./out
inputs:
  - cif: XXXI-MMFF-R00001.cif
    layer: 5
    atoms_in_mol: [32]
cif_engine:
  engine: legacy
fragment:
  cutmode: around
  solutes: [0]
  criteria: 6.0
  molname: [UNK]
  pieda: true
  cmm: true
  getmode: rfile
  template_ajf: ./UNK.ajf
fmo:
  method: MP2
  basis_set: 6-31Gdag
  memory: 6000
  cpfflag: true
  abinit_ver: rev23
  npro: 1
  is_xyz: true
hpc:
  scheduler: PJM
  queue: small
  group: hp190133
  nodes: 12
  proc_per_node: 2
  omp_threads: 24
  elapse: '24:00:00'
  abinit_dir: /data/hp190133/programs/ABINIT-MP
  binary_name: abinitmp_smp
postproc:
  enable: false
""")

    out = _run_cli(
        "pipeline", "--config", str(cfg_yaml),
        cwd=str(work),
    )
    assert out.returncode == 0, f"pipeline failed:\n{out.stdout}\n{out.stderr}"

    out_for_abmp = (
        work / "out" / "XXXI-MMFF-R00001"
        / "cifout" / "layer5" / "pdb" / "for_abmp"
    )
    ref_for_abmp = os.path.join(fixture_dir, "for_abmp")
    _compare_output_dir(str(out_for_abmp), ref_for_abmp)
