# -*- coding: utf-8 -*-
"""Regression tests for the legacy crystal-FMO pipeline (Phase B fixture).

These tests pin the **current** (v1.22.0) behaviour of the
``abmptools.readcif`` + ``abmptools.pdb2fmo`` toolchain so that the Phase
C refactor (which moves the logic into ``abmptools.crystal.*``) cannot
silently change byte-equivalent outputs.

The pipeline runs ``pdb2fmo -xyz`` (direct-coordinate AJF mode):
the FMO input AJF embeds ``&XYZ`` coordinates directly with full
floating-point precision (e.g. ``2.891228423347176``), avoiding the
``%8.3f`` truncation that an external PDB ReadGeom path would impose.
The XYZ file is read by ``pdb_io.exportardxyzfull`` from a sibling
``<base>.xyz`` next to the PDB; ``readcif`` already emits both files
under ``cifout/layer5/{pdb,xyz}/``, so the test stages copies the XYZ
file alongside the PDB before invoking ``pdb2fmo``.

Fixture: lives in **abmptools-sample** (the csp7 series is private and
must not be bundled in the public abmptools repo). The test points at
``$ABMPTOOLS_SAMPLE_DIR/sample/csp7_ciftest/crystal_reference/phase_b_layer5/``
where ``ABMPTOOLS_SAMPLE_DIR`` defaults to ``~/repos/abmptools-sample``.
The fixture covers 3 structures (small / mid / large/Z'>=2):

    R00001  small  (Natom=832, 1 mol/asym × 832 supercell atoms)
    R00004  mid    (Natom=520)
    R00002  large  (Natom=520 -- _cell_formula_units_Z=8 candidate)

Each structure ships:
    R000NN/XXXI-MMFF-R000NN.cif                            -- input
    R000NN/<base>layer5Zp1-around_ar6.0.{ajf,pdb}          -- expected

The shared driver files (``input_param``, ``segment_data.dat``,
``UNK.ajf``) live one level up in ``phase_b_layer5/`` since they are
identical for all three structures.

If the env var is unset or the directory is missing, the tests
``pytest.skip`` rather than fail.

Marked ``@pytest.mark.slow``: a single layer=5 invocation processes
1331 supercell copies and takes ~10-30 s per structure.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys

import pytest

from .test_regression import _compare_files, _compare_output_dir  # noqa: F401

TESTS_DIR = os.path.dirname(__file__)

_DEFAULT_SAMPLE_DIR = os.path.expanduser("~/repos/abmptools-sample")
_SAMPLE_DIR = os.environ.get("ABMPTOOLS_SAMPLE_DIR", _DEFAULT_SAMPLE_DIR)
CRYSTAL_FIXTURE = os.path.join(
    _SAMPLE_DIR, "sample", "csp7_ciftest", "crystal_reference", "phase_b_layer5"
)

# Each entry: (structure_id, atom_count_per_molecule)
_STRUCTURES = [
    ("R00001", 32),
    ("R00002", 32),
    ("R00004", 32),
]

_HAS_FIXTURE = os.path.isdir(CRYSTAL_FIXTURE) and all(
    os.path.isdir(os.path.join(CRYSTAL_FIXTURE, sid))
    for sid, _ in _STRUCTURES
)

_skip_no_fixture = pytest.mark.skipif(
    not _HAS_FIXTURE,
    reason=(
        f"crystal_csp7 fixture not available at {CRYSTAL_FIXTURE!r}. "
        "Set ABMPTOOLS_SAMPLE_DIR to your abmptools-sample checkout."
    ),
)


def _run(args, cwd, timeout=300):
    """Invoke ``python -m <args>`` in *cwd*."""
    return subprocess.run(
        [sys.executable, "-m"] + args,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


@pytest.mark.slow
@_skip_no_fixture
@pytest.mark.parametrize("structure_id,atom_count", _STRUCTURES)
def test_crystal_csp7_pipeline(tmp_path, structure_id, atom_count):
    """Run readcif -> pdb2fmo (-xyz) and check for_abmp/* matches fixture.

    Pipeline stages:

    1. ``python -m abmptools.readcif -i <cif> -an <N> -l 5``
        emits ``cifout/layer5/pdb/<base>layer5Zp1.pdb`` and
        ``cifout/layer5/xyz/<base>layer5Zp1.xyz`` (full-precision
        coordinates).
    2. Copy the XYZ file next to the PDB so
        ``pdb_io.exportardxyzfull`` (called when ``-xyz`` is set) can
        load coordinates from ``<base>.xyz``.
    3. ``python -m abmptools.pdb2fmo -i <pdb> -p input_param -xyz``
        emits ``for_abmp/<base>layer5Zp1-around_ar6.0.{ajf,pdb}`` --
        the ``.ajf`` carries ``Natom=<N>`` plus a ``&XYZ`` block with
        the full-precision coordinates embedded directly.

    The expected outputs are the byte-equivalent v1.22.0 outputs
    captured in Phase B — the fixture freezes both files (ajf + pdb)
    so the Phase C refactor cannot change them.
    """
    fixture_dir = os.path.join(CRYSTAL_FIXTURE, structure_id)
    cif_name = f"XXXI-MMFF-{structure_id}.cif"

    # Stage 1 inputs: cif at top level, driver files at top level too.
    # The drivers (`input_param` / `segment_data.dat` / `UNK.ajf`) and the
    # cif live in `<sample_dir>/sample/csp7_ciftest/crystal_reference/
    # phase_b_layer5/<...>` (private inputs).
    shutil.copy(os.path.join(fixture_dir, cif_name), tmp_path / cif_name)
    for driver in ("input_param", "segment_data.dat", "UNK.ajf"):
        shutil.copy(os.path.join(CRYSTAL_FIXTURE, driver), tmp_path / driver)

    # Stage 1: CIF -> supercell PDB + XYZ (legacy CLI).
    r1 = _run(
        ["abmptools.readcif", "-i", cif_name,
         "-an", str(atom_count), "-l", "5"],
        cwd=str(tmp_path),
    )
    assert r1.returncode == 0, f"readcif failed:\n{r1.stderr}"

    pdb_dir = tmp_path / "cifout" / "layer5" / "pdb"
    xyz_dir = tmp_path / "cifout" / "layer5" / "xyz"
    pdb_name = f"XXXI-MMFF-{structure_id}layer5Zp1.pdb"
    xyz_name = f"XXXI-MMFF-{structure_id}layer5Zp1.xyz"
    assert (pdb_dir / pdb_name).exists(), (
        f"readcif did not emit expected PDB: {pdb_name}\n"
        f"Found: {sorted(os.listdir(pdb_dir)) if pdb_dir.exists() else '<missing>'}"
    )
    assert (xyz_dir / xyz_name).exists(), (
        f"readcif did not emit expected XYZ: {xyz_name}"
    )

    # Stage 2 prep: drivers + XYZ must sit alongside the PDB for pdb2fmo -xyz.
    for driver in ("input_param", "segment_data.dat", "UNK.ajf"):
        shutil.copy(tmp_path / driver, pdb_dir / driver)
    shutil.copy(xyz_dir / xyz_name, pdb_dir / xyz_name)

    # Stage 2: PDB -> for_abmp/ (-xyz: full-precision direct coordinates).
    r2 = _run(
        ["abmptools.pdb2fmo", "-i", pdb_name, "-p", "input_param", "-xyz"],
        cwd=str(pdb_dir),
    )
    assert r2.returncode == 0, f"pdb2fmo failed:\n{r2.stderr}"

    # Compare the entire for_abmp/ output against the fixture.
    # In abmptools-sample, the expected ajf/pdb live flat in
    # `phase_b_layer5/<system>/`, so file-by-file compare instead of
    # whole-dir compare (the dir contains the cif too).
    out_dir = pdb_dir / "for_abmp"
    base = f"XXXI-MMFF-{structure_id}layer5Zp1-around_ar6.0"
    for ext in (".ajf", ".pdb"):
        out_file = out_dir / f"{base}{ext}"
        ref_file = os.path.join(fixture_dir, f"{base}{ext}")
        assert out_file.exists(), f"Missing output: {out_file}"
        _compare_files(str(out_file), ref_file)


@pytest.mark.slow
@_skip_no_fixture
def test_crystal_csp7_legacy_namespace_import():
    """``abmptools.crystal.legacy`` re-exports must match ``abmptools.*``.

    Phase A skeleton check: ensures the legacy namespace did not silently
    drop a module on Phase C refactor. Cheap import-only test, runs in
    well under a second.
    """
    from abmptools import (
        ajf2config as direct_ajf2config,
        getifiepieda as direct_getifiepieda,
        pdb2fmo as direct_pdb2fmo,
        pdbmodify as direct_pdbmodify,
        readcif as direct_readcif,
    )
    from abmptools.crystal import legacy

    # Identity check: legacy.<name> is the same module as abmptools.<name>.
    assert legacy.readcif is direct_readcif
    assert legacy.pdb2fmo is direct_pdb2fmo
    assert legacy.ajf2config is direct_ajf2config
    assert legacy.pdbmodify is direct_pdbmodify
    assert legacy.getifiepieda is direct_getifiepieda
