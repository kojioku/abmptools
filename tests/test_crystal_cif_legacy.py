# -*- coding: utf-8 -*-
"""Unit tests for ``abmptools.crystal.cif_engine_legacy``.

The Python API wrapper ``run_legacy`` must produce the same files as
``python -m abmptools.readcif`` for the Phase B fixture inputs. The
csp7 R00001 CIF lives in **abmptools-sample** (private inputs), so
the test resolves it via ``ABMPTOOLS_SAMPLE_DIR`` (default
``~/repos/abmptools-sample``) and ``pytest.skip`` when absent.
"""
from __future__ import annotations

import os
import shutil

import pytest

_DEFAULT_SAMPLE_DIR = os.path.expanduser("~/repos/abmptools-sample")
_SAMPLE_DIR = os.environ.get("ABMPTOOLS_SAMPLE_DIR", _DEFAULT_SAMPLE_DIR)
CRYSTAL_FIXTURE = os.path.join(
    _SAMPLE_DIR, "sample", "csp7_ciftest", "crystal_reference", "phase_b_layer5"
)

_R00001_CIF = os.path.join(CRYSTAL_FIXTURE, "R00001", "XXXI-MMFF-R00001.cif")
_HAS_FIXTURE = os.path.isfile(_R00001_CIF)


@pytest.mark.slow
@pytest.mark.skipif(
    not _HAS_FIXTURE,
    reason=(
        f"crystal_csp7 fixture not available at {_R00001_CIF!r}. "
        "Set ABMPTOOLS_SAMPLE_DIR to your abmptools-sample checkout."
    ),
)
def test_run_legacy_produces_pdb_and_xyz(tmp_path):
    """``run_legacy`` should emit layer5/pdb + layer5/xyz outputs.

    The PDB output must be byte-equivalent to the same input run via
    ``python -m abmptools.readcif`` (validated by Phase B regression).
    Here we just check filesystem layout + non-empty contents to catch
    refactor regressions without re-running the full ajf pipeline.
    """
    from abmptools.crystal.cif_engine_legacy import run_legacy

    cif_name = "XXXI-MMFF-R00001.cif"
    shutil.copy(_R00001_CIF, tmp_path / cif_name)

    out_dir = run_legacy(
        cif=cif_name,
        layer=5,
        atoms_in_mol=[32],
        odir="cifout",
        cwd=str(tmp_path),
    )

    pdb_dir = out_dir / "pdb"
    xyz_dir = out_dir / "xyz"
    assert pdb_dir.is_dir(), f"missing pdb dir: {pdb_dir}"
    assert xyz_dir.is_dir(), f"missing xyz dir: {xyz_dir}"

    pdb = pdb_dir / "XXXI-MMFF-R00001layer5Zp1.pdb"
    xyz = xyz_dir / "XXXI-MMFF-R00001layer5Zp1.xyz"
    assert pdb.is_file(), f"missing PDB output: {pdb}"
    assert xyz.is_file(), f"missing XYZ output: {xyz}"

    # PDB sanity: csp7 R00001 layer 5 should yield ~16003 atoms (16000 atoms
    # + 1 COMPND header + 1 AUTHOR + 1 END).
    pdb_lines = pdb.read_text().splitlines()
    hetatm_lines = [l for l in pdb_lines if l.startswith("HETATM")]
    assert len(hetatm_lines) == 16000, (
        f"unexpected HETATM count: got {len(hetatm_lines)}, want 16000"
    )

    # XYZ sanity: first line is the atom count (16000).
    first_line = xyz.read_text().splitlines()[0]
    assert first_line.strip() == "16000", (
        f"unexpected XYZ atom count header: {first_line!r}"
    )


def test_run_legacy_invalid_cif_raises(tmp_path):
    """A missing CIF should propagate as an OSError-class exception
    rather than silent success."""
    from abmptools.crystal.cif_engine_legacy import run_legacy

    with pytest.raises((FileNotFoundError, OSError, IndexError)):
        run_legacy(
            cif="does_not_exist.cif",
            layer=5,
            atoms_in_mol=[32],
            cwd=str(tmp_path),
        )
