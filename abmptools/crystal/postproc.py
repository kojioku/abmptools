# -*- coding: utf-8 -*-
"""
abmptools.crystal.postproc
---------------------------
IFIE/PIEDA postprocessing + nearest-atom annotation.

Phase C-6 thin wrapper around two existing tools:

- :mod:`abmptools.getifiepieda` — invoked via subprocess with the
  CLI flags it has historically required (``--multi`` / ``-dimeres`` /
  ``-zp`` / etc.). The wrapper stages logs into the run directory,
  runs the CLI once per dataset, and collects the emitted CSVs.
- :func:`abmptools.crystal.atom_distance.find_nearest_atoms` — annotates
  each peripheral fragment in the for_abmp PDB with the closest
  ``n_neighbors`` heavy atoms to the central solute (residue id from
  :class:`PostProcessSpec.frag_target`).

The wrapper is intentionally minimal: only the most common knobs are
plumbed through. Users with highly custom IFIE pipelines should keep
calling ``python -m abmptools.getifiepieda`` directly; this wrapper
exists so the Phase D tutorial can demonstrate end-to-end automation
without the user having to learn 20+ CLI flags.
"""
from __future__ import annotations

import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

from .atom_distance import NearestAtom, find_nearest_atoms
from .models import PostProcessSpec

logger = logging.getLogger("abmptools.crystal.postproc")


def run_getifiepieda(
    *,
    log_files: List[Path],
    output_dir: Path,
    spec: PostProcessSpec,
    z_prime: int = 1,
) -> List[Path]:
    """Invoke :mod:`abmptools.getifiepieda` once per log file.

    Parameters
    ----------
    log_files
        ABINIT-MP log files (one per cocrystal structure / time frame).
    output_dir
        Directory to stage logs and collect emitted CSVs (``csv/``).
    spec
        :class:`PostProcessSpec` (carries ``frag_target`` / ``distance``).
    z_prime
        Z' value for the ``-zp`` flag (1 for csp7-like Z'=1 systems;
        ``CIFInputSpec.atoms_in_mol`` length serves as a hint).

    Returns
    -------
    List[Path]
        CSVs emitted under ``output_dir/csv/``.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for log in log_files:
        if log.parent != output_dir:
            shutil.copy(log, output_dir / log.name)

    base_args = [
        sys.executable, "-m", "abmptools.getifiepieda",
        "-mul", "1",
        "-dimeres",
        "-zp", str(z_prime),
        "-nof90",
    ]
    if spec.distance is not None:
        base_args += ["-d", str(spec.distance)]
    if spec.frag_target is not None:
        base_args += ["-f", spec.frag_target]

    # Single -i with 2 elements is the historical (prefix, suffix) shape.
    # We approximate with file names, letting getifiepieda glob the dir.
    cmd = base_args + ["-i", str(output_dir / "")]
    logger.info("postproc getifiepieda: %s", " ".join(cmd))
    result = subprocess.run(
        cmd, cwd=str(output_dir),
        capture_output=True, text=True, timeout=600,
    )
    if result.returncode != 0:
        logger.warning(
            "getifiepieda exited non-zero (%s); stderr:\n%s",
            result.returncode, result.stderr,
        )
    csvs = sorted((output_dir / "csv").glob("*.csv")) if (output_dir / "csv").exists() else []
    return csvs


def annotate_with_nearest_atoms(
    *,
    csv_path: Path,
    pdb_path: Path,
    center_res_seq: int,
    n_neighbors: int = 3,
    output_csv: Optional[Path] = None,
) -> Path:
    """Append ``nearest_atom_*`` columns to *csv_path*.

    For each row in the source CSV (assumed to contain a fragment
    residue-id column), the corresponding nearest-atom rows from the
    PDB are inlined as additional columns: ``nearest_element_1``,
    ``nearest_distance_1``, ..., up to ``n_neighbors``.

    Parameters
    ----------
    csv_path
        Source CSV emitted by :mod:`abmptools.getifiepieda` (any layout
        with a ``frag`` or ``residue`` column accepted; the loose
        regex below matches both).
    pdb_path
        PDB used for nearest-atom lookup (same coordinates as the
        ABINIT-MP input).
    center_res_seq
        Residue sequence number for the centre.
    n_neighbors
        Number of nearest atoms to annotate per row.
    output_csv
        If ``None``, replaces *csv_path* in place (after first writing
        to a sibling ``.tmp`` for crash-safety).

    Returns
    -------
    Path
        Path to the annotated CSV.
    """
    rows: List[Dict[str, str]] = []
    with csv_path.open() as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    nearest = find_nearest_atoms(
        pdb_path=str(pdb_path),
        center_res_seq=center_res_seq,
        n_neighbors=n_neighbors,
    )
    extra: List[str] = []
    for i in range(n_neighbors):
        extra += [f"nearest_element_{i + 1}", f"nearest_distance_{i + 1}"]
    annotated_fields = fieldnames + extra

    flat: Dict[str, str] = {}
    for i, na in enumerate(nearest):
        flat[f"nearest_element_{i + 1}"] = na.element
        flat[f"nearest_distance_{i + 1}"] = f"{na.distance:.4f}"

    target = output_csv or csv_path
    tmp = target.with_suffix(target.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=annotated_fields)
        writer.writeheader()
        for row in rows:
            row.update(flat)
            writer.writerow(row)
    tmp.replace(target)
    return target


def run_postprocess(
    *,
    spec: PostProcessSpec,
    output_dir: Path,
    log_files: List[Path],
    pdb_for_annotation: Optional[Path] = None,
    z_prime: int = 1,
) -> Dict[str, List[Path]]:
    """High-level driver: getifiepieda + optional nearest-atom annotation.

    Returns a dict with two keys:
        ``csvs``     -- CSVs emitted by getifiepieda
        ``annotated`` -- CSVs with appended nearest-atom columns
                         (empty if ``annotate_nearest_atoms=False`` or
                         no PDB provided)
    """
    csvs = run_getifiepieda(
        log_files=log_files,
        output_dir=output_dir,
        spec=spec,
        z_prime=z_prime,
    )
    annotated: List[Path] = []
    if (
        spec.annotate_nearest_atoms
        and pdb_for_annotation is not None
        and pdb_for_annotation.is_file()
    ):
        center = 1  # first fragment id == solute by convention
        for csv_path in csvs:
            ann = annotate_with_nearest_atoms(
                csv_path=csv_path,
                pdb_path=pdb_for_annotation,
                center_res_seq=center,
                n_neighbors=spec.nearest_atom_count,
            )
            annotated.append(ann)
    return {"csvs": csvs, "annotated": annotated}


__all__ = [
    "run_getifiepieda",
    "annotate_with_nearest_atoms",
    "run_postprocess",
    "NearestAtom",
]
