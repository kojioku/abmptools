# -*- coding: utf-8 -*-
"""
abmptools.amorphous.packing
-----------------------------
Generate Packmol input and run Packmol to create a packed mixture PDB.
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


def _find_packmol(packmol_path: str = "packmol") -> str:
    """Locate the Packmol binary."""
    resolved = shutil.which(packmol_path)
    if resolved is None:
        raise FileNotFoundError(
            f"Packmol not found at '{packmol_path}'. "
            "Install Packmol or specify --packmol_path."
        )
    return resolved


def generate_packmol_input(
    pdb_paths: List[str],
    counts: List[int],
    box_size_nm: float,
    output_pdb: str,
    tolerance: float = 2.0,
    seed: Optional[int] = None,
    cluster_pdb: Optional[str] = None,
) -> str:
    """Generate a Packmol input string.

    Parameters
    ----------
    pdb_paths : list of str
        Paths to single-molecule PDB files for each component.
    counts : list of int
        Number of molecules for each component.
    box_size_nm : float
        Cubic box edge length [nm].
    output_pdb : str
        Path for the output mixture PDB.
    tolerance : float
        Minimum distance between molecules [Angstrom].
    seed : int or None
        Random seed for Packmol.
    cluster_pdb : str, optional
        Path to a pre-built rigid cluster PDB (e.g. water trimer with
        H-bond triangle). When set, the cluster is placed FIRST at box
        center with packmol's ``fixed`` constraint (no rotation), and
        ``counts`` should already exclude the cluster's atoms from the
        corresponding component (e.g. subtract 3 from the water count for
        a water-trimer cluster). The cluster contributes the same mol-type
        as one of ``pdb_paths`` (a trimer = 3 water mols), so total system
        atom count = cluster_atoms + sum(per-mol-atoms * counts).

    Returns
    -------
    str
        Packmol input file content.
    """
    box_ang = box_size_nm * 10.0  # nm → Angstrom
    margin = tolerance  # keep molecules away from box edges
    low = margin
    high = box_ang - margin

    lines = [
        f"tolerance {tolerance:.1f}",
        "filetype pdb",
        f"output {output_pdb}",
    ]
    if seed is not None:
        lines.append(f"seed {seed}")
    lines.append("")

    # Optional pre-built cluster: pass the cluster as a `fixed` rigid body
    # at box center, then place remaining mols around it. Used for UDF-route
    # `cluster_file` equivalence (e.g. water trimer with H-bond geometry that
    # packmol's per-molecule placement can't construct).
    if cluster_pdb is not None:
        center = box_ang / 2.0
        lines.extend([
            f"structure {cluster_pdb}",
            f"  number 1",
            f"  fixed {center:.2f} {center:.2f} {center:.2f} 0. 0. 0.",
            "end structure",
            "",
        ])

    for pdb_path, n in zip(pdb_paths, counts):
        if n <= 0:
            continue
        lines.extend([
            f"structure {pdb_path}",
            f"  number {n}",
            f"  inside box {low:.2f} {low:.2f} {low:.2f} {high:.2f} {high:.2f} {high:.2f}",
            "end structure",
            "",
        ])

    return "\n".join(lines)


def run_packmol(
    pdb_paths: List[str],
    counts: List[int],
    box_size_nm: float,
    output_pdb: str,
    build_dir: str,
    tolerance: float = 2.0,
    seed: Optional[int] = None,
    packmol_path: str = "packmol",
    cluster_pdb: Optional[str] = None,
) -> str:
    """Run Packmol to create a packed mixture PDB.

    Parameters
    ----------
    pdb_paths : list of str
        Paths to single-molecule PDB files.
    counts : list of int
        Number of molecules for each component.
    box_size_nm : float
        Cubic box edge length [nm].
    output_pdb : str
        Path for the output mixture PDB.
    build_dir : str
        Directory for Packmol input file.
    tolerance : float
        Minimum distance between molecules [Angstrom].
    seed : int or None
        Random seed.
    packmol_path : str
        Path to the Packmol binary.

    Returns
    -------
    str
        Absolute path to the generated mixture PDB.

    Raises
    ------
    RuntimeError
        If Packmol fails.
    """
    packmol_bin = _find_packmol(packmol_path)
    os.makedirs(build_dir, exist_ok=True)

    abs_pdb_paths = [str(Path(p).resolve()) for p in pdb_paths]
    abs_output_pdb = str(Path(output_pdb).resolve())

    abs_cluster_pdb = str(Path(cluster_pdb).resolve()) if cluster_pdb else None
    inp_text = generate_packmol_input(
        abs_pdb_paths, counts, box_size_nm, abs_output_pdb,
        tolerance=tolerance, seed=seed,
        cluster_pdb=abs_cluster_pdb,
    )
    inp_path = os.path.join(build_dir, "packmol.inp")
    Path(inp_path).write_text(inp_text)
    logger.info("Packmol input written to %s", inp_path)

    with open(inp_path, "rb") as inp_fh:
        result = subprocess.run(
            [packmol_bin],
            stdin=inp_fh,
            capture_output=True,
            text=True,
            cwd=build_dir,
        )

    log_path = os.path.join(build_dir, "packmol.log")
    Path(log_path).write_text(result.stdout + "\n" + result.stderr)

    # Packmol exit-code 173 means "ENDED WITHOUT PERFECT PACKING": the
    # solver couldn't satisfy the tolerance everywhere but still wrote
    # the best-found mixture.pdb to disk. That's a usable starting
    # configuration — the downstream EM stage in the 5-stage protocol
    # cleans up residual close contacts within the first few hundred
    # steps. Treat it as a warning and continue, only failing on
    # ``returncode > 0 and missing output_pdb`` or on hard exit codes
    # (segfaults, file errors, etc.). Exit code 0 is the "perfect
    # packing achieved" path; 173 is the "best found" path; anything
    # else is genuinely fatal.
    output_ok = os.path.isfile(output_pdb)
    if result.returncode == 0 and output_ok:
        pass
    elif result.returncode == 173 and output_ok:
        logger.warning(
            "Packmol returned 173 (ENDED WITHOUT PERFECT PACKING) but "
            "wrote %s. Continuing with the best-found configuration; "
            "the EM stage typically resolves any residual overlap. "
            "If you need a tighter pack, raise tolerance or nloop. "
            "Log: %s", output_pdb, log_path)
    else:
        raise RuntimeError(
            f"Packmol failed (returncode={result.returncode}).\n"
            f"See log: {log_path}\n"
            f"stderr: {result.stderr[:500]}"
        )

    logger.info("Packmol finished: %s", output_pdb)
    return str(Path(output_pdb).resolve())
