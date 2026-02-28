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

    for pdb_path, n in zip(pdb_paths, counts):
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

    inp_text = generate_packmol_input(
        pdb_paths, counts, box_size_nm, output_pdb,
        tolerance=tolerance, seed=seed,
    )
    inp_path = os.path.join(build_dir, "packmol.inp")
    Path(inp_path).write_text(inp_text)
    logger.info("Packmol input written to %s", inp_path)

    result = subprocess.run(
        [packmol_bin],
        input=inp_text,
        capture_output=True,
        text=True,
        cwd=build_dir,
    )

    log_path = os.path.join(build_dir, "packmol.log")
    Path(log_path).write_text(result.stdout + "\n" + result.stderr)

    if result.returncode != 0 or not os.path.isfile(output_pdb):
        raise RuntimeError(
            f"Packmol failed (returncode={result.returncode}).\n"
            f"See log: {log_path}\n"
            f"stderr: {result.stderr[:500]}"
        )

    logger.info("Packmol finished: %s", output_pdb)
    return str(Path(output_pdb).resolve())
