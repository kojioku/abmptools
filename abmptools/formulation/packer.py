# -*- coding: utf-8 -*-
"""packmol multi-component cubic-box packing for formulation builds.

Thin wrapper over :mod:`abmptools.amorphous.packing` — the same
``generate_packmol_input`` + ``run_packmol`` machinery handles
multi-component mixtures already (tolerance 2.0 Å, packmol exit-173
treated as recoverable). We only add a convenience entry point that
takes ``(component_pdbs, counts)`` and a packmol-inside box edge
(``box_size_nm - inner_box_margin_nm`` so the outer solvatebox shell
fits).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import os
import subprocess

from ..amorphous.packing import _find_packmol, run_packmol

logger = logging.getLogger(__name__)


def pack_mixed_solution(
    *,
    component_pdbs: Sequence[str],
    counts: Sequence[int],
    box_size_nm: float,
    out_pdb: str,
    tolerance_A: float = 2.0,
    inner_box_margin_nm: float = 0.5,
    seed: Optional[int] = None,
    packmol_path: str = "packmol",
) -> str:
    """Pack components into a cubic ``(box_size_nm - inner_box_margin_nm)`` box.

    Returns the absolute path to the resulting mixture PDB.
    """
    if len(component_pdbs) != len(counts):
        raise ValueError(
            f"pack_mixed_solution: component_pdbs ({len(component_pdbs)}) "
            f"and counts ({len(counts)}) length mismatch."
        )
    inner_edge = max(box_size_nm - inner_box_margin_nm, 0.5)
    out = str(Path(out_pdb).resolve())
    build_dir = str(Path(out).parent)
    Path(build_dir).mkdir(parents=True, exist_ok=True)

    abs_pdbs = [str(Path(p).resolve()) for p in component_pdbs]

    logger.info(
        "packmol: %d species into %.2f nm cubic (inner edge), output=%s",
        len(abs_pdbs), inner_edge, out,
    )
    run_packmol(
        pdb_paths=abs_pdbs,
        counts=list(counts),
        box_size_nm=inner_edge,
        output_pdb=out,
        build_dir=build_dir,
        tolerance=tolerance_A,
        seed=seed,
        packmol_path=packmol_path,
    )
    return out


def pack_cluster_templated(
    *,
    peptide_pdbs: Sequence[str],
    peptide_counts: Sequence[int],
    enhancer_pdbs: Sequence[str],
    enhancer_counts: Sequence[int],
    bile_pdbs: Sequence[str],
    bile_counts: Sequence[int],
    box_size_nm: float,
    out_pdb: str,
    core_radius_nm: float = 1.0,
    shell_radius_nm: float = 3.0,
    tolerance_A: float = 2.0,
    inner_box_margin_nm: float = 0.5,
    seed: Optional[int] = None,
    packmol_path: str = "packmol",
) -> str:
    """Pack peptide(s) at the box centre + enhancers in shell + bile outside.

    Produces a pre-formed micelle-like aggregate at t=0:

    - Peptides   ∈ sphere(centre, core_radius_nm)
    - Enhancers  ∈ sphere(centre, shell_radius_nm) − sphere(centre, core_radius_nm)
    - Bile salts ∈ inner box − sphere(centre, shell_radius_nm)

    All radii are in nm. The packmol input uses Å units internally.
    """
    out = str(Path(out_pdb).resolve())
    build_dir = str(Path(out).parent)
    Path(build_dir).mkdir(parents=True, exist_ok=True)

    inner_edge_nm = max(box_size_nm - inner_box_margin_nm, 0.5)
    box_ang = inner_edge_nm * 10.0
    centre = box_ang / 2.0
    core_ang = core_radius_nm * 10.0
    shell_ang = shell_radius_nm * 10.0
    margin = tolerance_A
    low = margin
    high = box_ang - margin

    abs_peptides = [str(Path(p).resolve()) for p in peptide_pdbs]
    abs_enhancers = [str(Path(p).resolve()) for p in enhancer_pdbs]
    abs_bile = [str(Path(p).resolve()) for p in bile_pdbs]

    lines: List[str] = [
        f"tolerance {tolerance_A:.1f}",
        "filetype pdb",
        f"output {out}",
    ]
    if seed is not None:
        lines.append(f"seed {seed}")
    lines.append("")

    # Peptide core: inside sphere(centre, core)
    for pdb, n in zip(abs_peptides, peptide_counts):
        lines.extend([
            f"structure {pdb}",
            f"  number {n}",
            f"  inside sphere {centre:.2f} {centre:.2f} {centre:.2f} {core_ang:.2f}",
            "end structure",
            "",
        ])
    # Enhancer shell: inside sphere(shell) outside sphere(core)
    for pdb, n in zip(abs_enhancers, enhancer_counts):
        lines.extend([
            f"structure {pdb}",
            f"  number {n}",
            f"  inside sphere {centre:.2f} {centre:.2f} {centre:.2f} {shell_ang:.2f}",
            f"  outside sphere {centre:.2f} {centre:.2f} {centre:.2f} {core_ang:.2f}",
            "end structure",
            "",
        ])
    # Bile salts: inner box outside shell sphere
    for pdb, n in zip(abs_bile, bile_counts):
        lines.extend([
            f"structure {pdb}",
            f"  number {n}",
            f"  inside box {low:.2f} {low:.2f} {low:.2f} {high:.2f} {high:.2f} {high:.2f}",
            f"  outside sphere {centre:.2f} {centre:.2f} {centre:.2f} {shell_ang:.2f}",
            "end structure",
            "",
        ])

    inp_text = "\n".join(lines)
    inp_path = os.path.join(build_dir, "packmol.inp")
    Path(inp_path).write_text(inp_text)
    logger.info(
        "cluster-templated packmol: peptide@core(<%.1f nm) + enhancer@shell + bile@outside "
        "(box %.2f nm, %d peptide species, %d enhancer, %d bile)",
        core_radius_nm, inner_edge_nm,
        len(abs_peptides), len(abs_enhancers), len(abs_bile),
    )

    packmol_bin = _find_packmol(packmol_path)
    with open(inp_path, "rb") as inp_fh:
        result = subprocess.run(
            [packmol_bin], stdin=inp_fh, capture_output=True, text=True,
            cwd=build_dir,
        )
    log_path = os.path.join(build_dir, "packmol.log")
    Path(log_path).write_text((result.stdout or "") + "\n" + (result.stderr or ""))
    output_ok = os.path.isfile(out)
    if result.returncode in (0, 173) and output_ok:
        if result.returncode == 173:
            logger.warning("packmol returncode 173 (best-found packing); continuing.")
        return out
    raise RuntimeError(
        f"cluster-templated packmol failed (rc={result.returncode}). "
        f"See {log_path}."
    )


__all__ = ["pack_mixed_solution", "pack_cluster_templated"]
