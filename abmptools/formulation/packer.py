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

from ..amorphous.packing import run_packmol

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


__all__ = ["pack_mixed_solution"]
