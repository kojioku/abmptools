# -*- coding: utf-8 -*-
"""Per-residue contact map (peptide ↔ small molecule).

Heavy-atom cutoff (0.5 nm by default, matching Hossain 2023). Outputs
a NumPy array of shape (n_residues, n_smallmol_species) with normalized
contact frequencies.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np

logger = logging.getLogger(__name__)


def _lazy_imports():
    try:
        import MDAnalysis  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "MDAnalysis is required for contact_map. Install via "
            "`pip install abmptools[formulation-analysis]`."
        ) from exc
    return MDAnalysis


def compute_contact_map(
    *,
    traj: str,
    top: str,
    out_dir: str,
    enhancer_resnames: Sequence[str] = (),
    bile_salt_resnames: Sequence[str] = (),
    cutoff_nm: float = 0.5,
    stride: int = 10,
) -> dict:
    """Compute residue → small-molecule contact frequencies.

    Outputs into *out_dir*:
        contact_map.npy        — (n_residues, n_species) float array
        contact_map_labels.txt — row + column labels
    """
    MDAnalysis = _lazy_imports()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    u = MDAnalysis.Universe(top, traj)
    protein = u.select_atoms("protein")
    residues = list(protein.residues)
    if not residues:
        raise RuntimeError("No 'protein' selection in topology.")

    species_groups: List[tuple] = []
    for res in enhancer_resnames:
        ag = u.select_atoms(f"resname {res}")
        if len(ag):
            species_groups.append((f"Enhancer:{res}", ag))
    for res in bile_salt_resnames:
        ag = u.select_atoms(f"resname {res}")
        if len(ag):
            species_groups.append((f"BileSalt:{res}", ag))

    n_res = len(residues)
    n_spec = len(species_groups)
    matrix = np.zeros((n_res, n_spec), dtype=np.uint32)
    cutoff_A = cutoff_nm * 10.0

    frames_used = 0
    for ts in u.trajectory[::stride]:
        frames_used += 1
        for j, (_lab, ag) in enumerate(species_groups):
            if len(ag) == 0:
                continue
            # Pairwise distances residue-by-residue (slow but simple)
            for i, res in enumerate(residues):
                res_atoms = res.atoms
                if len(res_atoms) == 0:
                    continue
                d2 = _min_dist2(res_atoms.positions, ag.positions)
                if d2 <= cutoff_A * cutoff_A:
                    matrix[i, j] += 1

    freq = matrix.astype(np.float64) / max(frames_used, 1)
    np.save(out / "contact_map.npy", freq)
    labels_path = out / "contact_map_labels.txt"
    with open(labels_path, "w") as f:
        f.write("# rows (residues):\n")
        for res in residues:
            f.write(f"{res.resnum}\t{res.resname}\n")
        f.write("# columns (species):\n")
        for lab, _ag in species_groups:
            f.write(f"{lab}\n")
    return {
        "contact_map_npy": str(out / "contact_map.npy"),
        "contact_map_labels": str(labels_path),
        "n_residues": n_res,
        "n_species": n_spec,
        "frames_used": frames_used,
    }


def _min_dist2(a: np.ndarray, b: np.ndarray) -> float:
    """Minimum squared distance between two coordinate sets."""
    if len(a) == 0 or len(b) == 0:
        return float("inf")
    # broadcasting: (n,1,3) - (1,m,3) → (n,m,3)
    diff = a[:, None, :] - b[None, :, :]
    d2 = (diff * diff).sum(axis=-1)
    return float(d2.min())


__all__ = ["compute_contact_map"]
