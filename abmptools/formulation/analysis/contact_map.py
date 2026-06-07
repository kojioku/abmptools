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


def compute_per_residue_contacts(
    *,
    traj: str,
    top: str,
    out_dir: str,
    n_peptides: int,
    enhancer_resnames: Sequence[str] = (),
    bile_salt_resnames: Sequence[str] = (),
    skip_cap_resnames: Sequence[str] = ("ACE", "NME", "NMe"),
    cutoff_nm: float = 0.5,
    stride: int = 1,
) -> dict:
    """論文 Hossain 2023 Fig 4 準拠の per-residue contact 解析。

    各 peptide で **「同じ位置の残基 (例: 第 i 番目の natural residue)」 を
    全 peptide で平均化**、 cap (ACE/NME) は除外。 ``-> AtomGroup`` の
    heavy-atom min distance を PBC 補正付きで判定。

    Output:
        per_residue_contacts.npy   — shape (n_natural_residues, 3),
                                     列 = [Peptide↔Enhancer, ↔BileSalt, ↔OtherPeptides]
                                     値 = "残基 i の atom 1 つでも cutoff 内に
                                            ある相手原子の数" の **per-peptide 平均** + per-frame 平均
        per_residue_contacts.csv   — 同 + residue label
        per_residue_contacts.json  — 平均 + meta
    """
    from MDAnalysis.lib import distances as mdadist
    import json
    import csv

    MDAnalysis = _lazy_imports()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    u = MDAnalysis.Universe(top, traj)
    protein = u.select_atoms("protein and not name H*")
    if protein.n_atoms == 0:
        raise ValueError("No heavy 'protein' atoms found.")
    n_per_pep = protein.n_atoms // n_peptides

    # peptide × n_peptides を atom index で分割、 各 peptide の residue 列を
    # 抽出。 cap (ACE/NME 等) を除く natural residue だけ残す。
    cap_set = set(skip_cap_resnames)
    per_peptide_res_groups: List[List] = []
    natural_resnames: List[str] = []
    for pep_i in range(n_peptides):
        pep_atoms = protein[pep_i * n_per_pep:(pep_i + 1) * n_per_pep]
        natural_res = [r for r in pep_atoms.residues
                       if r.resname not in cap_set]
        per_peptide_res_groups.append(
            [r.atoms.select_atoms("not name H*") for r in natural_res]
        )
        if pep_i == 0:
            natural_resnames = [r.resname for r in natural_res]
    n_res = len(natural_resnames)
    logger.info(
        "compute_per_residue_contacts: %d peptides × %d natural residues "
        "(cap %s excluded)", n_peptides, n_res, sorted(cap_set),
    )

    # species AtomGroup
    enh_sel = " or ".join(f"resname {r}" for r in enhancer_resnames)
    bs_sel = " or ".join(f"resname {r}" for r in bile_salt_resnames)
    enh_grp = u.select_atoms(f"({enh_sel}) and not name H*") if enh_sel else None
    bs_grp = u.select_atoms(f"({bs_sel}) and not name H*") if bs_sel else None

    cutoff_A = cutoff_nm * 10.0
    # contact_sum[i, k] = ∑frame ∑peptide  (#相手 atom that has any res-atom < cutoff_nm)
    contact_sum = np.zeros((n_res, 3), dtype=np.float64)  # 列: Enh, BS, OtherPep
    n_frames_used = 0
    for ts in u.trajectory[::stride]:
        n_frames_used += 1
        box = ts.dimensions
        for pep_i, res_list in enumerate(per_peptide_res_groups):
            other_pep_atoms = np.concatenate(
                [protein[op * n_per_pep:(op + 1) * n_per_pep].positions
                 for op in range(n_peptides) if op != pep_i],
                axis=0,
            )
            for res_i, res_grp in enumerate(res_list):
                if res_grp.n_atoms == 0:
                    continue
                # Enhancer
                if enh_grp is not None and enh_grp.n_atoms > 0:
                    d = mdadist.distance_array(
                        res_grp.positions, enh_grp.positions, box=box)
                    contact_sum[res_i, 0] += int((d < cutoff_A).any(axis=0).sum())
                # Bile salt
                if bs_grp is not None and bs_grp.n_atoms > 0:
                    d = mdadist.distance_array(
                        res_grp.positions, bs_grp.positions, box=box)
                    contact_sum[res_i, 1] += int((d < cutoff_A).any(axis=0).sum())
                # Other peptides
                d = mdadist.distance_array(
                    res_grp.positions, other_pep_atoms, box=box)
                contact_sum[res_i, 2] += int((d < cutoff_A).any(axis=0).sum())

    # per-peptide 平均 + per-frame 平均
    contact_mean = contact_sum / max(n_frames_used * n_peptides, 1)
    np.save(out / "per_residue_contacts.npy", contact_mean)

    species_cols = ["Peptide_Enhancer", "Peptide_BileSalt", "Peptide_OtherPeptides"]
    with open(out / "per_residue_contacts.csv", "w") as f:
        w = csv.writer(f)
        w.writerow(["res_idx", "resname"] + species_cols)
        for i, rn in enumerate(natural_resnames):
            w.writerow([i + 1, rn] + [f"{contact_mean[i, k]:.3f}" for k in range(3)])

    summary = {
        "n_peptides": int(n_peptides),
        "n_natural_residues": int(n_res),
        "n_frames_analyzed": int(n_frames_used),
        "cutoff_nm": float(cutoff_nm),
        "natural_resnames": natural_resnames,
        "per_residue_contacts": {
            species_cols[k]: [float(contact_mean[i, k]) for i in range(n_res)]
            for k in range(3)
        },
        "totals_per_peptide": {
            species_cols[k]: float(contact_mean[:, k].sum())
            for k in range(3)
        },
    }
    (out / "per_residue_contacts.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False)
    )
    return {
        "per_residue_contacts_npy": str(out / "per_residue_contacts.npy"),
        "per_residue_contacts_csv": str(out / "per_residue_contacts.csv"),
        "per_residue_contacts_json": str(out / "per_residue_contacts.json"),
        "n_residues": n_res,
        "summary": summary,
    }


__all__ = ["compute_contact_map", "compute_per_residue_contacts"]
