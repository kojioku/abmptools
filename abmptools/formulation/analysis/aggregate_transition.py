# -*- coding: utf-8 -*-
"""Aggregate transition matrix + cluster timeseries.

Hossain 2023 methodology:
- Per-frame heavy-atom contact graph between peptide copies
  (cutoff 0.5 nm)
- Connected components → cluster labels
- Track each peptide's cluster size + cluster id over time
- Emit transition matrix (peptide_idx_t → peptide_idx_t+1 same cluster?)

This module imports MDAnalysis + networkx lazily; both live in
``abmptools[formulation-analysis]``.
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
        import networkx  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "MDAnalysis + networkx are required for aggregate transition "
            "analysis. Install via "
            "`pip install abmptools[formulation-analysis]` or "
            "`pip install 'MDAnalysis>=2.6,<3.0' 'networkx>=3.0'`."
        ) from exc
    return MDAnalysis, networkx


def per_frame_clusters(
    peptide_coms: np.ndarray,
    *,
    cutoff_nm: float = 0.5,
) -> List[int]:
    """Connected components from a peptide×peptide COM-COM distance matrix.

    Returns one cluster id per peptide (0-indexed).
    """
    _, nx = _lazy_imports()
    n = peptide_coms.shape[0]
    g = nx.Graph()
    g.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(peptide_coms[i] - peptide_coms[j]))
            if d <= cutoff_nm:
                g.add_edge(i, j)
    labels = [-1] * n
    for cid, comp in enumerate(nx.connected_components(g)):
        for node in comp:
            labels[node] = cid
    return labels


def cluster_size_per_peptide(labels: Sequence[int]) -> List[int]:
    """Return per-peptide cluster size (peptide_i ∈ cluster of size k → k)."""
    counts: dict = {}
    for c in labels:
        counts[c] = counts.get(c, 0) + 1
    return [counts[c] for c in labels]


def compute_aggregate_transitions(
    *,
    traj: str,
    top: str,
    out_dir: str,
    cutoff_nm: float = 0.5,
    peptide_resname_first: str = "ALA",  # placeholder; passed externally
    stride: int = 1,
) -> dict:
    """Compute per-frame cluster labels + transition matrix.

    Outputs (written into *out_dir*):
        cluster_states.csv         — frame, peptide_idx, cluster_id, size
        aggregate_size_timeseries.csv — frame, max_cluster_size, n_clusters
        transitions.npy            — (n_peptides, n_peptides) co-cluster count

    Returns paths as a dict.
    """
    MDAnalysis, _ = _lazy_imports()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    u = MDAnalysis.Universe(top, traj)
    # Use protein residues as "peptide" — for AA work, every chain in
    # the protein selection is treated as one peptide copy.
    sel = u.select_atoms("protein")
    chains = sorted(set(sel.residues.segments.segids))
    peptide_groups = [
        u.select_atoms(f"protein and segid {sid}") for sid in chains
    ]
    n_pep = len(peptide_groups)
    if n_pep < 2:
        logger.warning(
            "compute_aggregate_transitions: only %d peptide copies — "
            "transition matrix is trivial.", n_pep,
        )

    cluster_rows: List[List] = []
    size_rows: List[List] = []
    co_cluster = np.zeros((n_pep, n_pep), dtype=np.uint32)

    for ts in u.trajectory[::stride]:
        coms = np.array([g.center_of_mass() / 10.0 for g in peptide_groups])  # Å → nm
        labels = per_frame_clusters(coms, cutoff_nm=cutoff_nm)
        sizes = cluster_size_per_peptide(labels)
        # accumulate co-cluster matrix
        for i in range(n_pep):
            for j in range(n_pep):
                if labels[i] == labels[j]:
                    co_cluster[i, j] += 1
        for i in range(n_pep):
            cluster_rows.append([ts.frame, i, labels[i], sizes[i]])
        size_rows.append([ts.frame, max(sizes) if sizes else 0,
                          len(set(labels))])

    np.save(out / "transitions.npy", co_cluster)
    _write_csv(out / "cluster_states.csv",
               ["frame", "peptide_idx", "cluster_id", "cluster_size"],
               cluster_rows)
    _write_csv(out / "aggregate_size_timeseries.csv",
               ["frame", "max_cluster_size", "n_clusters"],
               size_rows)
    return {
        "transitions_npy": str(out / "transitions.npy"),
        "cluster_states_csv": str(out / "cluster_states.csv"),
        "aggregate_size_csv": str(out / "aggregate_size_timeseries.csv"),
        "n_peptides": n_pep,
    }


def _write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence]) -> None:
    with open(path, "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(str(c) for c in row) + "\n")


__all__ = [
    "cluster_size_per_peptide",
    "compute_aggregate_transitions",
    "per_frame_clusters",
]
