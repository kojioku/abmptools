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


def per_frame_clusters_heavy_atom(
    peptide_atom_groups: Sequence,
    *,
    box,
    cutoff_nm: float = 0.5,
) -> List[int]:
    """Connected components by PBC-aware **heavy-atom min distance**.

    Hossain 2023 準拠の判定:
    任意の atom pair 1 つでも < cutoff_nm なら同 aggregate。
    NPT box は ``ts.dimensions`` を渡して minimum-image 補正。

    Parameters
    ----------
    peptide_atom_groups
        MDAnalysis ``AtomGroup`` の sequence、 各 peptide の heavy atoms
        を抽出済 (``select_atoms('protein and not name H*')`` 等)。
    box
        ``ts.dimensions`` (frame ごとの NPT box vector、 [a,b,c,α,β,γ] in Å)。
    cutoff_nm
        論文準拠 0.5 nm。

    Returns
    -------
    list of int
        ``peptide_atom_groups`` と同 index の cluster id。
    """
    from MDAnalysis.lib import distances as mdadist
    _, nx = _lazy_imports()
    n = len(peptide_atom_groups)
    g = nx.Graph()
    g.add_nodes_from(range(n))
    cutoff_A = cutoff_nm * 10.0
    for i in range(n):
        for j in range(i + 1, n):
            d = mdadist.distance_array(
                peptide_atom_groups[i].positions,
                peptide_atom_groups[j].positions,
                box=box,
            )
            if d.min() <= cutoff_A:
                g.add_edge(i, j)
    labels = [-1] * n
    for cid, comp in enumerate(nx.connected_components(g)):
        for node in comp:
            labels[node] = cid
    return labels


def per_frame_clusters(
    peptide_coms: np.ndarray,
    *,
    cutoff_nm: float = 0.5,
) -> List[int]:
    """Connected components from a peptide×peptide COM-COM distance matrix.

    Returns one cluster id per peptide (0-indexed).

    Note
    ----
    COM-COM 距離は heavy-atom min distance より粗い判定。 論文準拠の解析には
    :func:`per_frame_clusters_heavy_atom` を使うこと。 本関数は legacy / 単体
    test 用に維持。
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
    n_peptides: int,
    cutoff_nm: float = 0.5,
    stride: int = 1,
    use_heavy_atom: bool = True,
    peptide_selector: str = "protein",
) -> dict:
    """論文 Hossain 2023 Fig 1b/1c 準拠の cluster timeseries 解析。

    判定法 (``use_heavy_atom=True`` 推奨):
    任意の heavy atom pair が < ``cutoff_nm`` なら同 aggregate。
    NPT box は frame ごとに ``ts.dimensions`` で minimum-image 補正。

    Parameters
    ----------
    traj
        trajectory (.xtc / .trr)。 **raw `prod.xtc`** を推奨 (nojump xtc は
        first frame 基準 unwrap で cluster 解析に不適)。
    top
        topology (.gro 推奨。 tpr の version mismatch を避ける)。
    out_dir
        出力 dir。
    n_peptides
        peptide copy 数。 ``protein`` selection を atom index で等分する
        (GROMACS が moltype ごとに resid を 1-N にリセットする慣習に対応)。
    cutoff_nm
        論文準拠 0.5 nm。
    stride
        frame skip。 100 frame サンプリングしたい場合は
        ``stride = traj.n_frames // 100`` を caller 側で計算して渡す。
    use_heavy_atom
        ``True`` (推奨): heavy-atom min distance + PBC で判定 (論文準拠)。
        ``False``: COM-COM 距離 (legacy、 厳しすぎる傾向)。

    Outputs (out_dir に書き出し):
        cluster_states.csv          — frame, peptide_idx, cluster_id, cluster_size
        aggregate_size_timeseries.csv — frame, max_size, n_clusters, n_aggregated, pct_aggregated
        transitions.npy             — (n_peptides, n_peptides) co-cluster 累積
        aggregate_summary.json      — 平均値 + 頻度分布
    """
    from collections import Counter
    import json

    MDAnalysis, _ = _lazy_imports()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    u = MDAnalysis.Universe(top, traj)

    # peptide × n_peptides を atom index で分割 (GROMACS の resid リセット
    # への対応; 各 peptide は同 atom 数前提)
    if use_heavy_atom:
        protein = u.select_atoms(f"({peptide_selector}) and not name H*")
    else:
        protein = u.select_atoms(peptide_selector)
    if protein.n_atoms == 0:
        raise ValueError(
            f"No atoms found for peptide_selector={peptide_selector!r}. "
            "For Amber route use 'protein', for whole-peptide GAFF use "
            "e.g. 'resname OCT'."
        )
    if protein.n_atoms % n_peptides != 0:
        logger.warning(
            "protein.n_atoms (%d) is not divisible by n_peptides (%d); "
            "atom-index split may misalign. Verify the topology.",
            protein.n_atoms, n_peptides,
        )
    n_per_pep = protein.n_atoms // n_peptides
    peptide_groups = [
        protein[i * n_per_pep:(i + 1) * n_per_pep] for i in range(n_peptides)
    ]
    logger.info(
        "compute_aggregate_transitions: %d peptides × %d heavy atoms (cutoff=%.2f nm)",
        n_peptides, n_per_pep, cutoff_nm,
    )

    cluster_rows: List[List] = []
    size_rows: List[List] = []
    co_cluster = np.zeros((n_peptides, n_peptides), dtype=np.uint32)
    max_size_history: List[int] = []

    for ts in u.trajectory[::stride]:
        if use_heavy_atom:
            labels = per_frame_clusters_heavy_atom(
                peptide_groups, box=ts.dimensions, cutoff_nm=cutoff_nm,
            )
        else:
            coms = np.array([g.center_of_mass() / 10.0 for g in peptide_groups])
            labels = per_frame_clusters(coms, cutoff_nm=cutoff_nm)
        sizes = cluster_size_per_peptide(labels)
        for i in range(n_peptides):
            for j in range(n_peptides):
                if labels[i] == labels[j]:
                    co_cluster[i, j] += 1
            cluster_rows.append([ts.frame, i, labels[i], sizes[i]])
        max_size = max(sizes) if sizes else 0
        n_clusters = len(set(labels))
        n_aggregated = sum(1 for s in sizes if s >= 2)
        pct_aggregated = 100.0 * n_aggregated / n_peptides
        size_rows.append([ts.frame, max_size, n_clusters,
                          n_aggregated, pct_aggregated])
        max_size_history.append(max_size)

    np.save(out / "transitions.npy", co_cluster)
    _write_csv(out / "cluster_states.csv",
               ["frame", "peptide_idx", "cluster_id", "cluster_size"],
               cluster_rows)
    _write_csv(out / "aggregate_size_timeseries.csv",
               ["frame", "max_size", "n_clusters",
                "n_aggregated", "pct_aggregated"],
               size_rows)

    # Summary statistics + 頻度分布 (論文 Fig 1c 対応)
    pct_agg_arr = np.array([r[4] for r in size_rows])
    n_agg_arr = np.array([r[3] for r in size_rows])
    size_freq = Counter(max_size_history)
    total = sum(size_freq.values())
    size_distribution = {
        str(s): {"count": int(c), "fraction": float(c / total)}
        for s, c in sorted(size_freq.items())
    }
    summary = {
        "n_peptides": int(n_peptides),
        "n_frames_analyzed": int(len(size_rows)),
        "cutoff_nm": float(cutoff_nm),
        "use_heavy_atom": bool(use_heavy_atom),
        "max_size": {
            "mean": float(np.mean(max_size_history)),
            "max": int(max(max_size_history)),
            "mode": int(max(size_freq, key=size_freq.get)),
        },
        "pct_aggregated": {
            "mean": float(pct_agg_arr.mean()),
            "min": float(pct_agg_arr.min()),
            "max": float(pct_agg_arr.max()),
        },
        "n_aggregated_per_frame": {
            "mean": float(n_agg_arr.mean()),
        },
        "max_size_distribution": size_distribution,
    }
    (out / "aggregate_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False)
    )

    return {
        "transitions_npy": str(out / "transitions.npy"),
        "cluster_states_csv": str(out / "cluster_states.csv"),
        "aggregate_size_csv": str(out / "aggregate_size_timeseries.csv"),
        "summary_json": str(out / "aggregate_summary.json"),
        "n_peptides": n_peptides,
        "summary": summary,
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
