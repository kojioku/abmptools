# -*- coding: utf-8 -*-
"""matplotlib helpers for formulation analysis outputs."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)


def plot_aggregate_size(
    csv_path: str,
    out_png: str,
) -> str:
    """Stack the max_cluster_size + n_clusters timeseries."""
    import matplotlib.pyplot as plt
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    frames = data[:, 0]
    max_size = data[:, 1]
    n_clusters = data[:, 2]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(frames, max_size, label="max cluster size", color="tab:blue")
    ax.plot(frames, n_clusters, label="n_clusters", color="tab:orange")
    ax.set_xlabel("frame")
    ax.set_ylabel("count")
    ax.set_title("Peptide aggregation timeseries")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    return out_png


def plot_contact_heatmap(
    npy_path: str,
    out_png: str,
    *,
    row_labels: Optional[Sequence[str]] = None,
    col_labels: Optional[Sequence[str]] = None,
) -> str:
    """Plot a residue × species contact-frequency heatmap."""
    import matplotlib.pyplot as plt
    mat = np.load(npy_path)
    fig, ax = plt.subplots(figsize=(max(6, mat.shape[1] * 1.5), max(4, mat.shape[0] * 0.18)))
    im = ax.imshow(mat, aspect="auto", cmap="viridis")
    if row_labels:
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=7)
    if col_labels:
        ax.set_xticks(range(len(col_labels)))
        ax.set_xticklabels(col_labels, rotation=45, ha="right")
    fig.colorbar(im, ax=ax, label="contact frequency")
    ax.set_title("Per-residue contact frequency")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    return out_png


def plot_aggregate_timeseries(
    csv_path: str,
    out_png: str,
    *,
    title_prefix: str = "",
) -> str:
    """論文 Fig 1b/1c 対応の 2-panel plot (max size + % aggregated)。

    CSV format: ``frame,max_size,n_clusters,n_aggregated,pct_aggregated``。
    """
    import matplotlib.pyplot as plt
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    frame = data[:, 0]
    max_size = data[:, 1]
    pct_agg = data[:, 4]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axes[0].plot(frame, max_size, "o-", color="tab:blue", ms=3)
    axes[0].set_ylabel("Max aggregate size")
    axes[0].set_title(f"{title_prefix}Max aggregate size + % aggregated (Hossain 2023 Fig 1b/c)")
    axes[0].grid(alpha=0.3)
    axes[1].plot(frame, pct_agg, "o-", color="tab:red", ms=3)
    axes[1].set_ylabel("% aggregated peptides")
    axes[1].set_xlabel("frame")
    axes[1].set_ylim(0, 105)
    axes[1].axhspan(60, 80, alpha=0.15, color="gray",
                    label="Hossain steady-state range (60-80%)")
    axes[1].legend(loc="lower right")
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    return out_png


def plot_max_size_distribution(
    csv_path: str,
    out_png: str,
    *,
    title_prefix: str = "",
) -> str:
    """論文 Fig 2 node-size 確率に対応する max aggregate size の histogram。"""
    import matplotlib.pyplot as plt
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    max_size = data[:, 1].astype(int)
    sizes = sorted(set(max_size))
    counts = [int((max_size == s).sum()) for s in sizes]
    total = sum(counts)
    fractions = [100 * c / total for c in counts]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(sizes, fractions, color="tab:purple", edgecolor="black")
    for b, c in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1,
                f"{c}f", ha="center", fontsize=8)
    ax.set_xlabel("Max aggregate size")
    ax.set_ylabel("frame %")
    ax.set_title(f"{title_prefix}Max aggregate size distribution (Fig 2 node prob.)")
    ax.set_xticks(sizes)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    return out_png


def plot_per_residue_contacts(
    csv_path: str,
    out_png: str,
    *,
    title_prefix: str = "",
) -> str:
    """論文 Fig 4 対応の per-residue contact stacked bar。

    CSV format: ``res_idx,resname,Peptide_Enhancer,Peptide_BileSalt,Peptide_OtherPeptides``。
    """
    import matplotlib.pyplot as plt
    import csv as _csv
    rows = []
    with open(csv_path) as f:
        r = _csv.reader(f)
        next(r)
        for row in r:
            rows.append((int(row[0]), row[1],
                         float(row[2]), float(row[3]), float(row[4])))
    labels = [f"{rn}{idx}" for idx, rn, *_ in rows]
    enh = [r[2] for r in rows]
    bs = [r[3] for r in rows]
    pep = [r[4] for r in rows]
    x = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(max(7, len(rows) * 0.7), 5))
    w = 0.27
    ax.bar(x - w, enh, w, label="↔ Enhancer", color="tab:blue")
    ax.bar(x,     bs,  w, label="↔ BileSalt", color="tab:red")
    ax.bar(x + w, pep, w, label="↔ Other Peptide", color="tab:gray")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0)
    ax.set_xlabel("Residue (per peptide average)")
    ax.set_ylabel("Mean # contacts (heavy-atom < 0.5 nm)")
    ax.set_title(f"{title_prefix}Per-residue contacts (Hossain 2023 Fig 4)")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    return out_png


def plot_workflow_outputs(
    *,
    out_dir: str,
    results: dict,
) -> dict:
    """``run_analysis`` の各 sub-result から plot 一式を生成。"""
    out = Path(out_dir) / "plots"
    out.mkdir(parents=True, exist_ok=True)
    paths: dict = {}

    if "aggregate" in results and isinstance(results["aggregate"], dict):
        csv = results["aggregate"].get("aggregate_size_csv")
        if csv:
            paths["aggregate_timeseries"] = plot_aggregate_timeseries(
                csv, str(out / "aggregate_timeseries.png"),
            )
            paths["max_size_distribution"] = plot_max_size_distribution(
                csv, str(out / "max_size_distribution.png"),
            )
    if "contacts" in results and isinstance(results["contacts"], dict):
        csv = results["contacts"].get("per_residue_contacts_csv")
        if csv:
            paths["per_residue_contacts"] = plot_per_residue_contacts(
                csv, str(out / "per_residue_contacts.png"),
            )
    return paths


__all__ = [
    "plot_aggregate_size",
    "plot_contact_heatmap",
    "plot_aggregate_timeseries",
    "plot_max_size_distribution",
    "plot_per_residue_contacts",
    "plot_workflow_outputs",
]
