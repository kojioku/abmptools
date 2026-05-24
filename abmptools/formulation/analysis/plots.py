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


__all__ = ["plot_aggregate_size", "plot_contact_heatmap"]
