"""
distance_dist.py
----------------
H-bond distance / angle distribution analysis.

Three flavours of plot, all driven from the per-frame ``HBond`` records
already collected by :class:`abmptools.hbond.analyzer.Analyzer`:

A. Overall d(D...A) histogram (single curve, all detected H-bonds).
B. Per-class d(D...A) histogram (overlaid).

   - In ``classify_mode='imc'``: ``COOH-COOH (dual)``, ``COOH-COOH
     (chain/single)``, ``COOH-amide``.
   - In ``classify_mode='generic'``: one curve per
     ``(donor_type, acceptor_type)`` pair.

C. Joint (d, angle) 2D heatmap covering every H-bond regardless of class.

A summary CSV (``mean / median / std / peak / p25 / p75`` per class) and a
long-form histogram CSV (``label, bin_center, count``) are also written so
the same data can be re-plotted downstream.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


@dataclass
class DistanceStats:
    """1D histogram + summary statistics for one class."""
    label: str
    n: int
    mean: float
    median: float
    std: float
    peak: float        # bin center of the most populated bin
    p25: float
    p75: float
    bin_edges: np.ndarray
    counts: np.ndarray


def default_bin_edges(
    d_min: float = 2.0, d_max: float = 3.6, bin_width: float = 0.05,
) -> np.ndarray:
    """Default 2.0–3.6 Å with 0.05 Å bins (32 bins)."""
    n = int(round((d_max - d_min) / bin_width)) + 1
    return np.linspace(d_min, d_max, n)


def compute_distance_stats(
    distances: Sequence[float],
    label: str,
    bin_edges: np.ndarray,
) -> Optional[DistanceStats]:
    """Build a :class:`DistanceStats` from a 1-D sample. ``None`` if empty."""
    arr = np.asarray(list(distances), dtype=float)
    if arr.size == 0:
        return None
    counts, edges = np.histogram(arr, bins=bin_edges)
    if counts.sum() > 0:
        peak_bin = int(np.argmax(counts))
        peak = 0.5 * (edges[peak_bin] + edges[peak_bin + 1])
    else:
        peak = float("nan")
    return DistanceStats(
        label=label,
        n=int(arr.size),
        mean=float(arr.mean()),
        median=float(np.median(arr)),
        std=float(arr.std(ddof=0)),
        peak=float(peak),
        p25=float(np.percentile(arr, 25)),
        p75=float(np.percentile(arr, 75)),
        bin_edges=edges,
        counts=counts,
    )


def _dual_pair_set(cls) -> set:
    """Return ``{frozenset({mol_i, mol_j})}`` of dual COOH-COOH pairs."""
    pairs: set = set()
    if cls is None:
        return pairs
    for r in cls.carboxyl_roles:
        if r.role == "dual":
            for (pm, _pg) in r.dual_partners:
                pairs.add(frozenset({r.mol_index, pm}))
    return pairs


def aggregate_distances_imc(
    frame_results: Iterable,
    bin_edges: np.ndarray,
) -> Dict[str, DistanceStats]:
    """IMC-mode: ``cc`` split into ``dual`` vs ``chain/single`` + ``ca``."""
    buckets: Dict[str, List[float]] = {
        "all": [],
        "COOH-COOH": [],
        "COOH-COOH (dual)": [],
        "COOH-COOH (chain/single)": [],
        "COOH-amide": [],
    }
    for fr in frame_results:
        dual_set = _dual_pair_set(fr.classification)
        for hb in fr.hbonds_cc:
            buckets["all"].append(hb.d_da)
            buckets["COOH-COOH"].append(hb.d_da)
            key = frozenset({hb.donor_mol, hb.acceptor_mol})
            if key in dual_set:
                buckets["COOH-COOH (dual)"].append(hb.d_da)
            else:
                buckets["COOH-COOH (chain/single)"].append(hb.d_da)
        for hb in fr.hbonds_ca:
            buckets["all"].append(hb.d_da)
            buckets["COOH-amide"].append(hb.d_da)
    out: Dict[str, DistanceStats] = {}
    for label, dists in buckets.items():
        s = compute_distance_stats(dists, label, bin_edges)
        if s is not None:
            out[label] = s
    return out


def aggregate_distances_generic(
    frame_results: Iterable,
    bin_edges: np.ndarray,
) -> Dict[str, DistanceStats]:
    """Generic mode: one bucket per ``(donor_type, acceptor_type)``."""
    buckets: Dict[str, List[float]] = {"all": []}
    for fr in frame_results:
        for (dt, at), hbs in fr.hbonds_by_pair_type.items():
            label = f"{dt}->{at}"
            buckets.setdefault(label, [])
            for hb in hbs:
                buckets["all"].append(hb.d_da)
                buckets[label].append(hb.d_da)
    out: Dict[str, DistanceStats] = {}
    for label, dists in buckets.items():
        s = compute_distance_stats(dists, label, bin_edges)
        if s is not None:
            out[label] = s
    return out


def aggregate_distance_angle(
    frame_results: Iterable, mode: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(d_DA[N], angle[N])`` over every detected H-bond."""
    d_list: List[float] = []
    a_list: List[float] = []
    for fr in frame_results:
        if mode == "imc":
            for hb in fr.hbonds_cc:
                d_list.append(hb.d_da); a_list.append(hb.angle)
            for hb in fr.hbonds_ca:
                d_list.append(hb.d_da); a_list.append(hb.angle)
        else:
            for hbs in fr.hbonds_by_pair_type.values():
                for hb in hbs:
                    d_list.append(hb.d_da); a_list.append(hb.angle)
    return np.asarray(d_list, dtype=float), np.asarray(a_list, dtype=float)


# ---- CSV writers ------------------------------------------------------------

def write_distance_stats_csv(
    stats_map: Dict[str, DistanceStats], out_path: str,
) -> None:
    """Per-class summary table."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "label", "n", "mean_DA", "median_DA", "std_DA",
            "peak_DA", "p25_DA", "p75_DA",
        ])
        for s in stats_map.values():
            w.writerow([
                s.label, s.n,
                f"{s.mean:.4f}", f"{s.median:.4f}", f"{s.std:.4f}",
                f"{s.peak:.4f}", f"{s.p25:.4f}", f"{s.p75:.4f}",
            ])


def write_distance_histogram_csv(
    stats_map: Dict[str, DistanceStats], out_path: str,
) -> None:
    """Long-form histogram (``label, bin_center, count``) for re-plotting."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "bin_center_DA", "count"])
        for s in stats_map.values():
            centers = 0.5 * (s.bin_edges[:-1] + s.bin_edges[1:])
            for c, n in zip(centers, s.counts):
                w.writerow([s.label, f"{c:.4f}", int(n)])


# ---- plots ------------------------------------------------------------------

DEFAULT_COLORS_BY_LABEL = {
    "all": "black",
    "COOH-COOH": "red",
    "COOH-COOH (dual)": "red",
    "COOH-COOH (chain/single)": "magenta",
    "COOH-amide": "blue",
}


def plot_distance_histogram(
    stats: DistanceStats, out_path: str, criteria_label: str,
    color: str = "steelblue",
) -> Optional[str]:
    """A: overall d(D...A) histogram with mean / peak annotations."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    centers = 0.5 * (stats.bin_edges[:-1] + stats.bin_edges[1:])
    width = float(stats.bin_edges[1] - stats.bin_edges[0])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(centers, stats.counts, width=width * 0.95, color=color,
           edgecolor="black", linewidth=0.3)
    ax.axvline(stats.mean, color="black", linestyle="--", linewidth=0.8,
               label=f"mean {stats.mean:.2f} Å")
    ax.axvline(stats.peak, color="orange", linestyle="--", linewidth=0.8,
               label=f"peak {stats.peak:.2f} Å")
    ax.set_xlabel("d(D...A) [Å]")
    ax.set_ylabel("count")
    ax.set_title(
        f"H-bond donor-acceptor distance "
        f"({criteria_label}, N={stats.n})"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=100)
    plt.close(fig)
    return out_path


def plot_distance_histogram_classified(
    stats_map: Dict[str, DistanceStats], out_path: str,
    criteria_label: str,
    color_map: Optional[Dict[str, str]] = None,
    exclude: Sequence[str] = ("all", "COOH-COOH"),
) -> Optional[str]:
    """B: overlaid d(D...A) histograms per class (step + fill)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    color_map = color_map or DEFAULT_COLORS_BY_LABEL
    labels = [lbl for lbl in stats_map if lbl not in exclude]
    if not labels:
        return None
    fig, ax = plt.subplots(figsize=(7, 4))
    for lbl in labels:
        s = stats_map[lbl]
        centers = 0.5 * (s.bin_edges[:-1] + s.bin_edges[1:])
        color = color_map.get(lbl)
        ax.step(centers, s.counts, where="mid", linewidth=1.5,
                color=color, label=f"{lbl} (N={s.n})")
        ax.fill_between(centers, s.counts, step="mid", alpha=0.18,
                        color=color)
    ax.set_xlabel("d(D...A) [Å]")
    ax.set_ylabel("count")
    ax.set_title(
        f"H-bond donor-acceptor distance, per class ({criteria_label})"
    )
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=100)
    plt.close(fig)
    return out_path


def plot_distance_angle_2d(
    distances: np.ndarray, angles: np.ndarray, out_path: str,
    criteria_label: str,
    d_bin_width: float = 0.05, angle_bin_width: float = 5.0,
    d_min: float = 2.0, d_max: float = 3.6,
    angle_min: float = 110.0, angle_max: float = 180.0,
) -> Optional[str]:
    """C: 2D heatmap of (d_DA, ∠(D-H...A))."""
    if distances.size == 0:
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    # Expand range to cover any samples beyond the default window.
    d_lo = float(min(d_min, distances.min()))
    d_hi = float(max(d_max, distances.max() + d_bin_width))
    a_lo = float(min(angle_min, angles.min()))
    a_hi = float(max(angle_max, angles.max() + angle_bin_width))
    d_edges = np.arange(d_lo, d_hi + d_bin_width, d_bin_width)
    a_edges = np.arange(a_lo, a_hi + angle_bin_width, angle_bin_width)
    h, xed, yed = np.histogram2d(distances, angles, bins=(d_edges, a_edges))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    pcm = ax.pcolormesh(xed, yed, h.T, cmap="viridis", shading="auto")
    ax.set_xlabel("d(D...A) [Å]")
    ax.set_ylabel("∠(D-H...A) [deg]")
    ax.set_title(
        f"H-bond (distance, angle) joint distribution "
        f"({criteria_label}, N={int(h.sum())})"
    )
    fig.colorbar(pcm, ax=ax, label="count")
    fig.tight_layout()
    fig.savefig(out_path, dpi=100)
    plt.close(fig)
    return out_path
