"""Plot NMR (Yuan et al. 2015, Mol. Pharm. 12, 4518 Figure 5 + Table 1) vs
MD per-COOH H-bond classification (4 species).

Layout (3 rows):
  Top:    Yuan Figure 5 image (CPMAS 13C NMR of amorphous IMC, carboxyl region)
  Middle: Yuan Table 1 deconvolution peak areas (dimer / chain / amide / free)
  Bottom: MD per-COOH percentages at the same tentative chemical-shift positions

All three panels share the same chemical-shift x-axis so the deconvolution
peaks line up vertically.

Reference:
  Yuan X. et al., "Hydrogen Bonding Interactions in Amorphous Indomethacin
  ... studied using 13C Solid-State NMR", Mol. Pharm. 2015, 12, 4518-4528.
  doi: 10.1021/acs.molpharmaceut.5b00705

Usage:
  cd sample/hbond/imc_amorphous && python plot_nmr_comparison.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread

HERE = Path(__file__).resolve().parent
SUMMARY_CSV = HERE / "output" / "imc_hbond_summary.csv"
NMR_IMG = Path("/home/okuwaki/llm-project/SI/imc-bond-nmr.png")
OUT_PNG = HERE / "output" / "imc_hbond_nmr_comparison.png"

# Yuan 2015 Table 1 — amorphous IMC deconvolution
YUAN_2015 = {
    "dual":   {"ppm": 179.3, "pct": 58.5, "label": "cyclic dimer"},
    "chain":  {"ppm": 176.3, "pct": 15.2, "label": "carboxylic acid chain"},
    "single": {"ppm": 172.4, "pct": 18.9, "label": "COOH-amide"},
    "free":   {"ppm": 170.4, "pct":  7.5, "label": "free COOH"},
}

color_map = {
    "dual":   "#e74c3c",   # red
    "chain":  "#c0398a",   # magenta (matches Mol_Name DEFAULT_COLORS "Magenta")
    "single": "#3498db",   # blue
    "free":   "#7f7f7f",   # gray
}

# ---- Read MD summary (last record) ---------------------------------------
with open(SUMMARY_CSV) as f:
    r = list(csv.DictReader(f))[-1]
n_total = int(r["n_carboxyls"])
md_counts = {
    "dual":   int(r["n_carb_dual"]),
    "chain":  int(r["n_carb_chain"]),
    "single": int(r["n_carb_single"]),
    "free":   int(r["n_carb_free"]),
}
md_ratios = {
    "dual":   float(r["ratio_carb_dual"]),
    "chain":  float(r["ratio_carb_chain"]),
    "single": float(r["ratio_carb_single"]),
    "free":   float(r["ratio_carb_free"]),
}

# ---- Plot -----------------------------------------------------------------
fig, (ax_img, ax_nmr, ax_md) = plt.subplots(
    3, 1, figsize=(10, 10),
    gridspec_kw={"height_ratios": [1.4, 1.0, 1.0]},
)

# Row 1: NMR Figure 5 image
ax_img.imshow(imread(NMR_IMG))
ax_img.axis("off")
ax_img.set_title(
    "Yuan et al. 2015, Figure 5 — CPMAS 13C NMR of amorphous IMC carboxyl region",
    fontsize=11, pad=8,
)

# Row 2: Yuan Table 1 deconvolution bars
ax_nmr.set_xlim(184, 168)
ax_nmr.set_ylim(0, 65)
for role, info in YUAN_2015.items():
    ppm, pct = info["ppm"], info["pct"]
    ax_nmr.bar(
        ppm, pct, width=0.9, color=color_map[role], alpha=0.85,
        edgecolor="black", linewidth=0.5,
        label="{lbl} @ {ppm} ppm: {pct:.1f}%".format(
            lbl=info["label"], ppm=ppm, pct=pct,
        ),
    )
    ax_nmr.text(
        ppm, pct + 1.2, "{:.1f}%".format(pct),
        ha="center", va="bottom", fontsize=10, fontweight="bold",
        color=color_map[role],
    )
ax_nmr.set_ylabel("% peak area (Yuan Table 1)", fontsize=10)
ax_nmr.set_title(
    "Experiment: Yuan 2015 Table 1 deconvolution (Gaussian + LSQ fit)",
    fontsize=11, pad=8,
)
ax_nmr.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
ax_nmr.grid(True, alpha=0.25, axis="y")
ax_nmr.set_xticklabels([])

# Row 3: MD per-COOH bars at the same chemical-shift positions
ax_md.set_xlim(184, 168)
ax_md.set_ylim(0, 65)
for role, info in YUAN_2015.items():
    ppm = info["ppm"]
    pct = md_ratios[role] * 100.0
    ax_md.bar(
        ppm, pct, width=0.9, color=color_map[role], alpha=0.85,
        edgecolor="black", linewidth=0.5,
        label="{lbl}: {n}/{N} ({pct:.1f}%)".format(
            lbl=info["label"], n=md_counts[role], N=n_total, pct=pct,
        ),
    )
    ax_md.text(
        ppm, pct + 1.2, "{:.1f}%".format(pct),
        ha="center", va="bottom", fontsize=10, fontweight="bold",
        color=color_map[role],
    )
ax_md.set_xlabel("Chemical shift / ppm (Yuan 2015 Table 1 positions)", fontsize=10)
ax_md.set_ylabel("% of COOH (MD)", fontsize=10)
ax_md.set_title(
    "MD: per-COOH H-bond state at T=450 K (N={} COOH, Luzar-Chandler)".format(
        n_total,
    ),
    fontsize=11, pad=8,
)
ax_md.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
ax_md.grid(True, alpha=0.25, axis="y")

# Tight layout
fig.tight_layout()
fig.subplots_adjust(hspace=0.32)
fig.savefig(OUT_PNG, dpi=120, bbox_inches="tight")
print("saved:", OUT_PNG)
print()
print("MD vs NMR (Yuan 2015):")
print("  {:8s}  {:>6s}  {:>6s}".format("species", "NMR%", "MD%"))
for role in ("dual", "chain", "single", "free"):
    print("  {:8s}  {:6.1f}  {:6.1f}".format(
        role, YUAN_2015[role]["pct"], md_ratios[role] * 100,
    ))
