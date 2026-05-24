# -*- coding: utf-8 -*-
"""gmx dssp wrapper — secondary structure timeseries.

gmx 2021+ ships ``gmx dssp`` which writes a colour-coded .xpm and a
per-residue/per-frame text format. We parse the xpm into a numpy
array of H/E/T/C fractions over time.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional

import numpy as np

from .._subprocess import run_command

logger = logging.getLogger(__name__)


def run_gmx_dssp(
    *,
    traj: str,
    tpr: str,
    out_xpm: str,
    selection: str = "protein",
    gmx_path: str = "gmx",
    workdir: Optional[str] = None,
) -> str:
    """Run ``gmx dssp`` and return the xpm path."""
    Path(out_xpm).parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        gmx_path, "dssp",
        "-f", traj, "-s", tpr,
        "-o", out_xpm,
        "-sel", selection,
    ]
    run_command(cmd, cwd=workdir, input_text="\n")
    return out_xpm


_XPM_LINE_RE = re.compile(r'"(?P<row>[A-Za-z. ]+)"')


def parse_dssp_xpm(xpm_path: str) -> dict:
    """Parse a gmx dssp xpm into per-frame state counts.

    Returns
    -------
    dict
        ``"H_fraction"``, ``"E_fraction"``, ``"T_fraction"``,
        ``"C_fraction"`` keys → 1D numpy arrays (per frame).
    """
    text = Path(xpm_path).read_text()
    # Crude parse: each data line is a string of single-letter codes per
    # residue, transposed so rows = time, cols = residue.
    rows: list = []
    in_data = False
    for line in text.splitlines():
        if 'static char' in line:
            in_data = True
            continue
        if not in_data:
            continue
        m = _XPM_LINE_RE.search(line)
        if m:
            rows.append(m.group("row"))
    if not rows:
        return {
            "H_fraction": np.array([]),
            "E_fraction": np.array([]),
            "T_fraction": np.array([]),
            "C_fraction": np.array([]),
        }
    # Transpose: rows[i] = "HHHEE..." per frame? actually per residue.
    # Convention: gmx dssp xpm stores residues as rows, frames as cols.
    n_res = len(rows)
    n_frames = max(len(r) for r in rows) if rows else 0
    h = np.zeros(n_frames, dtype=np.uint32)
    e = np.zeros(n_frames, dtype=np.uint32)
    t = np.zeros(n_frames, dtype=np.uint32)
    c = np.zeros(n_frames, dtype=np.uint32)
    for r in rows:
        for j, ch in enumerate(r[:n_frames]):
            if ch == "H":
                h[j] += 1
            elif ch == "E":
                e[j] += 1
            elif ch in ("T", "B", "G", "S"):
                t[j] += 1
            else:
                c[j] += 1
    denom = float(max(n_res, 1))
    return {
        "H_fraction": h.astype(np.float64) / denom,
        "E_fraction": e.astype(np.float64) / denom,
        "T_fraction": t.astype(np.float64) / denom,
        "C_fraction": c.astype(np.float64) / denom,
    }


__all__ = ["parse_dssp_xpm", "run_gmx_dssp"]
