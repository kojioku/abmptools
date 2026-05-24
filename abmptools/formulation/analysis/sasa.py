# -*- coding: utf-8 -*-
"""gmx sasa wrapper."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

from .._subprocess import run_command

logger = logging.getLogger(__name__)


def run_gmx_sasa(
    *,
    traj: str,
    tpr: str,
    out_xvg: str,
    selection: str = "protein",
    surface: str = "protein",
    output_groups: Sequence[str] = ("protein",),
    gmx_path: str = "gmx",
    workdir: Optional[str] = None,
) -> str:
    """Run ``gmx sasa`` and return the xvg path."""
    Path(out_xvg).parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        gmx_path, "sasa",
        "-f", traj, "-s", tpr,
        "-o", out_xvg,
        "-surface", surface,
        "-output", *output_groups,
    ]
    run_command(cmd, cwd=workdir, input_text=f"{selection}\n")
    return out_xvg


__all__ = ["run_gmx_sasa"]
