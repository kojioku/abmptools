# -*- coding: utf-8 -*-
"""gmx hbond wrapper."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from .._subprocess import run_command

logger = logging.getLogger(__name__)


def run_gmx_hbond(
    *,
    traj: str,
    tpr: str,
    out_xvg: str,
    donor_group: str,
    acceptor_group: str,
    ndx: Optional[str] = None,
    gmx_path: str = "gmx",
    workdir: Optional[str] = None,
) -> str:
    """Run ``gmx hbond`` between *donor_group* and *acceptor_group*.

    The two groups are selected from the input ndx (or default groups
    if ndx is None).
    """
    Path(out_xvg).parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        gmx_path, "hbond",
        "-f", traj, "-s", tpr,
        "-num", out_xvg,
    ]
    if ndx:
        cmd += ["-n", ndx]
    sel = f"{donor_group}\n{acceptor_group}\n"
    run_command(cmd, cwd=workdir, input_text=sel)
    return out_xvg


__all__ = ["run_gmx_hbond"]
