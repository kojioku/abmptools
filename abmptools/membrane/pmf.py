# -*- coding: utf-8 -*-
"""
abmptools.membrane.pmf
----------------------
PMF analysis: WHAM (gmx wham) and / or PyMBAR (optional).

Inputs (from MD output)
-----------------------
- ``windows/win_*/pullx.xvg``  — pull-coord values per window
- ``windows/win_*/pullf.xvg``  — pull forces per window
- ``windows/win_*/window.tpr`` — used by gmx wham for k / x0

Outputs (to ``analysis/``)
--------------------------
- ``pmf.xvg``                  — PMF(z) [kJ/mol]
- ``histo.xvg``                — per-window histograms
- ``bsResult.xvg`` (optional)  — bootstrap error estimate
"""
from __future__ import annotations

import logging
from typing import Any, Dict

from .models import MembraneConfig

logger = logging.getLogger(__name__)


def run_wham(
    *, windows_dir: str, analysis_dir: str, config: MembraneConfig,
) -> Dict[str, Any]:
    """Run ``gmx wham`` on completed window output.

    Returns
    -------
    dict
        Keys ``"pmf"``, ``"histo"``, optionally ``"bs_result"``.
    """
    raise NotImplementedError(
        "Phase A: PMF analysis not yet implemented. "
        "Phase B will invoke 'gmx wham' with the per-window tpr/pullx/pullf "
        "lists, then optionally re-run via PyMBAR for cross-check."
    )


def run_pymbar(
    *, windows_dir: str, analysis_dir: str, config: MembraneConfig,
) -> Dict[str, Any]:
    """Optional MBAR analysis (cross-check vs WHAM)."""
    raise NotImplementedError("Phase B (optional)")
