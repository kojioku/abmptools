# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.pmf
--------------------------
PMF analysis via ``gmx wham``. CG version of :mod:`abmptools.membrane.pmf`.

WHAM is force-field agnostic: it reads ``window.tpr`` (which encodes the
window's pull k and x0) and ``pullx.xvg`` (the pull coord time series).
The temperature is the only physics-dependent parameter.

Both :class:`MembraneCGBuildConfig` (CG) and ``MembraneConfig`` (AA)
expose ``config.gmx_path`` and ``config.equilibration.temperature_K``,
so the AA implementation is duck-compatible. We delegate to it.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

from abmptools.membrane.pmf import run_wham as _aa_run_wham

from .models import MembraneCGBuildConfig

logger = logging.getLogger(__name__)


def run_wham(
    *, windows_dir: Path, analysis_dir: Path,
    config: MembraneCGBuildConfig,
    bootstrap_n: int = 0,
    temperature_K: Optional[float] = None,
) -> Dict[str, Any]:
    """Run ``gmx wham`` on completed CG window output.

    See :func:`abmptools.membrane.pmf.run_wham` for full semantics. CG
    delegates to the AA implementation; only difference is the *config*
    type. Field names match exactly: AA's ``MembraneConfig`` and CG's
    :class:`MembraneCGBuildConfig` both expose ``gmx_path`` and
    ``equilibration.temperature_K``.

    Parameters
    ----------
    windows_dir
        Path to the ``windows/`` directory containing ``win_NNN/``.
    analysis_dir
        Output directory for PMF files (``pmf.xvg`` / ``histo.xvg`` etc).
    config
        :class:`MembraneCGBuildConfig`.
    bootstrap_n
        Number of bootstrap samples for error estimation (default 0 disables).
    temperature_K
        Temperature override; defaults to
        ``config.equilibration.temperature_K``.

    Returns
    -------
    dict
        Keys ``"pmf"``, ``"histo"``, ``"tpr_list"``, ``"pullx_list"``,
        and (when ``bootstrap_n > 0``) ``"bs_result"``.
    """
    return _aa_run_wham(
        windows_dir=str(windows_dir),
        analysis_dir=str(analysis_dir),
        config=config,           # duck-typed: only ``.gmx_path`` and
                                  # ``.equilibration.temperature_K`` accessed
        bootstrap_n=bootstrap_n,
        temperature_K=temperature_K,
    )
