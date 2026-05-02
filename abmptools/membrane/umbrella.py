# -*- coding: utf-8 -*-
"""
abmptools.membrane.umbrella
---------------------------
Per-window MDP and run-script generation for umbrella sampling.

For each window i ∈ [0, n_windows), write
``windows/win_{i:03d}/window.mdp`` containing a static pull-coord at
``z_min + i * dz``. The actual MD invocation (``gmx grompp`` + ``gmx
mdrun``) is encoded in ``run.sh`` written by :func:`write_run_script`.
"""
from __future__ import annotations

import logging
from typing import Dict

from .models import MembraneConfig

logger = logging.getLogger(__name__)


def write_window_mdps(
    *, config: MembraneConfig, windows_dir: str,
) -> Dict[int, str]:
    """Write per-window MDPs.

    Returns
    -------
    dict
        Mapping window index → absolute path to that window's MDP.
    """
    raise NotImplementedError(
        "Phase A: window MDP writer not yet implemented. "
        "Phase B will iterate config.umbrella.n_windows, compute "
        "z_init for each, and emit window.mdp via render_pull_block."
    )


def write_run_script(
    *,
    config: MembraneConfig,
    top: str, gro: str, ndx: str,
    equil_mdps: Dict[str, str],
    pull_mdp: str,
    window_mdps: Dict[int, str],
    output_dir: str,
) -> str:
    """Write top-level run.sh tying all stages together.

    Returns the absolute path to run.sh.

    The script (POSIX bash) does:
      1. em / nvt / npt equilibration
      2. pulling (pull.tpr → pull.xtc)
      3. extract per-window starting GRO via gmx trjconv
      4. per-window grompp + mdrun (sequential or GNU-parallel-ready)
      5. ``gmx wham`` to produce PMF (or hand off to ``analyze_pmf``)
    """
    raise NotImplementedError(
        "Phase A: run.sh writer not yet implemented. "
        "Phase B will emit a posix bash script with placeholders for "
        "GROMACS executable path / MPI / GPU."
    )
