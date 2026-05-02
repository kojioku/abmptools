# -*- coding: utf-8 -*-
"""
abmptools.membrane.pulling
--------------------------
Reaction-coordinate generation: peptide pulled along z through the bilayer.

The output of stage 4 is a single MD trajectory (``pull.xtc``) plus
snapshots at uniform z-spacing that seed each US window.

Why pulling
-----------
Each US window needs an initial structure where the peptide centre-of-mass
sits near the window's reference z. The cheapest way to get those is a
single steered-MD run with constant-velocity pulling:

    pull-coord1-init = <initial z>
    pull-coord1-rate = <nm/ps>
    pull-coord1-k    = <large k>

Then we extract frames at z = z_min, z_min + dz, … z_max from the
trajectory.
"""
from __future__ import annotations

import logging

from .models import MembraneConfig

logger = logging.getLogger(__name__)


def write_pulling_mdp(*, config: MembraneConfig, pull_dir: str) -> str:
    """Write pull.mdp and return its absolute path."""
    raise NotImplementedError(
        "Phase A: pulling MDP writer not yet implemented. "
        "Phase B will compose pull.mdp from npt.mdp + render_pull_block "
        "with constant pull-rate."
    )


def extract_window_frames(
    *, pull_xtc: str, pull_tpr: str, config: MembraneConfig, out_dir: str,
) -> dict[int, str]:
    """Extract per-window starting GRO files from the pulling trajectory.

    Returns a dict mapping window index → starting .gro path.
    """
    raise NotImplementedError("Phase B")
