# -*- coding: utf-8 -*-
"""
abmptools.membrane.mdp_us_protocol
----------------------------------
MDP file writers for the membrane-US pipeline.

Stages emitted
--------------
- ``em.mdp``           — energy minimisation (steep)
- ``nvt.mdp``          — NVT heating to T (V-rescale, semiisotropic-irrelevant)
- ``npt.mdp``          — NPT relaxation, **semiisotropic** pressure coupling
- ``pull.mdp``         — pulled NPT (handed off to pulling.py for pull-code)
- ``window_<i>.mdp``   — per-window NPT US (handed off to umbrella.py)

The ``npt`` / ``pull`` / ``window`` MDPs all use semiisotropic pressure
coupling so the bilayer can independently relax in the lipid plane
versus the membrane normal:

    Pcoupl = C-rescale         (or Berendsen during equilibration)
    Pcoupltype = semiisotropic
    ref_p = 1.0 1.0            ; (xy, z) bar
    compressibility = 4.5e-5 4.5e-5

This module reuses the temperature-coupling string utilities from
:mod:`abmptools.amorphous.mdp_protocol` where applicable; semiisotropic
and pull-code blocks are membrane-specific and live here.
"""
from __future__ import annotations

import logging
from typing import Dict

from .models import MembraneConfig

logger = logging.getLogger(__name__)


def write_equilibration_mdps(
    *, config: MembraneConfig, equil_dir: str,
) -> Dict[str, str]:
    """Write em / nvt / npt MDPs for pre-pulling equilibration.

    Returns
    -------
    dict
        Keys ``"em"``, ``"nvt"``, ``"npt"`` mapped to absolute file paths.
    """
    raise NotImplementedError(
        "Phase A: equilibration MDP writers not yet implemented. "
        "Phase B will write em.mdp (steep), nvt.mdp (V-rescale T-couple), "
        "and npt.mdp (semiisotropic Pcoupl)."
    )


def render_em_mdp(*, config: MembraneConfig) -> str:
    """Return the contents of em.mdp as a single string."""
    raise NotImplementedError("Phase B")


def render_nvt_mdp(*, config: MembraneConfig) -> str:
    """Return the contents of nvt.mdp as a single string."""
    raise NotImplementedError("Phase B")


def render_npt_mdp(*, config: MembraneConfig) -> str:
    """Return the contents of npt.mdp (semiisotropic) as a single string."""
    raise NotImplementedError("Phase B")


def render_pull_block(*, config: MembraneConfig,
                      pull_init_nm: float | None = None,
                      pull_rate: float | None = None) -> str:
    """Render the [pull] / pull-coord block for membrane US.

    Common code shared between pulling.py (rate-pulling) and
    umbrella.py (per-window static pull at init = z_min + i * dz).
    """
    raise NotImplementedError("Phase B")
