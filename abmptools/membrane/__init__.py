# -*- coding: utf-8 -*-
"""
abmptools.membrane
==================
Lipid bilayer + peptide umbrella sampling (US) builder for GROMACS.

Companion to :mod:`abmptools.amorphous`. Where ``amorphous`` builds 3D
random boxes via OpenFF + Packmol, ``membrane`` builds 2D-periodic
bilayer systems for permeability PMF calculations.

Pipeline overview
-----------------
1. Bilayer construction (lipids + water + ions + peptide)
2. Force-field parameterisation (AMBER or CHARMM36 backend)
3. Equilibration MDP series (em → nvt → npt-semiisotropic)
4. Pull MDP for reaction-coordinate generation (z-axis)
5. Umbrella-sampling window MDP series
6. PMF analysis (gmx wham / pymbar)

Commercial-use license rule
---------------------------
This subpackage is designed to produce systems usable for **commercial
research** without separate licensing. To preserve that property,
the following are **NOT permitted**:

- CGenFF web server output (academic-only; commercial requires Silcsbio)
- CHARMM-GUI automated topology / system generation (academic-only)

The following **are** permitted:

- AMBER force fields (ff19SB, ff14SB, Lipid21, Lipid17, GAFF/GAFF2)
  bundled with AmberTools — fully free, including commercial use.
- CHARMM36 force-field **parameter values** sourced from the MacKerell
  group's free download or Klauda lab's GROMACS port — commercial OK
  per MacKerell lab's stated license ("free of charge to academic and
  industrial researchers"), provided no CGenFF / CHARMM-GUI services
  are used.
- Packmol / packmol-memgen (AmberTools) — free, commercial OK.
- GROMACS (LGPL).
"""
from __future__ import annotations

from .models import (
    MembraneConfig,
    LipidSpec,
    PeptideSpec,
    IonSpec,
    USProtocol,
    Backend,
)
from .builder import MembraneUSBuilder

__all__ = [
    "MembraneConfig",
    "LipidSpec",
    "PeptideSpec",
    "IonSpec",
    "USProtocol",
    "Backend",
    "MembraneUSBuilder",
]
