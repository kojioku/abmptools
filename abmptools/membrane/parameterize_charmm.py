# -*- coding: utf-8 -*-
"""
abmptools.membrane.parameterize_charmm
--------------------------------------
CHARMM36 backend: CHARMM36m + CHARMM36 lipids → GROMACS top/gro.

Pipeline
--------
1. Stage the CHARMM36 force-field directory in the build dir
   (``forcefield.itp``, ``ffnonbonded.itp``, lipid/protein .rtp, …).
   Source: Klauda lab's GROMACS port of CHARMM36
   (charmm36-jul2022.ff or newer). Path is supplied via
   :attr:`MembraneConfig.charmm_ff_dir` and **must already be present
   on disk** — this package does not download it (network pull would be
   a license-policy decision; left to the user / system admin).

2. Run ``gmx pdb2gmx`` on the peptide+bilayer PDB with
   ``-ff charmm36-jul2022 -water tip3p`` to obtain the topology.
   Note: pdb2gmx handles standard residues (proteins, CHARMM36 lipids,
   ions, water) without CGenFF.

3. Build a single ``system.top`` that ``#include``s the staged
   force-field directory's master itp.

4. Generate index groups (Bilayer / Peptide / Solvent / Ions / System).

License
-------
Uses only **CHARMM36 force-field parameter values** (free for academic
and industrial use per MacKerell lab). Does **NOT**:
- call CGenFF web server (commercial use forbidden by Silcsbio license)
- consume CHARMM-GUI generated files (commercial use requires subscription)

Therefore peptide + standard lipid systems are fully commercial-OK on
this route, provided no novel small molecules requiring CGenFF appear.
For novel small molecules, switch to backend="amber" and use GAFF2.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Tuple

from .models import MembraneConfig

logger = logging.getLogger(__name__)


def parameterize(
    *, config: MembraneConfig, input_pdb: str, build_dir: str,
) -> Tuple[str, str, str]:
    """Run CHARMM36 parameterisation and emit GROMACS files.

    Parameters
    ----------
    config : MembraneConfig
        Build configuration. Must have ``backend="charmm36"`` and
        a non-empty ``charmm_ff_dir``.
    input_pdb : str
        Bilayer + peptide PDB from stage 1.
    build_dir : str
        Directory to stage the force field, run pdb2gmx, and write
        the GROMACS files.

    Returns
    -------
    (top, gro, ndx) : tuple[str, str, str]
        Absolute paths to system.top, system.gro, system.ndx.
    """
    raise NotImplementedError(
        "Phase A: CHARMM36 parameterisation not yet implemented. "
        "Phase C will: (a) stage charmm_ff_dir into build_dir, "
        "(b) run gmx pdb2gmx with -ff charmm36 -water tip3p, "
        "(c) post-process top to remove pdb2gmx-injected ifdef blocks, "
        "(d) generate index groups."
    )


def stage_forcefield(*, charmm_ff_dir: str, build_dir: str) -> str:
    """Copy / symlink the CHARMM36 .ff directory into build_dir.

    Returns the in-build path so pdb2gmx finds it via ``-ff``.
    """
    raise NotImplementedError("Phase C")


def run_pdb2gmx(
    *, input_pdb: str, ff_name: str, build_dir: str, gmx_path: str = "gmx",
) -> Tuple[str, str]:
    """Run ``gmx pdb2gmx`` and return (top_path, gro_path)."""
    raise NotImplementedError("Phase C")
