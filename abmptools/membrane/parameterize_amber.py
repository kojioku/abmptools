# -*- coding: utf-8 -*-
"""
abmptools.membrane.parameterize_amber
-------------------------------------
AMBER backend: ff19SB + Lipid21 + TIP3P → GROMACS top/gro via parmed.

Pipeline
--------
1. Run *tleap* on the bilayer+peptide PDB:
     - load force fields (leaprc.protein.ff19SB, leaprc.lipid21,
       leaprc.water.tip3p, ions1lm_iod_tip3p)
     - solvate (already done by packmol-memgen, so just charge-balance)
     - save AMBER prmtop + inpcrd
2. Convert to GROMACS top/gro via *parmed*::

     parm = parmed.load_file(prmtop, inpcrd)
     parm.save("system.top", format="gromacs")
     parm.save("system.gro")

3. Generate index (.ndx) groups: ``Bilayer``, ``Peptide``, ``Solvent``,
   ``Ions``, ``System``.

License
-------
Pure-AMBER route. AmberTools is free for academic and commercial use,
including redistribution. No CGenFF, no CHARMM-GUI.
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
    """Run AMBER parameterisation and emit GROMACS files.

    Parameters
    ----------
    config : MembraneConfig
        Build configuration.
    input_pdb : str
        Bilayer + peptide PDB from stage 1 (packmol-memgen output).
    build_dir : str
        Directory to write tleap inputs/outputs and the GROMACS files.

    Returns
    -------
    (top, gro, ndx) : tuple[str, str, str]
        Absolute paths to system.top, system.gro, system.ndx.
    """
    raise NotImplementedError(
        "Phase A: AMBER parameterisation not yet implemented. "
        "Phase B will: (a) write tleap input, (b) run tleap, "
        "(c) parmed-convert prmtop/inpcrd to GROMACS, "
        "(d) generate index groups."
    )


def write_tleap_input(
    *, config: MembraneConfig, input_pdb: str, build_dir: str,
) -> str:
    """Write a tleap input script to *build_dir*/tleap.in.

    Returns the path to the script.
    """
    raise NotImplementedError("Phase B")


def run_tleap(*, tleap_input: str, tleap_path: str = "tleap") -> Tuple[str, str]:
    """Run tleap, returning (prmtop_path, inpcrd_path)."""
    raise NotImplementedError("Phase B")


def amber_to_gromacs(
    *, prmtop: str, inpcrd: str, top_out: str, gro_out: str,
) -> None:
    """Convert AMBER prmtop/inpcrd to GROMACS top/gro via parmed."""
    raise NotImplementedError("Phase B")
