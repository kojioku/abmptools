# -*- coding: utf-8 -*-
"""
abmptools.membrane.bilayer
--------------------------
Bilayer + peptide construction via *packmol-memgen* (AmberTools).

Why packmol-memgen
------------------
packmol-memgen (Schott-Verdugo & Gohlke, 2019) is the most convenient
license-OK route for assembling a bilayer + water + ions + peptide PDB
without CHARMM-GUI:

- Bundled with AmberTools — fully commercial OK
- Knows AMBER Lipid21 / Lipid17 residue conventions out-of-the-box
  (head-tail split for POPC/POPE/POPG/etc.)
- Optional ``--charmm`` flag for CHARMM-formatted output (see Phase C)
- Handles peptide insertion via ``--solute <pep.pdb>``

Output is a single PDB with ``LIP``, ``WAT``, ``Na+``/``Cl-``, and
peptide residues. Stage 2 (parameterize_amber / parameterize_charmm)
consumes this PDB directly.
"""
from __future__ import annotations

import logging
from pathlib import Path

from .models import MembraneConfig

logger = logging.getLogger(__name__)


def build_bilayer(*, config: MembraneConfig, build_dir: str) -> str:
    """Build a bilayer (+ optional peptide) PDB and return its path.

    Parameters
    ----------
    config : MembraneConfig
        Build configuration.
    build_dir : str
        Directory in which to invoke packmol-memgen.

    Returns
    -------
    str
        Absolute path to the assembled PDB
        (``{build_dir}/bilayer_peptide.pdb``).
    """
    raise NotImplementedError(
        "Phase A: bilayer construction not yet implemented. "
        "Phase B will assemble the packmol-memgen command line from "
        "config.lipids / config.peptide / config.box_xy_nm / "
        "config.water_thickness_nm and run it."
    )


def write_peptide_pdb(*, config: MembraneConfig, out_path: str) -> str:
    """Build an initial peptide PDB from sequence (or copy pdb_path).

    Returns the absolute output path.
    """
    raise NotImplementedError("Phase B")


def assemble_packmol_memgen_cmd(*, config: MembraneConfig,
                                peptide_pdb: str, output_pdb: str) -> list[str]:
    """Compose the packmol-memgen argv from a config object.

    Pure function for unit testing — no side effects.
    """
    raise NotImplementedError("Phase B")
