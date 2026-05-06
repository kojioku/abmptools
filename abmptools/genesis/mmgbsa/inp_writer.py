# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.inp_writer
-----------------------------------
Render GENESIS ``atdyn`` ``.inp`` control files for the three MM/GBSA
systems (``complex`` / ``ligand`` / ``receptor``).

Mirrors POC ``3_rungbsa.py`` ``content_template``: ``[ENERGY]
implicit_solvent=GBSA`` + NOBC + ``CUTOFF`` electrostatic + 1-step
``[MINIMIZE]`` for single-point evaluation.

Each renderer returns a string; the orchestrator writes it to disk
(``<name>.inp`` in the per-target build dir).
"""
from __future__ import annotations

from textwrap import dedent
from typing import List

from .models import EnergyProtocol, MinimizationProtocol


SYSTEM_NAMES: List[str] = ["complex", "ligand", "receptor"]


def render_gbsa_inp(
    name: str,
    energy: EnergyProtocol,
    minimize: MinimizationProtocol,
) -> str:
    """Render an atdyn ``.inp`` for ``<name>`` (= complex/ligand/receptor).

    File expects ``<name>.prmtop`` + ``<name>.inpcrd`` to be in cwd at
    run time, and writes ``<name>.dcd`` + ``<name>.rst``.

    POC defaults: ``electrostatic=CUTOFF`` + ``cutoffdist=99.9`` +
    ``implicit_solvent=GBSA`` + ``MINIMIZE method=SD nsteps=1``.
    """
    if name not in SYSTEM_NAMES:
        raise ValueError(
            f"name must be one of {SYSTEM_NAMES}, got {name!r}"
        )

    # Render GBSA-specific lines only when implicit_solvent=GBSA.
    if energy.implicit_solvent == "GBSA":
        gbsa_lines = (
            f"implicit_solvent = GBSA\n"
            f"gbsa_salt_cons   = {energy.gbsa_salt_cons:.4f}\n"
            f"gbsa_surf_tens   = {energy.gbsa_surf_tens:.6f}\n"
            f"gbsa_vdw_offset  = {energy.gbsa_vdw_offset:.4f}"
        )
    else:
        gbsa_lines = "implicit_solvent = NONE"

    return dedent(f"""\
        [INPUT]
        prmtopfile = {name}.prmtop   # AMBER parameter topology file
        ambcrdfile = {name}.inpcrd   # AMBER coordinates file

        [OUTPUT]
        dcdfile = {name}.dcd          # DCD trajectory file
        rstfile = {name}.rst          # restart file

        [ENERGY]
        forcefield       = AMBER
        electrostatic    = {energy.electrostatic}
        switchdist       = {energy.switchdist_A:.4f}
        cutoffdist       = {energy.cutoffdist_A:.4f}
        pairlistdist     = {energy.pairlistdist_A:.4f}
        {gbsa_lines}

        [MINIMIZE]
        method           = {minimize.method}
        nsteps           = {minimize.nsteps}
        eneout_period    = {minimize.eneout_period}
        crdout_period    = {minimize.crdout_period}
        rstout_period    = {minimize.rstout_period}
        nbupdate_period  = {minimize.nbupdate_period}

        [BOUNDARY]
        type             = NOBC      # no periodic boundary; required for GBSA
        """).rstrip() + "\n"


def render_three_inps(
    energy: EnergyProtocol,
    minimize: MinimizationProtocol,
) -> dict:
    """Render all three ``.inp`` files.

    Returns ``{"complex": str, "ligand": str, "receptor": str}``.
    """
    return {name: render_gbsa_inp(name, energy, minimize) for name in SYSTEM_NAMES}


__all__ = [
    "SYSTEM_NAMES",
    "render_gbsa_inp",
    "render_three_inps",
]
