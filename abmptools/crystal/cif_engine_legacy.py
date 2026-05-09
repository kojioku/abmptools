# -*- coding: utf-8 -*-
"""
abmptools.crystal.cif_engine_legacy
------------------------------------
Python API wrapper around :func:`abmptools.readcif.run_legacy_cif_pipeline`.

Phase C-2 thin adapter: lets the new :class:`CrystalOrchestrator` invoke
the legacy CIF parser by-passing the CLI argv path. Behaviour is
identical to the historical ``python -m abmptools.readcif`` invocation;
the Phase B regression fixture
(``tests/regression/reference/main/crystal_csp7/``) covers byte-
equivalence.

The underlying parser is the self-rolled implementation in
:mod:`abmptools.readcif` (982 lines). This adapter does **not** change
its symmetry-operation table or layer-count handling; users hitting
unsupported space groups should switch to ``--engine ase`` (Phase C-3).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import List, Optional

from abmptools.readcif import run_legacy_cif_pipeline


def run_legacy(
    cif: str,
    *,
    layer: int = 5,
    atoms_in_mol: Optional[List[int]] = None,
    odir: str = "cifout",
    asymmetric_only: bool = False,
    write_pdb: bool = True,
    cwd: Optional[str] = None,
) -> Path:
    """Run the legacy readcif pipeline programmatically.

    Parameters
    ----------
    cif
        Path to the CIF file (absolute, or relative to *cwd*).
    layer
        Supercell expansion layer (``layer=N`` -> :math:`(2N+1)^3` cells).
    atoms_in_mol
        Atoms per asymmetric-unit molecule (``[32]`` for Z'=1, ``[32,32]``
        for Z'=2). Defaults to ``[32]`` to match the csp7 demo.
    odir
        Output directory (``readcif`` writes
        ``<odir>/layer<L>/{pdb,xyz}/``).
    asymmetric_only
        If True, emit only the asymmetric unit (skip layer expansion).
    write_pdb
        If True (default), readcif emits PDB files alongside XYZ.
    cwd
        Working directory; defaults to current directory. The legacy
        parser writes into ``odir`` relative to the **current** dir,
        so callers that need a sandbox should ``os.chdir`` (the
        orchestrator does this internally).

    Returns
    -------
    Path
        Directory containing the generated ``layer<L>/pdb/*.pdb`` and
        ``layer<L>/xyz/*.xyz`` files.
    """
    if atoms_in_mol is None:
        atoms_in_mol = [32]

    args = argparse.Namespace(
        input=[[cif]],
        atomnum=[atoms_in_mol],
        layer=layer,
        odir=odir,
        asymonly=asymmetric_only,
        nopdb=write_pdb,        # readcif uses store_false; True means "emit PDB"
        noout=False,
        calcdist=False,
        intra=False,
        dist=2.0,
        tgtatom=[],
        min=False,
    )

    if cwd is not None:
        prev = os.getcwd()
        os.chdir(cwd)
        try:
            run_legacy_cif_pipeline(args)
        finally:
            os.chdir(prev)
        out_root = Path(cwd) / odir
    else:
        run_legacy_cif_pipeline(args)
        out_root = Path(odir)

    return out_root / f"layer{layer}"


__all__ = ["run_legacy"]
