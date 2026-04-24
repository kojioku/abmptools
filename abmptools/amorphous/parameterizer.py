# -*- coding: utf-8 -*-
"""
abmptools.amorphous.parameterizer
-----------------------------------
Parameterize a multi-component system with OpenFF and export GROMACS files
using OpenFF Interchange.

All OpenFF / OpenMM imports are lazy.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


def _check_interchange_available() -> None:
    """Raise RuntimeError if openff-interchange is not installed."""
    try:
        import openff.interchange  # noqa: F401
    except ImportError:
        raise RuntimeError(
            "openff-interchange is required for GROMACS parameterization.\n"
            "Install: conda install -c conda-forge openff-interchange"
        )


def _check_openmm_available() -> None:
    """Raise RuntimeError if OpenMM is not installed."""
    try:
        import openmm  # noqa: F401
    except ImportError:
        raise RuntimeError(
            "OpenMM is required for parameterization.\n"
            "Install: conda install -c conda-forge openmm"
        )


def create_interchange(
    molecules: List[Any],
    counts: List[int],
    box_size_nm: float,
    mixture_pdb: str,
    forcefield_name: str = "openff_unconstrained-2.1.0.offxml",
    use_precomputed_charges: bool = False,
) -> Any:
    """Build an OpenFF Interchange for a packed mixture.

    Parameters
    ----------
    molecules : list of openff.toolkit.Molecule
        One OpenFF Molecule per component (with conformer, optionally
        with pre-assigned partial charges — see *use_precomputed_charges*).
    counts : list of int
        Number of molecules for each component.
    box_size_nm : float
        Cubic box edge length [nm].
    mixture_pdb : str
        Path to the Packmol-generated mixture PDB.
    forcefield_name : str
        OpenFF force field OFFXML name.
    use_precomputed_charges : bool
        When ``True``, pass *molecules* as ``charge_from_molecules`` to
        :meth:`Interchange.from_smirnoff`, telling Interchange to reuse
        the pre-assigned ``partial_charges`` on each Molecule instead of
        invoking AM1-BCC via ``sqm``. Required for the
        ``charge_method='nagl'`` / ``'gasteiger'`` path (Windows-native
        support, since AmberTools' ``sqm`` has no Windows build).

    Returns
    -------
    openff.interchange.Interchange
        The parameterized system.
    """
    _check_interchange_available()
    _check_openmm_available()

    from openff.toolkit import ForceField, Topology, Molecule
    from openff.interchange import Interchange
    from openff.units import unit as off_unit
    import numpy as np

    ff = ForceField(forcefield_name)

    # Build OpenFF Topology from the mixture PDB with molecule templates
    topology = Topology.from_pdb(
        mixture_pdb,
        unique_molecules=molecules,
    )

    # Set box vectors
    box_nm = box_size_nm
    box_vectors = np.eye(3) * box_nm * off_unit.nanometer
    topology.box_vectors = box_vectors

    logger.info("Creating Interchange with %s (precomputed charges=%s) ...",
                forcefield_name, use_precomputed_charges)
    kwargs = {"force_field": ff, "topology": topology}
    if use_precomputed_charges:
        kwargs["charge_from_molecules"] = molecules
    interchange = Interchange.from_smirnoff(**kwargs)

    return interchange


def export_gromacs(
    interchange: Any,
    gro_path: str,
    top_path: str,
) -> Dict[str, str]:
    """Export Interchange to GROMACS .gro and .top files.

    Parameters
    ----------
    interchange : openff.interchange.Interchange
        The parameterized system.
    gro_path : str
        Output path for the .gro file.
    top_path : str
        Output path for the .top file.

    Returns
    -------
    dict
        {"gro": abs_path, "top": abs_path}
    """
    os.makedirs(os.path.dirname(gro_path) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(top_path) or ".", exist_ok=True)

    interchange.to_gro(gro_path)
    interchange.to_top(top_path)

    logger.info("Wrote GROMACS files: %s, %s", gro_path, top_path)
    return {
        "gro": str(Path(gro_path).resolve()),
        "top": str(Path(top_path).resolve()),
    }
