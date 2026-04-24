# -*- coding: utf-8 -*-
"""
abmptools.amorphous.molecule_prep
-----------------------------------
SMILES / SDF → OpenFF Molecule → single-molecule PDB.

All OpenFF / RDKit imports are lazy.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Tuple

logger = logging.getLogger(__name__)


def _check_openff_available() -> None:
    """Raise RuntimeError if openff-toolkit is not installed."""
    try:
        import openff.toolkit  # noqa: F401
    except ImportError:
        raise RuntimeError(
            "openff-toolkit is required for molecule preparation.\n"
            "Install: conda install -c conda-forge openff-toolkit"
        )


def _load_molecule_from_smiles(smiles: str, name: str = "") -> Any:
    """Create an OpenFF Molecule from a SMILES string."""
    _check_openff_available()
    from openff.toolkit import Molecule

    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    mol.generate_conformers(n_conformers=1)
    if name:
        mol.name = name
    return mol


def _load_molecule_from_sdf(sdf_path: str, name: str = "") -> Any:
    """Create an OpenFF Molecule from an SDF file."""
    _check_openff_available()
    from openff.toolkit import Molecule

    mol = Molecule.from_file(sdf_path, allow_undefined_stereo=True)
    if mol.n_conformers == 0:
        mol.generate_conformers(n_conformers=1)
    if name:
        mol.name = name
    return mol


def prepare_molecule(
    smiles: str = "",
    sdf_path: str = "",
    name: str = "",
    charge_method: str = "",
    nagl_model: str = "openff-gnn-am1bcc-0.1.0-rc.3.pt",
) -> Any:
    """Load or create an OpenFF Molecule from SMILES or SDF.

    Parameters
    ----------
    smiles : str
        SMILES string (mutually exclusive with *sdf_path*).
    sdf_path : str
        Path to an SDF/MOL file (mutually exclusive with *smiles*).
    name : str
        Optional molecule name.
    charge_method : str
        ``""`` / ``"am1bcc"`` — leave charges unassigned so
        :class:`openff.interchange.Interchange` invokes AM1-BCC via
        AmberTools ``sqm`` (Linux/macOS default, fails on Windows).
        ``"nagl"``            — pre-assign using openff-nagl's ML
        AM1-BCC approximation (Windows-compatible).
        ``"gasteiger"``       — pre-assign Gasteiger charges
        (fast, lower fidelity; useful as a sanity check).
    nagl_model : str
        OpenFF NAGL model identifier (used only when
        ``charge_method='nagl'``).

    Returns
    -------
    openff.toolkit.Molecule
        The prepared molecule (with at least one conformer, and partial
        charges if *charge_method* requested pre-assignment).
    """
    if smiles and sdf_path:
        raise ValueError("Provide either 'smiles' or 'sdf_path', not both.")
    if not smiles and not sdf_path:
        raise ValueError("Either 'smiles' or 'sdf_path' must be provided.")

    if smiles:
        mol = _load_molecule_from_smiles(smiles, name)
    else:
        mol = _load_molecule_from_sdf(sdf_path, name)

    if not name:
        mol.name = mol.name or "MOL"

    _maybe_assign_partial_charges(mol, charge_method, nagl_model)

    logger.info("Prepared molecule '%s': %d atoms, MW=%.2f g/mol",
                mol.name, mol.n_atoms,
                _molecular_weight(mol))
    return mol


def _maybe_assign_partial_charges(mol: Any, method: str, nagl_model: str) -> None:
    """Pre-assign partial charges on *mol* when *method* is not AM1-BCC.

    AM1-BCC is handled downstream by Interchange itself (when we pass the
    topology without `charge_from_molecules`), so we deliberately skip it
    here. Anything else is assigned up-front on the Molecule object; the
    caller (Interchange) must then be told to use those charges via
    ``charge_from_molecules=[...]``.
    """
    if not method or method == "am1bcc":
        return
    if method == "nagl":
        try:
            mol.assign_partial_charges(partial_charge_method=nagl_model)
        except Exception as e:
            raise RuntimeError(
                f"Failed to assign NAGL charges using model '{nagl_model}'. "
                "Ensure openff-nagl + openff-nagl-models are installed "
                "(conda install -c conda-forge openff-nagl openff-nagl-models). "
                f"Original error: {e}"
            ) from e
        return
    if method == "gasteiger":
        mol.assign_partial_charges(partial_charge_method="gasteiger")
        return
    raise ValueError(
        f"Unknown charge_method '{method}'. "
        "Choose from: '', 'am1bcc', 'nagl', 'gasteiger'."
    )


def _molecular_weight(mol: Any) -> float:
    """Compute molecular weight [g/mol] from an OpenFF Molecule."""
    from openff.units import unit as off_unit
    try:
        # openff-toolkit >= 0.15
        mw = sum(
            atom.mass.m_as(off_unit.dalton) for atom in mol.atoms
        )
    except (AttributeError, TypeError):
        # fallback: manual element mass lookup
        from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
        mw = sum(atom.mass.magnitude for atom in mol.atoms)
    return float(mw)


def get_molecular_weight(mol: Any) -> float:
    """Return molecular weight [g/mol] for an OpenFF Molecule."""
    return _molecular_weight(mol)


def write_single_mol_pdb(mol: Any, output_path: str) -> str:
    """Write a single-molecule PDB from an OpenFF Molecule.

    Parameters
    ----------
    mol : openff.toolkit.Molecule
        Molecule with at least one conformer.
    output_path : str
        Path to write the PDB file.

    Returns
    -------
    str
        Absolute path to the written PDB.
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    mol.to_file(output_path, file_format="pdb")
    return str(Path(output_path).resolve())
