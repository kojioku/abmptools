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
    pdb_path: str = "",
    name: str = "",
) -> Any:
    """Load or create an OpenFF Molecule from SMILES, SDF, or PDB.

    Parameters
    ----------
    smiles : str
        SMILES string (mutually exclusive with the other two).
    sdf_path : str
        Path to an SDF/MOL file.
    pdb_path : str
        Path to a PDB file (Phase 9-a; use for pre-built oligomer /
        polymer structures). The OpenFF Toolkit's
        ``Molecule.from_polymer_pdb`` is tried first so chains that
        match its supported chemistry get proper bond perception;
        otherwise ``Molecule.from_file`` (legacy path) is the fallback.
        OpenFF Toolkit's polymer support is improving but still limited
        for arbitrary chemistries — if FF assignment fails downstream,
        either supply a SMILES-derived monomer instead or use the legacy
        UDF route.
    name : str
        Optional molecule name.

    Returns
    -------
    openff.toolkit.Molecule
        The prepared molecule (with at least one conformer).
    """
    sources = [bool(smiles), bool(sdf_path), bool(pdb_path)]
    n_sources = sum(sources)
    if n_sources == 0:
        raise ValueError("One of 'smiles', 'sdf_path', or 'pdb_path' must be provided.")
    if n_sources > 1:
        raise ValueError("Provide exactly one of 'smiles', 'sdf_path', or 'pdb_path'.")

    if smiles:
        mol = _load_molecule_from_smiles(smiles, name)
    elif sdf_path:
        mol = _load_molecule_from_sdf(sdf_path, name)
    else:
        mol = _load_molecule_from_pdb(pdb_path, name)

    if not name:
        mol.name = mol.name or "MOL"
    logger.info("Prepared molecule '%s': %d atoms, MW=%.2f g/mol",
                mol.name, mol.n_atoms,
                _molecular_weight(mol))
    return mol


def _load_molecule_from_pdb(pdb_path: str, name: str = "") -> Any:
    """Load an OpenFF Molecule from a PDB file (oligomer / polymer).

    OpenFF Toolkit's PDB support is uneven:

    1. ``Molecule.from_polymer_pdb`` works for biopolymers covered by
       the bundled substructure library (proteins, DNA), not for
       generic chain oligomers like propane×3.
    2. ``Molecule.from_file`` for ``.pdb`` raises NotImplementedError
       in current releases — RDKit can't safely infer bond orders
       from a PDB alone.
    3. ``Molecule.from_pdb_and_smiles`` works but needs an oligomer
       SMILES the caller would have to construct.

    The reliable path for generic oligomer PDBs is OpenBabel: it
    perceives bonds from CONECT + chemistry, writes an SDF with the
    *same atom order*, and OpenFF Toolkit reads SDF natively. We
    therefore try (1) first, then route through SDF.
    """
    from openff.toolkit import Molecule

    # 1. Polymer-aware loader (biopolymer chemistry only).
    try:
        mol = Molecule.from_polymer_pdb(pdb_path)
        if name:
            mol.name = name
        return mol
    except Exception:
        pass

    # 2. OpenBabel PDB → SDF round-trip preserves atom order and adds
    #    bond-order information that OpenFF Toolkit needs.
    sdf_path = _pdb_to_sdf_via_openbabel(pdb_path)
    try:
        mol = _load_molecule_from_sdf(sdf_path, name)
    except Exception as e:
        raise RuntimeError(
            f"OpenFF Toolkit could not load '{pdb_path}' even after "
            f"OpenBabel-assisted SDF conversion ({sdf_path}). "
            "For complex polymers, fall back to the legacy UDF route "
            "(backend='udf')."
        ) from e
    return mol


def _pdb_to_sdf_via_openbabel(pdb_path: str) -> str:
    """Convert ``pdb_path`` to ``<pdb_path>.sdf`` using OpenBabel.

    OpenBabel preserves the input atom order and adds inferred
    bond-order information so OpenFF Toolkit can ingest the result.
    Raises ``RuntimeError`` if the ``obabel`` binary is missing or
    fails.
    """
    import shutil
    import subprocess
    from pathlib import Path

    if shutil.which("obabel") is None:
        raise RuntimeError(
            "obabel (OpenBabel CLI) was not found on PATH. It is "
            "required to convert oligomer PDBs into a form OpenFF "
            "Toolkit can ingest. Install via "
            "'micromamba install -c conda-forge openbabel'."
        )

    pdb = Path(pdb_path)
    sdf = pdb.with_suffix(".sdf")
    cmd = ["obabel", str(pdb), "-O", str(sdf)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"obabel failed to convert {pdb_path} → {sdf}: "
            f"{e.stderr.strip() or e.stdout.strip() or 'no diagnostic'}"
        ) from e
    if not sdf.is_file():
        raise RuntimeError(
            f"obabel returned 0 but {sdf} was not produced."
        )
    return str(sdf)


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
