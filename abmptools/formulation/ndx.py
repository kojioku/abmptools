# -*- coding: utf-8 -*-
"""GROMACS .ndx generation for formulation systems.

Build the index file in pure Python (no ``gmx make_ndx``) by reading
residue names from the .gro and classifying against the spec-supplied
resname table.

Groups produced:
    System         — all atoms
    Peptide        — protein residues + ACE/NME caps
    Enhancer       — all enhancer copies (charged + neutral)
    BileSalt       — all bile-salt copies
    Solvent        — TIP3P (WAT/HOH/SOL/TIP3)
    Ions           — Na+/Cl- aliases
    Peptide_Solute — Peptide ∪ Enhancer ∪ BileSalt (= "solute" for
                     2-group thermostat by default)
    Non_Peptide    — System minus Peptide
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

logger = logging.getLogger(__name__)


# Static resname tables (extend if a force field uses non-standard names).
PROTEIN_RESNAMES: Set[str] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "CYM",
    "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "ACE", "NME", "NHE",
}
SOLVENT_RESNAMES: Set[str] = {"WAT", "HOH", "SOL", "TIP3", "TIP3P", "T3P"}
ION_RESNAMES: Set[str] = {
    "NA", "NA+", "Na+", "SOD", "CL", "CL-", "Cl-", "CLA",
    "K", "K+", "POT", "MG", "MG2+", "CA", "CA2+",
}


def parse_gro_residues(gro_path: str) -> List[Tuple[int, str, str]]:
    """Parse .gro → list of (atom_index_1based, residue_name, atom_name).

    Identical layout to :func:`abmptools.membrane.parameterize_amber.parse_gro_residues`.
    """
    out: List[Tuple[int, str, str]] = []
    with open(gro_path) as f:
        f.readline()  # title
        n_atoms = int(f.readline())
        for _ in range(n_atoms):
            line = f.readline()
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            atom_idx = int(line[15:20])
            out.append((atom_idx, resname, atomname))
    return out


def classify_atoms(
    atoms: Sequence[Tuple[int, str, str]],
    *,
    enhancer_resnames: Iterable[str],
    bile_salt_resnames: Iterable[str],
    peptide_resnames: Iterable[str] = (),
) -> Dict[str, List[int]]:
    """Group atom indices (1-based) into named GROMACS index groups.

    Parameters
    ----------
    atoms
        Output of :func:`parse_gro_residues`.
    enhancer_resnames, bile_salt_resnames
        3-letter tags as declared in :class:`EnhancerSpec.resname` /
        :class:`BileSaltSpec.resname`.
    """
    enh = set(enhancer_resnames)
    bile = set(bile_salt_resnames)
    pep_extra = set(peptide_resnames)  # GAFF whole-peptide tags (e.g. "OCT")
    groups: Dict[str, List[int]] = {
        "System": [],
        "Peptide": [],
        "Enhancer": [],
        "BileSalt": [],
        "Solvent": [],
        "Ions": [],
    }
    for atom_idx, resname, _atom in atoms:
        groups["System"].append(atom_idx)
        if resname in PROTEIN_RESNAMES or resname in pep_extra:
            groups["Peptide"].append(atom_idx)
        elif resname in enh:
            groups["Enhancer"].append(atom_idx)
        elif resname in bile:
            groups["BileSalt"].append(atom_idx)
        elif resname in SOLVENT_RESNAMES:
            groups["Solvent"].append(atom_idx)
        elif resname in ION_RESNAMES:
            groups["Ions"].append(atom_idx)
        else:
            logger.debug(
                "unclassified residue %s (atom %d) — left in System only",
                resname, atom_idx,
            )
    # composite groups
    groups["Peptide_Solute"] = (
        groups["Peptide"] + groups["Enhancer"] + groups["BileSalt"]
    )
    # Non_Peptide_Solute = System − Peptide_Solute (solvent + ions only).
    # We expose this under the name ``Non_Peptide`` for the default
    # 2-group thermostat to avoid atom overlap between tc-grps.
    solute_set = set(groups["Peptide_Solute"])
    groups["Non_Peptide"] = [
        idx for idx in groups["System"] if idx not in solute_set
    ]
    return groups


def write_ndx(groups: Dict[str, List[int]], ndx_path: str) -> str:
    """Write a GROMACS .ndx for the given group mapping."""
    Path(ndx_path).parent.mkdir(parents=True, exist_ok=True)
    order = [
        "System", "Peptide", "Enhancer", "BileSalt",
        "Solvent", "Ions", "Peptide_Solute", "Non_Peptide",
    ]
    lines: List[str] = []
    for name in order:
        idxs = groups.get(name, [])
        lines.append(f"[ {name} ]")
        if not idxs:
            lines.append("")
            continue
        # 15 atoms per line, matching gmx make_ndx output
        for i in range(0, len(idxs), 15):
            lines.append(" ".join(f"{x:>5d}" for x in idxs[i : i + 15]))
        lines.append("")
    Path(ndx_path).write_text("\n".join(lines))
    logger.info("wrote index: %s", ndx_path)
    return ndx_path


def write_index_from_gro(
    *,
    gro_path: str,
    ndx_path: str,
    enhancer_resnames: Iterable[str],
    bile_salt_resnames: Iterable[str],
    peptide_resnames: Iterable[str] = (),
) -> str:
    """High-level: parse the .gro, classify, write the .ndx."""
    atoms = parse_gro_residues(gro_path)
    groups = classify_atoms(
        atoms,
        enhancer_resnames=enhancer_resnames,
        bile_salt_resnames=bile_salt_resnames,
        peptide_resnames=peptide_resnames,
    )
    return write_ndx(groups, ndx_path)
