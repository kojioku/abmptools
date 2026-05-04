# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.system_packer
------------------------------------
GROMACS subprocess wrappers + Python-only index writer for the CG
peptide-membrane builder.

Most wrappers reuse :mod:`abmptools.cg.peptide.system_packer` directly
(``add_ions`` is FF-agnostic). The only new functionality is the
membrane-specific index file writer that produces named groups
(``Bilayer`` / ``Peptide`` / ``W`` / ``NA`` / ``CL`` / ``Non_Bilayer``)
purely from the ``.gro`` -- no ``gmx make_ndx`` stdin choreography.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple

# Re-export AA membrane's parser (gro residue parsing is FF-agnostic).
from abmptools.membrane.parameterize_amber import parse_gro_residues

# Re-export cg.peptide.add_ions; v1 of the membrane build inherits the
# same NA/CL ion convention and W solvent group.
from abmptools.cg.peptide.system_packer import (  # noqa: F401  (re-export)
    PackedSystem,
    add_ions as add_ions_cg,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Residue classification (CG-flavoured)
# ---------------------------------------------------------------------------

#: Martini 3 lipid residue names. v1 supports POPC explicitly; the rest are
#: included so that future mixture builds work when the dataclass relaxes
#: the single-species restriction.
CG_LIPID_RESNAMES = frozenset({
    "POPC", "DOPC", "POPE", "DOPE", "POPG", "DPPC", "DSPC", "DLPC",
    "POPA", "POPS", "DOPS", "DAPC", "POPI",
    # Cholesterol family (placeholders; M3 chol parameters published 2024+)
    "CHOL", "CHL1",
})

#: Standard 20 amino-acid 3-letter codes (Martini residue names match
#: full atomistic 3-letter; insane preserves them from the input PDB).
CG_PROTEIN_RESNAMES = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
    # Histidine variants (rare in M3 PDBs but harmless to recognise)
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
})

#: Martini 3 solvent bead names.
CG_SOLVENT_RESNAMES = frozenset({
    "W",        # Martini 3 standard
    "WN",       # negative-charge water variant (rarely used as solvent)
    "WP",       # positive-charge variant
})

#: Ion residue names. After :func:`topology_composer.normalize_ion_atom_names_gro`,
#: ``NA+``/``CL-`` are converted to ``NA``/``CL``; both forms are listed
#: here so that an un-normalized .gro is still classified correctly.
CG_NA_RESNAMES = frozenset({"NA", "NA+"})
CG_CL_RESNAMES = frozenset({"CL", "CL-"})


def classify_residue_cg(resname: str) -> str:
    """Classify a CG residue name.

    Returns one of: ``'lipid'``, ``'peptide'``, ``'solvent'``,
    ``'na'``, ``'cl'``, ``'other'``.
    """
    r = resname.strip().upper()
    if r in CG_LIPID_RESNAMES:
        return "lipid"
    if r in CG_PROTEIN_RESNAMES:
        return "peptide"
    if r in CG_SOLVENT_RESNAMES:
        return "solvent"
    if r in CG_NA_RESNAMES:
        return "na"
    if r in CG_CL_RESNAMES:
        return "cl"
    return "other"


# ---------------------------------------------------------------------------
# Index writer (Python-only, no gmx make_ndx)
# ---------------------------------------------------------------------------

def write_ndx_from_gro_cg(
    *, gro_path: str, ndx_path: str,
) -> str:
    """Write a GROMACS .ndx with CG-membrane groups derived from residue names.

    Groups (in this order, omitted if empty):

    - ``System``       -- all atoms
    - ``Bilayer``      -- POPC/DOPC/... (CG_LIPID_RESNAMES)
    - ``Peptide``      -- 20 standard amino-acid 3-letter codes
    - ``W``            -- Martini W solvent bead
    - ``NA``           -- sodium (post-normalization)
    - ``CL``           -- chloride (post-normalization)
    - ``Non_Bilayer``  -- System minus Bilayer (= Peptide + W + NA + CL +
                          any unclassified). Used for 2-group thermostatting.

    Returns the absolute path to the written .ndx.
    """
    atoms: List[Tuple[int, str, str]] = parse_gro_residues(gro_path)

    groups: Dict[str, List[int]] = {
        "System":      [],
        "Bilayer":     [],
        "Peptide":     [],
        "W":           [],
        "NA":          [],
        "CL":          [],
        "Non_Bilayer": [],
    }
    other_residues: set = set()
    for idx, resname, _atomname in atoms:
        groups["System"].append(idx)
        cls = classify_residue_cg(resname)
        if cls == "lipid":
            groups["Bilayer"].append(idx)
        else:
            groups["Non_Bilayer"].append(idx)

        if cls == "peptide":
            groups["Peptide"].append(idx)
        elif cls == "solvent":
            groups["W"].append(idx)
        elif cls == "na":
            groups["NA"].append(idx)
        elif cls == "cl":
            groups["CL"].append(idx)
        elif cls != "lipid":
            other_residues.add(resname)

    if other_residues:
        logger.warning(
            "write_ndx_from_gro_cg: unclassified residues land in System "
            "and Non_Bilayer only: %s", sorted(other_residues),
        )

    out_lines: List[str] = []
    for name, indices in groups.items():
        if not indices:
            continue
        out_lines.append(f"[ {name} ]")
        # GROMACS convention: 15 indices per line, right-aligned width 5.
        for i in range(0, len(indices), 15):
            chunk = indices[i:i + 15]
            out_lines.append(" ".join(f"{j:>5d}" for j in chunk))
        out_lines.append("")

    Path(ndx_path).write_text("\n".join(out_lines))
    return str(Path(ndx_path).resolve())
