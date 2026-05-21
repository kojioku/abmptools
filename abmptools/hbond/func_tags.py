"""
func_tags.py
------------
Force-field-agnostic functional atom tagging.

Different force fields use different atom-type names for the same chemical
role (e.g. carbonyl C is ``c`` in GAFF2, ``C`` in OPLS-AA, ``C`` in CHARMM,
``c3``/``c1`` etc.). To make H-bond detection FF-independent, we introduce
**functional tags** -- canonical labels for chemical roles:

    carbonyl_C    sp2 C in C=O (carboxyl, amide, ester, ketone, urea, ...)
    aromatic_C    aromatic sp2 C
    sp3_C         sp3 C (aliphatic)
    sp2_C         other sp2 C (alkene)
    carbonyl_O    sp2 O in C=O
    hydroxyl_O    sp3 O in O-H (carboxyl OH, alcohol)
    ether_O       sp3 O in ether/ester (no H)
    hydroxyl_H    H on hydroxyl_O (donor H)
    amide_N       sp2 N in amide (or amine sp2)
    amine_N       sp3 N in amine
    amide_H       H on amide_N (donor H, only present in 1°/2° amide)
    aromatic_H    H on aromatic C
    aliphatic_H   H on sp3 C
    halogen_F/Cl/Br/I
    sulfide_S
    phosphate_P

Built-in mappings: GAFF2, OPLS-AA-2005, CHARMM36, OpenFF-Sage-2.

For atom types not in the built-in table, users can extend at runtime:
    GAFF2.add_mapping("my_custom_type", "carbonyl_C")

Or rely on the **element + bond-graph fallback** (``fallback_tag_by_element``)
which kicks in automatically for any atom whose ``atom_type`` is unknown
(``None`` or a per-atom unique name like OpenFF SMIRNOFF's ``MOL0_0``). This
makes the package usable with OpenFF Sage trajectories without needing an
extra antechamber pass to assign GAFF type names.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set


# Canonical functional tags. Used as keys throughout the package.
TAG_CARBONYL_C = "carbonyl_C"
TAG_AROMATIC_C = "aromatic_C"
TAG_SP3_C = "sp3_C"
TAG_SP2_C = "sp2_C"
TAG_CARBONYL_O = "carbonyl_O"
TAG_HYDROXYL_O = "hydroxyl_O"
TAG_ETHER_O = "ether_O"
TAG_HYDROXYL_H = "hydroxyl_H"
TAG_AMIDE_N = "amide_N"
TAG_AMINE_N = "amine_N"
TAG_AMIDE_H = "amide_H"
TAG_AROMATIC_H = "aromatic_H"
TAG_ALIPHATIC_H = "aliphatic_H"
TAG_HALOGEN_F = "halogen_F"
TAG_HALOGEN_CL = "halogen_Cl"
TAG_HALOGEN_BR = "halogen_Br"
TAG_HALOGEN_I = "halogen_I"
TAG_SULFIDE_S = "sulfide_S"
TAG_PHOSPHATE_P = "phosphate_P"

ALL_TAGS: Set[str] = {
    TAG_CARBONYL_C, TAG_AROMATIC_C, TAG_SP3_C, TAG_SP2_C,
    TAG_CARBONYL_O, TAG_HYDROXYL_O, TAG_ETHER_O,
    TAG_HYDROXYL_H, TAG_AMIDE_N, TAG_AMINE_N, TAG_AMIDE_H,
    TAG_AROMATIC_H, TAG_ALIPHATIC_H,
    TAG_HALOGEN_F, TAG_HALOGEN_CL, TAG_HALOGEN_BR, TAG_HALOGEN_I,
    TAG_SULFIDE_S, TAG_PHOSPHATE_P,
}


@dataclass
class FunctionalTagMapping:
    """Atom-type to functional-tag mapping for a single force field."""
    force_field: str
    type_to_tag: Dict[str, str] = field(default_factory=dict)
    description: str = ""

    def get_tag(self, atom_type: str) -> Optional[str]:
        return self.type_to_tag.get(atom_type)

    def add_mapping(self, atom_type: str, tag: str) -> None:
        if tag not in ALL_TAGS:
            raise ValueError(f"Unknown tag {tag!r}. Must be one of {sorted(ALL_TAGS)}")
        self.type_to_tag[atom_type] = tag

    def known_types(self) -> Set[str]:
        return set(self.type_to_tag)


# ============================================================================
# Built-in mappings
# ============================================================================

GAFF2 = FunctionalTagMapping(
    force_field="GAFF2",
    description="Wang et al. (2004) generalized amber FF v2",
    type_to_tag={
        # Carbons
        "c":  TAG_CARBONYL_C,    # sp2 C in C=O (carboxyl/amide/ester/ketone)
        "c1": TAG_SP2_C,         # sp1 alkyne (lump into sp2-like)
        "c2": TAG_SP2_C,         # sp2 alkene
        "c3": TAG_SP3_C,
        "ca": TAG_AROMATIC_C,
        "cc": TAG_AROMATIC_C,    # aromatic C in pure ring
        "cd": TAG_AROMATIC_C,
        "ce": TAG_SP2_C,
        "cf": TAG_SP2_C,
        # Oxygens
        "o":  TAG_CARBONYL_O,    # sp2 O in C=O
        "oh": TAG_HYDROXYL_O,    # O in O-H (carboxyl OH, alcohol)
        "os": TAG_ETHER_O,       # ester/ether O (no H)
        # Nitrogens
        "n":  TAG_AMIDE_N,       # amide N (sp2)
        "n2": TAG_AMIDE_N,       # double-bonded N
        "n3": TAG_AMINE_N,       # sp3 amine N
        "n4": TAG_AMINE_N,       # quaternary
        "na": TAG_AMIDE_N,       # N in aromatic ring with H
        "nb": TAG_AROMATIC_C,    # N in aromatic ring no H (treat as aromatic)
        "nh": TAG_AMIDE_N,       # amine adjacent to aromatic
        # Hydrogens
        "ho": TAG_HYDROXYL_H,    # H on hydroxyl_O
        "hn": TAG_AMIDE_H,       # H on N
        "ha": TAG_AROMATIC_H,
        "hc": TAG_ALIPHATIC_H,
        "h1": TAG_ALIPHATIC_H,
        "h2": TAG_ALIPHATIC_H,
        "h3": TAG_ALIPHATIC_H,
        "h4": TAG_AROMATIC_H,    # H next to ring N
        "h5": TAG_AROMATIC_H,
        # Halogens
        "f":  TAG_HALOGEN_F,
        "cl": TAG_HALOGEN_CL,
        "br": TAG_HALOGEN_BR,
        "i":  TAG_HALOGEN_I,
        # Sulfur / Phosphorus
        "s":  TAG_SULFIDE_S,
        "sh": TAG_SULFIDE_S,
        "ss": TAG_SULFIDE_S,
        "p":  TAG_PHOSPHATE_P,
        "p3": TAG_PHOSPHATE_P,
        "p4": TAG_PHOSPHATE_P,
        "p5": TAG_PHOSPHATE_P,
    },
)


OPLS_AA = FunctionalTagMapping(
    force_field="OPLS-AA",
    description="Jorgensen OPLS-AA (1996), GROMACS oplsaa.ff style names",
    type_to_tag={
        # GROMACS oplsaa.ff naming convention
        "opls_267": TAG_CARBONYL_C,   # carboxylic acid C
        "opls_268": TAG_HYDROXYL_O,   # acid OH
        "opls_269": TAG_CARBONYL_O,   # acid C=O
        "opls_270": TAG_HYDROXYL_H,   # acid H
        "opls_154": TAG_HYDROXYL_O,   # alcohol O
        "opls_155": TAG_HYDROXYL_H,   # alcohol H
        "opls_157": TAG_SP3_C,        # C bonded to alcohol O
        "opls_235": TAG_CARBONYL_C,   # amide C
        "opls_236": TAG_CARBONYL_O,   # amide O
        "opls_237": TAG_AMIDE_N,      # amide N (primary)
        "opls_238": TAG_AMIDE_H,      # amide H
        "opls_239": TAG_AMIDE_N,      # amide N (secondary)
        "opls_240": TAG_AMIDE_H,      # H on sec amide N
        "opls_241": TAG_AMIDE_N,      # amide N (tertiary)
        "opls_135": TAG_SP3_C,
        "opls_140": TAG_ALIPHATIC_H,
        "opls_145": TAG_AROMATIC_C,
        "opls_146": TAG_AROMATIC_H,
        # Common alternative names without the opls_ prefix
        "CT": TAG_SP3_C,
        "CA": TAG_AROMATIC_C,
        "HA": TAG_AROMATIC_H,
        "HC": TAG_ALIPHATIC_H,
        "OH": TAG_HYDROXYL_O,
        "HO": TAG_HYDROXYL_H,
        "OA": TAG_HYDROXYL_O,    # alcohol O
        "O":  TAG_CARBONYL_O,
        "N":  TAG_AMIDE_N,
        "H":  TAG_AMIDE_H,
    },
)


CHARMM36 = FunctionalTagMapping(
    force_field="CHARMM36",
    description="MacKerell CHARMM36 (Klauda variant via GROMACS port)",
    type_to_tag={
        # Carbons
        "CT1": TAG_SP3_C,         # sp3 with 1 H (CH)
        "CT2": TAG_SP3_C,         # CH2
        "CT3": TAG_SP3_C,         # CH3
        "CC": TAG_CARBONYL_C,     # carboxyl C
        "C":  TAG_CARBONYL_C,     # amide / carbonyl C
        "CD": TAG_CARBONYL_C,     # carboxyl C (alt)
        "CA": TAG_AROMATIC_C,     # aromatic C
        "CPH1": TAG_AROMATIC_C,   # histidine ring
        "CPH2": TAG_AROMATIC_C,
        # Oxygens
        "OC":  TAG_CARBONYL_O,    # carboxyl C=O (charged form)
        "OB":  TAG_CARBONYL_O,    # carboxyl C=O (neutral)
        "O":   TAG_CARBONYL_O,    # amide / carbonyl O
        "OH1": TAG_HYDROXYL_O,    # alcohol/carboxyl OH
        "OS":  TAG_ETHER_O,       # ester O
        # Nitrogens
        "NH1": TAG_AMIDE_N,       # peptide N (secondary amide)
        "NH2": TAG_AMIDE_N,       # primary amide N
        "NH3": TAG_AMINE_N,       # cationic amine (Lys etc.)
        "NP":  TAG_AMIDE_N,       # proline-like (tertiary amide)
        "NC2": TAG_AMIDE_N,       # arginine guanidinium-like
        # Hydrogens
        "H":   TAG_AMIDE_H,       # polar H on N (amide donor)
        "HC":  TAG_AMIDE_H,       # polar H (charged)
        "HA1": TAG_ALIPHATIC_H,   # sp3 C-H (1 H type)
        "HA2": TAG_ALIPHATIC_H,
        "HA3": TAG_ALIPHATIC_H,
        "HP":  TAG_AROMATIC_H,    # aromatic ring H
        "HR1": TAG_AROMATIC_H,
        "HR2": TAG_AROMATIC_H,
        "HR3": TAG_AROMATIC_H,
        "HO": TAG_HYDROXYL_H,     # alcohol/carboxyl OH
        "HS": TAG_SULFIDE_S,      # actually wrong, just placeholder
        # Halogens
        "FA": TAG_HALOGEN_F,
        "CLA": TAG_HALOGEN_CL,
        "BRA": TAG_HALOGEN_BR,
        "IA": TAG_HALOGEN_I,
        # Sulfur
        "S": TAG_SULFIDE_S,
        "SM": TAG_SULFIDE_S,
    },
)


OPENFF_SAGE = FunctionalTagMapping(
    force_field="OpenFF-Sage-2",
    description="OpenFF Sage 2.x (SMIRNOFF) -- uses GAFF-like names after parmed conversion",
    # OpenFF -> GROMACS via openff-interchange typically emits names like
    # `c`, `oh`, `o`, etc. (same as GAFF), so reuse GAFF2 mapping below
    type_to_tag=dict(GAFF2.type_to_tag),
)


# Registry of built-in FFs
BUILTIN_MAPPINGS: Dict[str, FunctionalTagMapping] = {
    "GAFF2": GAFF2,
    "GAFF":  GAFF2,           # GAFF1 has same atom-type names for our purposes
    "OPLS-AA": OPLS_AA,
    "OPLS":  OPLS_AA,
    "CHARMM36": CHARMM36,
    "CHARMM": CHARMM36,
    "OpenFF": OPENFF_SAGE,
    "OpenFF-Sage": OPENFF_SAGE,
    "OpenFF-Sage-2": OPENFF_SAGE,
}


def get_mapping(force_field: str) -> FunctionalTagMapping:
    """Look up a built-in mapping by force-field name (case-insensitive)."""
    for key, mapping in BUILTIN_MAPPINGS.items():
        if key.lower() == force_field.lower():
            return mapping
    raise KeyError(
        f"Unknown force field {force_field!r}. "
        f"Built-in: {sorted(set(BUILTIN_MAPPINGS))}"
    )


def detect_force_field(
    atom_types: Iterable[str],
    candidates: Optional[Iterable[FunctionalTagMapping]] = None,
) -> str:
    """Heuristically detect the FF from a sample of atom types.

    Picks the FF whose `type_to_tag` keys overlap the most with the input
    atom types. Ties broken by first-listed FF in `candidates` (default
    GAFF2 > OPLS-AA > CHARMM36).

    Parameters
    ----------
    atom_types : iterable of str
        Atom types observed in the system (e.g. from a BDF topology).
    candidates : iterable of FunctionalTagMapping, optional
        FFs to consider. Default: [GAFF2, OPLS_AA, CHARMM36].

    Returns
    -------
    str : detected FF name (the ``force_field`` attribute)
    """
    if candidates is None:
        candidates = [GAFF2, OPLS_AA, CHARMM36]
    types = set(atom_types)
    best_match = -1
    best_ff = None
    for mapping in candidates:
        overlap = len(types & mapping.known_types())
        if overlap > best_match:
            best_match = overlap
            best_ff = mapping.force_field
    if best_ff is None or best_match == 0:
        raise RuntimeError(
            f"Could not detect force field from atom types: {sorted(types)[:10]}..."
        )
    return best_ff


def tag_atoms(
    atom_types: List[str], mapping: FunctionalTagMapping,
) -> List[Optional[str]]:
    """Apply a mapping to a list of atom types.

    Returns a list of tags (same length as input); ``None`` for unmapped types.
    """
    return [mapping.get_tag(t) for t in atom_types]


def fallback_tag_by_element(
    atom_names: List[str],
    bonds: List[tuple],
    initial_tags: List[Optional[str]],
) -> List[Optional[str]]:
    """Element + bond-graph fallback for atoms whose ``atom_type`` is unknown.

    Used when the force field (e.g. OpenFF SMIRNOFF) writes per-atom unique
    names like ``MOL0_0`` into the UDF ``Atom_Type_Name`` field, leaving the
    canonical-type lookup with no hit. Tags only the atoms whose entry in
    ``initial_tags`` is ``None``; entries already tagged by the FF mapping
    are left untouched (mapping wins over fallback).

    Parameters
    ----------
    atom_names : ``Atom_Name`` strings (element symbol prefix, e.g. ``"O"``,
                 ``"C1"``, ``"HO"``).
    bonds : list of ``(atom1_local_index, atom2_local_index)`` pairs.
    initial_tags : per-atom tag list; ``None`` for atoms needing fallback.

    Rules (element symbol = first uppercase letter of ``Atom_Name``):

    - **O**: bonded to ≥1 H → ``hydroxyl_O`` (alcohol or carboxyl OH side).
             No H → ``carbonyl_O``.
    - **H**: bonded to O → ``hydroxyl_H``; bonded to N → ``amide_H``.
    - **N**: bonded to C → ``amide_N`` (amide / amine distinction is done by
             the higher-level ``detect_amides`` / ``detect_amine_donors``).
    - **C**: bonded to a ``carbonyl_O`` neighbor (= O without H, tagged in
             pass 1) → ``carbonyl_C``.
    """
    from collections import defaultdict

    adj: Dict[int, List[int]] = defaultdict(list)
    for (a, b) in bonds:
        adj[a].append(b)
        adj[b].append(a)

    def _elem(name: Optional[str]) -> str:
        s = (name or "").strip()
        return s[0].upper() if s else ""

    out_tags = list(initial_tags)
    # Pass 1: tag O / H / N by their immediate-neighbour element symbols.
    for i, nm in enumerate(atom_names):
        if out_tags[i] is not None:
            continue
        e = _elem(nm)
        ne = [_elem(atom_names[j]) for j in adj[i]]
        if e == "O":
            if "H" in ne:
                out_tags[i] = TAG_HYDROXYL_O
            elif ne:  # bonded but no H → carbonyl-like
                out_tags[i] = TAG_CARBONYL_O
        elif e == "H":
            if "O" in ne:
                out_tags[i] = TAG_HYDROXYL_H
            elif "N" in ne:
                out_tags[i] = TAG_AMIDE_H
        elif e == "N":
            if "C" in ne:
                out_tags[i] = TAG_AMIDE_N

    # Pass 2: tag C using the O tags assigned in pass 1.
    # A C atom bonded to a TAG_CARBONYL_O (O with no H) is a carbonyl_C.
    for i, nm in enumerate(atom_names):
        if out_tags[i] is not None:
            continue
        e = _elem(nm)
        if e == "C":
            for j in adj[i]:
                if out_tags[j] == TAG_CARBONYL_O:
                    out_tags[i] = TAG_CARBONYL_C
                    break
    return out_tags
