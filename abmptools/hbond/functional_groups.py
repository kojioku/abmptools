"""
functional_groups.py
--------------------
Force-field-agnostic functional group detection using canonical tags.

Detects:
- carboxyl  COOH:  carbonyl_C bonded to hydroxyl_O (+hydroxyl_H) AND carbonyl_O
- amide     C(=O)-N: carbonyl_C bonded to amide_N AND carbonyl_O
              tert=True iff N has no hydroxyl_H/amide_H neighbor (= tertiary amide)
- hydroxyl  O-H:   hydroxyl_O bonded to hydroxyl_H (excluding carboxyl OH by default)
- amine     N-H:   amine_N/amide_N bonded to amide_H (donor for secondary amides)

Atom indices in returned dataclasses are 0-based LOCAL atom indices within
each molecule (i.e. the position in molecule.atoms list). This matches the
convention in Set_of_Molecules.molecule[].bond[k].atom1/atom2 (which uses
0-based local indices, NOT global Atom_ID).

The previous GAFF2-specific implementation has been refactored to dispatch
through ``func_tags`` so OPLS-AA / CHARMM36 / OpenFF systems can use the
same detection logic.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .bdf_reader import MoleculeTopology
from .func_tags import (
    FunctionalTagMapping, GAFF2,
    TAG_AMIDE_H, TAG_AMIDE_N, TAG_AMINE_N,
    TAG_CARBONYL_C, TAG_CARBONYL_O,
    TAG_HYDROXYL_H, TAG_HYDROXYL_O,
    detect_force_field, fallback_tag_by_element, get_mapping,
)


# Module-level flag (mirrors ``AnalyzerConfig.use_element_fallback``); the
# analyzer flips it per run, but it can also be toggled directly for callers
# that bypass the analyzer.
USE_ELEMENT_FALLBACK = True


@dataclass
class CarboxylGroup:
    """COOH group: c=carbonyl C, o=C=O, oh=hydroxyl O, ho=hydroxyl H."""
    mol_index: int
    c_atom: int     # local index in molecule.atoms
    o_atom: int
    oh_atom: int
    ho_atom: int


@dataclass
class AmideGroup:
    """C(=O)-N: c=carbonyl C, o=C=O, n=amide N.

    ``tert=True`` means N has no donor H attached (= tertiary amide).
    ``nh_atom`` is the local index of the N-H hydrogen if present, else None.
    """
    mol_index: int
    c_atom: int
    o_atom: int
    n_atom: int
    tert: bool
    nh_atom: Optional[int] = None


@dataclass
class HydroxylGroup:
    """Free O-H (typically alcohol, optionally excluding carboxyl OH)."""
    mol_index: int
    oh_atom: int
    ho_atom: int


@dataclass
class EtherGroup:
    """Ether-type oxygen acceptor: ``o_atom`` = O bonded to two carbons
    (C-O-C) with no hydrogen. Covers dialkyl / aryl ethers and the
    single-bonded ester oxygen."""
    mol_index: int
    o_atom: int


@dataclass
class AmineDonorGroup:
    """N-H donor site (typically secondary/primary amide, primary/secondary amine).

    Each H on the N is represented by a separate AmineDonorGroup
    (i.e. primary amide -NH2 yields TWO AmineDonorGroup entries).
    """
    mol_index: int
    n_atom: int
    h_atom: int
    from_amide: bool   # True if N is amide_N; False if amine_N


# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

def _tag_atoms_of_mol(
    topo: MoleculeTopology, mapping: FunctionalTagMapping,
) -> List[Optional[str]]:
    """Return per-atom tag list (local-index ordered).

    Falls back to element + bond-graph heuristics for atoms whose
    ``atom_type`` is not in the mapping (e.g. OpenFF SMIRNOFF's per-atom
    unique ``MOL0_0`` names). The fallback can be disabled by setting
    ``USE_ELEMENT_FALLBACK = False`` at module scope (typically set by
    ``AnalyzerConfig.use_element_fallback`` for a single run).
    """
    tags = [mapping.get_tag(a.atom_type) for a in topo.atoms]
    if USE_ELEMENT_FALLBACK and any(t is None for t in tags):
        atom_names = [a.atom_name for a in topo.atoms]
        bonds = [(b.atom1, b.atom2) for b in topo.bonds]
        tags = fallback_tag_by_element(atom_names, bonds, tags)
    return tags


def _build_adjacency(
    topo: MoleculeTopology, tags: List[Optional[str]],
) -> Dict[int, List[Tuple[int, Optional[str]]]]:
    """Build adj[local_atom_idx] = [(neighbor_local_idx, neighbor_tag), ...].

    Bond table atom1/atom2 are 0-based local atom indices within the molecule.
    """
    adj: Dict[int, List[Tuple[int, Optional[str]]]] = defaultdict(list)
    for bond in topo.bonds:
        a, b = bond.atom1, bond.atom2
        adj[a].append((b, tags[b] if 0 <= b < len(tags) else None))
        adj[b].append((a, tags[a] if 0 <= a < len(tags) else None))
    return adj


def _resolve_mapping(
    molecules: List[MoleculeTopology],
    mapping: Optional[FunctionalTagMapping] = None,
) -> FunctionalTagMapping:
    """If mapping is None, auto-detect from the first molecule's atom types."""
    if mapping is not None:
        return mapping
    sample_types = {a.atom_type for m in molecules for a in m.atoms}
    try:
        ff_name = detect_force_field(sample_types)
        return get_mapping(ff_name)
    except (KeyError, RuntimeError):
        return GAFF2  # fall back to GAFF2


# ----------------------------------------------------------------------------
# Public detection functions
# ----------------------------------------------------------------------------

def detect_carboxyls(
    molecules: List[MoleculeTopology],
    mapping: Optional[FunctionalTagMapping] = None,
) -> List[CarboxylGroup]:
    """Find all COOH groups across all molecules.

    Pattern: ``carbonyl_C`` bonded to **both** ``hydroxyl_O`` (with attached
    ``hydroxyl_H``) and ``carbonyl_O``.
    """
    mapping = _resolve_mapping(molecules, mapping)
    out: List[CarboxylGroup] = []
    for mi, topo in enumerate(molecules):
        tags = _tag_atoms_of_mol(topo, mapping)
        adj = _build_adjacency(topo, tags)
        for c_atom, c_tag in enumerate(tags):
            if c_tag != TAG_CARBONYL_C:
                continue
            o_neigh = [n for n, t in adj[c_atom] if t == TAG_CARBONYL_O]
            oh_neigh = [n for n, t in adj[c_atom] if t == TAG_HYDROXYL_O]
            if len(o_neigh) == 1 and len(oh_neigh) == 1:
                oh = oh_neigh[0]
                ho_neigh = [n for n, t in adj[oh] if t == TAG_HYDROXYL_H]
                if len(ho_neigh) == 1:
                    out.append(CarboxylGroup(
                        mol_index=mi,
                        c_atom=c_atom,
                        o_atom=o_neigh[0],
                        oh_atom=oh,
                        ho_atom=ho_neigh[0],
                    ))
    return out


def detect_amides(
    molecules: List[MoleculeTopology],
    mapping: Optional[FunctionalTagMapping] = None,
) -> List[AmideGroup]:
    """Find all C(=O)-N (amide) groups.

    ``tert=True`` iff N has no ``amide_H`` neighbor.
    ``nh_atom`` is the local index of one such H (or None for tertiary).
    """
    mapping = _resolve_mapping(molecules, mapping)
    out: List[AmideGroup] = []
    for mi, topo in enumerate(molecules):
        tags = _tag_atoms_of_mol(topo, mapping)
        adj = _build_adjacency(topo, tags)
        for c_atom, c_tag in enumerate(tags):
            if c_tag != TAG_CARBONYL_C:
                continue
            o_neigh = [n for n, t in adj[c_atom] if t == TAG_CARBONYL_O]
            n_neigh = [n for n, t in adj[c_atom]
                       if t in (TAG_AMIDE_N, TAG_AMINE_N)]
            if len(o_neigh) == 1 and len(n_neigh) == 1:
                n_atom = n_neigh[0]
                h_on_n = [n for n, t in adj[n_atom] if t == TAG_AMIDE_H]
                tert = (len(h_on_n) == 0)
                out.append(AmideGroup(
                    mol_index=mi,
                    c_atom=c_atom,
                    o_atom=o_neigh[0],
                    n_atom=n_atom,
                    tert=tert,
                    nh_atom=(h_on_n[0] if h_on_n else None),
                ))
    return out


def detect_hydroxyls(
    molecules: List[MoleculeTopology],
    exclude_carboxyl: bool = True,
    mapping: Optional[FunctionalTagMapping] = None,
) -> List[HydroxylGroup]:
    """Find O-H groups. Excludes carboxyl OH if ``exclude_carboxyl=True``."""
    mapping = _resolve_mapping(molecules, mapping)
    out: List[HydroxylGroup] = []
    carboxyl_oh = set()
    if exclude_carboxyl:
        for cg in detect_carboxyls(molecules, mapping):
            carboxyl_oh.add((cg.mol_index, cg.oh_atom))
    for mi, topo in enumerate(molecules):
        tags = _tag_atoms_of_mol(topo, mapping)
        adj = _build_adjacency(topo, tags)
        for oh_atom, oh_tag in enumerate(tags):
            if oh_tag != TAG_HYDROXYL_O:
                continue
            if (mi, oh_atom) in carboxyl_oh:
                continue
            ho_neigh = [n for n, t in adj[oh_atom] if t == TAG_HYDROXYL_H]
            if len(ho_neigh) == 1:
                out.append(HydroxylGroup(
                    mol_index=mi, oh_atom=oh_atom, ho_atom=ho_neigh[0]
                ))
    return out


def detect_ethers(
    molecules: List[MoleculeTopology],
    mapping: Optional[FunctionalTagMapping] = None,
) -> List[EtherGroup]:
    """Find ether-type oxygen acceptors: an O bonded to **exactly two carbons**
    (C-O-C) with **no hydrogen**.

    Covers dialkyl / aryl ethers and the single-bonded ester oxygen; excludes
    hydroxyl O (has an H neighbour) and carbonyl O (a single heavy neighbour,
    i.e. the C=O appears as degree-1 in the bond graph). Element + connectivity
    based (not tag based), so it works for force fields that have no dedicated
    ether atom type — e.g. OpenFF SMIRNOFF, whose per-atom-unique names make the
    element fallback tag an ether O as ``carbonyl_O``. ``mapping`` is accepted
    for signature parity with the other detectors but is not used.
    """
    from .diagram import element_from_name
    out: List[EtherGroup] = []
    for mi, topo in enumerate(molecules):
        elems = [element_from_name(a.atom_name) for a in topo.atoms]
        adj: Dict[int, List[int]] = defaultdict(list)
        for b in topo.bonds:
            adj[b.atom1].append(b.atom2)
            adj[b.atom2].append(b.atom1)
        for ai, e in enumerate(elems):
            if e != "O":
                continue
            nbr_elems = [elems[j] for j in adj[ai] if 0 <= j < len(elems)]
            if "H" in nbr_elems:
                continue   # hydroxyl / water oxygen, not an ether
            heavy = [x for x in nbr_elems if x != "H"]
            if len(heavy) == 2 and all(x == "C" for x in heavy):
                out.append(EtherGroup(mol_index=mi, o_atom=ai))
    return out


def detect_amine_donors(
    molecules: List[MoleculeTopology],
    include_amide: bool = True,
    include_amine: bool = True,
    mapping: Optional[FunctionalTagMapping] = None,
) -> List[AmineDonorGroup]:
    """Find N-H donor sites.

    Returns one entry per N-H bond (so primary amide -NH2 yields 2 entries).

    Parameters
    ----------
    include_amide : if True, include amide_N - amide_H pairs (secondary/primary amides)
    include_amine : if True, include amine_N - amide_H pairs (primary/secondary amines)
    """
    mapping = _resolve_mapping(molecules, mapping)
    out: List[AmineDonorGroup] = []
    for mi, topo in enumerate(molecules):
        tags = _tag_atoms_of_mol(topo, mapping)
        adj = _build_adjacency(topo, tags)
        for n_atom, n_tag in enumerate(tags):
            if n_tag == TAG_AMIDE_N:
                if not include_amide:
                    continue
                from_amide = True
            elif n_tag == TAG_AMINE_N:
                if not include_amine:
                    continue
                from_amide = False
            else:
                continue
            for h_neighbor, h_tag in adj[n_atom]:
                if h_tag == TAG_AMIDE_H:
                    out.append(AmineDonorGroup(
                        mol_index=mi,
                        n_atom=n_atom,
                        h_atom=h_neighbor,
                        from_amide=from_amide,
                    ))
    return out


def summarize_groups(
    molecules: List[MoleculeTopology],
    mapping: Optional[FunctionalTagMapping] = None,
) -> dict:
    """High-level summary of functional groups detected.

    Now includes secondary amide donor counts (N-H).
    """
    mapping = _resolve_mapping(molecules, mapping)
    carboxyls = detect_carboxyls(molecules, mapping)
    amides = detect_amides(molecules, mapping)
    hydroxyls = detect_hydroxyls(molecules, mapping=mapping)
    n_donors = detect_amine_donors(molecules, mapping=mapping)
    return {
        "force_field": mapping.force_field,
        "carboxyl": carboxyls,
        "amide": amides,
        "hydroxyl": hydroxyls,
        "amine_donor": n_donors,
        "n_carboxyls": len(carboxyls),
        "n_amides": len(amides),
        "n_amides_tert": sum(1 for a in amides if a.tert),
        "n_amides_nonsec": sum(1 for a in amides if not a.tert),
        "n_hydroxyls": len(hydroxyls),
        "n_amine_donors": len(n_donors),
        "n_mols_with_carboxyl": len(set(g.mol_index for g in carboxyls)),
        "n_mols_with_amide": len(set(g.mol_index for g in amides)),
    }
