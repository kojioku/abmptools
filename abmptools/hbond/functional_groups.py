"""
functional_groups.py
--------------------
Auto-detect functional groups from GAFF2 atom types + bond graph.

Detects:
- carboxyl  COOH:  c (sp2) - oh - ho  AND  c - o (C=O)
- amide     C(=O)-N: c - n  AND  c - o  (tertiary if n has no H)
- hydroxyl  O-H:   oh - ho

Atom indices in returned dataclasses are 0-based LOCAL atom indices within
each molecule (i.e. the position in molecule.atoms list). This matches the
convention in Set_of_Molecules.molecule[].bond[k].atom1/atom2 (which uses
0-based local indices, NOT global Atom_ID).

Use these local indices with TrajectoryFrame.positions[mol_index][local_idx]
to retrieve 3D coordinates.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .bdf_reader import MoleculeTopology


@dataclass
class CarboxylGroup:
    """COOH group: c=carbonyl C, o=C=O, oh=hydroxyl O, ho=hydroxyl H."""
    mol_index: int
    c_atom: int   # local index in molecule.atoms
    o_atom: int
    oh_atom: int
    ho_atom: int


@dataclass
class AmideGroup:
    """C(=O)-N: c=carbonyl C, o=C=O, n=amide N. tert=True means N has no H."""
    mol_index: int
    c_atom: int
    o_atom: int
    n_atom: int
    tert: bool


@dataclass
class HydroxylGroup:
    """Free O-H (not part of carboxyl): oh-ho."""
    mol_index: int
    oh_atom: int
    ho_atom: int


def _build_adjacency(topo: MoleculeTopology) -> Dict[int, List[Tuple[int, str]]]:
    """Build adj[local_atom_idx] = [(neighbor_local_idx, neighbor_type), ...].

    Bond table atom1/atom2 are 0-based local atom indices within the molecule.
    """
    adj: Dict[int, List[Tuple[int, str]]] = defaultdict(list)
    # local index = position in topo.atoms
    type_of = {idx: a.atom_type for idx, a in enumerate(topo.atoms)}
    for bond in topo.bonds:
        a, b = bond.atom1, bond.atom2
        adj[a].append((b, type_of.get(b, "?")))
        adj[b].append((a, type_of.get(a, "?")))
    return adj


def detect_carboxyls(
    molecules: List[MoleculeTopology],
) -> List[CarboxylGroup]:
    """Find all COOH groups in all molecules."""
    out = []
    for mi, topo in enumerate(molecules):
        adj = _build_adjacency(topo)
        type_of = {idx: a.atom_type for idx, a in enumerate(topo.atoms)}
        for c_atom, c_type in type_of.items():
            if c_type != "c":
                continue
            neighbors = adj[c_atom]
            o_neighbors = [n for n, t in neighbors if t == "o"]
            oh_neighbors = [n for n, t in neighbors if t == "oh"]
            if len(o_neighbors) == 1 and len(oh_neighbors) == 1:
                oh = oh_neighbors[0]
                ho_neighbors = [n for n, t in adj[oh] if t == "ho"]
                if len(ho_neighbors) == 1:
                    out.append(CarboxylGroup(
                        mol_index=mi,
                        c_atom=c_atom,
                        o_atom=o_neighbors[0],
                        oh_atom=oh,
                        ho_atom=ho_neighbors[0],
                    ))
    return out


def detect_amides(
    molecules: List[MoleculeTopology],
) -> List[AmideGroup]:
    """Find all C(=O)-N groups. tert=True iff N has no H neighbor."""
    out = []
    for mi, topo in enumerate(molecules):
        adj = _build_adjacency(topo)
        type_of = {idx: a.atom_type for idx, a in enumerate(topo.atoms)}
        for c_atom, c_type in type_of.items():
            if c_type != "c":
                continue
            neighbors = adj[c_atom]
            o_neighbors = [n for n, t in neighbors if t == "o"]
            n_neighbors = [n for n, t in neighbors if t == "n"]
            if len(o_neighbors) == 1 and len(n_neighbors) == 1:
                n_atom = n_neighbors[0]
                h_on_n = [n for n, t in adj[n_atom] if t in ("hn", "h", "ho")]
                tert = (len(h_on_n) == 0)
                out.append(AmideGroup(
                    mol_index=mi,
                    c_atom=c_atom,
                    o_atom=o_neighbors[0],
                    n_atom=n_atom,
                    tert=tert,
                ))
    return out


def detect_hydroxyls(
    molecules: List[MoleculeTopology],
    exclude_carboxyl: bool = True,
) -> List[HydroxylGroup]:
    """Find all O-H groups. Excludes carboxyl OH if exclude_carboxyl=True."""
    out = []
    carboxyl_oh = set()
    if exclude_carboxyl:
        for cg in detect_carboxyls(molecules):
            carboxyl_oh.add((cg.mol_index, cg.oh_atom))
    for mi, topo in enumerate(molecules):
        adj = _build_adjacency(topo)
        type_of = {idx: a.atom_type for idx, a in enumerate(topo.atoms)}
        for oh_atom, oh_type in type_of.items():
            if oh_type != "oh":
                continue
            if (mi, oh_atom) in carboxyl_oh:
                continue
            ho_neighbors = [n for n, t in adj[oh_atom] if t == "ho"]
            if len(ho_neighbors) == 1:
                out.append(HydroxylGroup(
                    mol_index=mi,
                    oh_atom=oh_atom,
                    ho_atom=ho_neighbors[0],
                ))
    return out


def summarize_groups(molecules: List[MoleculeTopology]) -> dict:
    """High-level summary of functional groups found."""
    carboxyls = detect_carboxyls(molecules)
    amides = detect_amides(molecules)
    hydroxyls = detect_hydroxyls(molecules)
    return {
        "carboxyl": carboxyls,
        "amide": amides,
        "hydroxyl": hydroxyls,
        "n_mols_with_carboxyl": len(set(g.mol_index for g in carboxyls)),
        "n_mols_with_amide": len(set(g.mol_index for g in amides)),
        "n_carboxyls": len(carboxyls),
        "n_amides": len(amides),
    }
