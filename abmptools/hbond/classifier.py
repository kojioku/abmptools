"""
classifier.py
-------------
Classify per-molecule role from a list of detected H-bonds.

Two classifications:
- dual: mol pair (i, j) where i's carboxyl OH→j's carboxyl O=C AND
        j's carboxyl OH→i's carboxyl O=C are BOTH present
- single: i's carboxyl OH→j's amide C=O (any direction)
- free: neither

Priority: a molecule participating in any dual gets "dual" role.
Otherwise, any single → "single". Else "free".
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple

from .functional_groups import AmideGroup, CarboxylGroup
from .hbond_detector import HBond


@dataclass
class MolRole:
    """Role assignment for one molecule."""
    mol_index: int
    role: str                       # "dual", "single", "free"
    partners: Set[int]              # mol indices it's bonded to (under chosen role)


@dataclass
class ClassificationResult:
    """Top-level classification result."""
    n_molecules: int
    roles: List[MolRole]
    dual_pairs: List[Tuple[int, int]]      # (i, j) with i < j
    single_pairs: List[Tuple[int, int]]    # (donor_mol, acceptor_mol)
    n_dual_mols: int
    n_single_mols: int
    n_free_mols: int


def classify(
    n_molecules: int,
    hbonds_carb_carb: List[HBond],
    hbonds_carb_amide: List[HBond],
    carboxyls: List[CarboxylGroup] = None,
    amides: List[AmideGroup] = None,
) -> ClassificationResult:
    """Classify each molecule based on detected H-bond network.

    Parameters
    ----------
    n_molecules : total molecule count
    hbonds_carb_carb : H-bonds where donor=carboxyl OH, acceptor=carboxyl O=C
    hbonds_carb_amide : H-bonds where donor=carboxyl OH, acceptor=amide O=C
    carboxyls, amides : functional groups (currently unused but kept for future
                        extension; e.g. multiple carboxyls per mol)
    """
    pair_donors: Dict[Tuple[int, int], List[HBond]] = defaultdict(list)
    for hb in hbonds_carb_carb:
        key = (hb.donor_mol, hb.acceptor_mol)
        pair_donors[key].append(hb)

    dual_pairs: Set[Tuple[int, int]] = set()
    for (i, j), _ in pair_donors.items():
        if (j, i) in pair_donors:
            pair = (min(i, j), max(i, j))
            dual_pairs.add(pair)

    single_pairs: List[Tuple[int, int]] = [
        (hb.donor_mol, hb.acceptor_mol) for hb in hbonds_carb_amide
    ]

    dual_mol_partners: Dict[int, Set[int]] = defaultdict(set)
    for i, j in dual_pairs:
        dual_mol_partners[i].add(j)
        dual_mol_partners[j].add(i)

    single_mol_partners: Dict[int, Set[int]] = defaultdict(set)
    for don, acc in single_pairs:
        single_mol_partners[don].add(acc)
        single_mol_partners[acc].add(don)

    roles: List[MolRole] = []
    n_dual = 0
    n_single = 0
    n_free = 0
    for m in range(n_molecules):
        if m in dual_mol_partners:
            roles.append(MolRole(
                mol_index=m, role="dual", partners=dual_mol_partners[m]
            ))
            n_dual += 1
        elif m in single_mol_partners:
            roles.append(MolRole(
                mol_index=m, role="single", partners=single_mol_partners[m]
            ))
            n_single += 1
        else:
            roles.append(MolRole(
                mol_index=m, role="free", partners=set()
            ))
            n_free += 1

    return ClassificationResult(
        n_molecules=n_molecules,
        roles=roles,
        dual_pairs=sorted(dual_pairs),
        single_pairs=single_pairs,
        n_dual_mols=n_dual,
        n_single_mols=n_single,
        n_free_mols=n_free,
    )
