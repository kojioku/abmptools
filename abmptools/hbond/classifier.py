"""
classifier.py
-------------
Per-functional-group H-bond role classification.

A single molecule can contribute to multiple roles through different functional
groups (e.g. one COOH joining a cyclic dimer while another donates to an amide).
So the primary unit of classification is the *functional group*, not the
molecule.

Per-COOH roles:
- dual:   participates in a cyclic dimer with another COOH
          (i.e. both i→j and j→i carboxyl-carboxyl H-bonds exist for this COOH)
- single: not in a dual; donates an O-H to one or more amide carbonyls
- free:   neither dual nor single donor

Per-amide roles:
- accept: receives at least one carboxyl OH→C=O H-bond
- free:   no incoming carboxyl H-bond

Per-molecule representative role (for coloring only):
- ``dual``   if the molecule has any dual carboxyl
- ``single`` else if the molecule has any single donor carboxyl
- ``free``   otherwise
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from .functional_groups import AmideGroup, CarboxylGroup
from .hbond_detector import HBond


@dataclass
class CarboxylRole:
    """Role assignment for one COOH group.

    ``carboxyl_index`` is the 0-based position of this COOH within its molecule
    (so a molecule with two COOH groups yields 0 and 1).
    ``dual_partners`` are ``(partner_mol_index, partner_carboxyl_index)`` tuples.
    ``single_acceptors`` are ``(acceptor_mol_index, acceptor_amide_index)``.
    """
    mol_index: int
    carboxyl_index: int
    role: str                                              # "dual" / "single" / "free"
    dual_partners: Set[Tuple[int, int]] = field(default_factory=set)
    single_acceptors: Set[Tuple[int, int]] = field(default_factory=set)


@dataclass
class AmideRole:
    """Role assignment for one amide C=O group."""
    mol_index: int
    amide_index: int
    role: str                                              # "accept" / "free"
    donor_carboxyls: Set[Tuple[int, int]] = field(default_factory=set)


@dataclass
class MolRole:
    """Per-molecule representative role (for coloring) + per-group counts."""
    mol_index: int
    role: str                                              # "dual" / "single" / "free"
    partners: Set[int] = field(default_factory=set)        # mol-level partners under role
    n_carboxyls_dual: int = 0
    n_carboxyls_single: int = 0
    n_carboxyls_free: int = 0
    n_amides_accept: int = 0
    n_amides_free: int = 0


@dataclass
class FunctionalGroupClassification:
    """Top-level classification result.

    The mol-level fields (``roles``, ``n_dual_mols`` ...) are kept for
    backward compatibility (they drive ``colorize_udf``); the new
    per-functional-group fields (``n_carboxyls_*``, ``ratio_*``,
    ``carboxyl_roles``, ``amide_roles``) are the primary outputs.
    """
    n_molecules: int
    roles: List[MolRole]
    dual_pairs: List[Tuple[int, int]]                      # mol-level (i, j) with i < j
    single_pairs: List[Tuple[int, int]]                    # mol-level (donor, acceptor)
    n_dual_mols: int
    n_single_mols: int
    n_free_mols: int
    # Per-functional-group counts (primary metric)
    n_carboxyls: int = 0
    n_carboxyls_dual: int = 0
    n_carboxyls_single: int = 0
    n_carboxyls_free: int = 0
    n_amides: int = 0
    n_amides_accept: int = 0
    n_amides_free: int = 0
    # Ratios over total group counts
    ratio_carboxyl_dual: float = 0.0
    ratio_carboxyl_single: float = 0.0
    ratio_carboxyl_free: float = 0.0
    ratio_amide_accept: float = 0.0
    ratio_amide_free: float = 0.0
    # Per-group role lists
    carboxyl_roles: List[CarboxylRole] = field(default_factory=list)
    amide_roles: List[AmideRole] = field(default_factory=list)


ClassificationResult = FunctionalGroupClassification


def _safe_ratio(num: int, den: int) -> float:
    return num / den if den > 0 else 0.0


def classify(
    n_molecules: int,
    hbonds_carb_carb: List[HBond],
    hbonds_carb_amide: List[HBond],
    carboxyls: Optional[List[CarboxylGroup]] = None,
    amides: Optional[List[AmideGroup]] = None,
) -> FunctionalGroupClassification:
    """Classify each functional group based on detected H-bond network.

    Parameters
    ----------
    n_molecules : total molecule count (drives ``roles`` length)
    hbonds_carb_carb : H-bonds donor=carboxyl OH, acceptor=carboxyl O=C
    hbonds_carb_amide : H-bonds donor=carboxyl OH, acceptor=amide O=C
    carboxyls, amides : functional groups. If omitted, only the legacy
                        mol-level fields are populated (per-group fields stay
                        empty / zero).
    """
    carboxyls = list(carboxyls or [])
    amides = list(amides or [])

    # Indexing helpers --------------------------------------------------------
    carb_idx_in_mol: Dict[int, int] = {}
    _carb_seen: Dict[int, int] = defaultdict(int)
    for i, cg in enumerate(carboxyls):
        carb_idx_in_mol[i] = _carb_seen[cg.mol_index]
        _carb_seen[cg.mol_index] += 1

    amide_idx_in_mol: Dict[int, int] = {}
    _amide_seen: Dict[int, int] = defaultdict(int)
    for i, ag in enumerate(amides):
        amide_idx_in_mol[i] = _amide_seen[ag.mol_index]
        _amide_seen[ag.mol_index] += 1

    carb_donor_lookup = {
        (cg.mol_index, cg.oh_atom, cg.ho_atom): i
        for i, cg in enumerate(carboxyls)
    }
    carb_acceptor_lookup = {
        (cg.mol_index, cg.o_atom): i for i, cg in enumerate(carboxyls)
    }
    amide_acceptor_lookup = {
        (ag.mol_index, ag.o_atom): i for i, ag in enumerate(amides)
    }

    # COOH-COOH relations -----------------------------------------------------
    cc_relations: Set[Tuple[int, int]] = set()
    for hb in hbonds_carb_carb:
        d = carb_donor_lookup.get((hb.donor_mol, hb.donor_d, hb.donor_h))
        a = carb_acceptor_lookup.get((hb.acceptor_mol, hb.acceptor_a))
        if d is not None and a is not None:
            cc_relations.add((d, a))

    dual_carb_pairs: Set[Tuple[int, int]] = set()
    for (i, j) in cc_relations:
        if (j, i) in cc_relations:
            dual_carb_pairs.add((min(i, j), max(i, j)))

    carb_dual_partners: Dict[int, Set[int]] = defaultdict(set)
    for (i, j) in dual_carb_pairs:
        carb_dual_partners[i].add(j)
        carb_dual_partners[j].add(i)

    # COOH→amide relations ----------------------------------------------------
    ca_relations: List[Tuple[int, int]] = []
    for hb in hbonds_carb_amide:
        d = carb_donor_lookup.get((hb.donor_mol, hb.donor_d, hb.donor_h))
        a = amide_acceptor_lookup.get((hb.acceptor_mol, hb.acceptor_a))
        if d is not None and a is not None:
            ca_relations.append((d, a))

    carb_to_amides: Dict[int, Set[int]] = defaultdict(set)
    amide_from_carbs: Dict[int, Set[int]] = defaultdict(set)
    for (d, a) in ca_relations:
        carb_to_amides[d].add(a)
        amide_from_carbs[a].add(d)

    # Per-COOH role -----------------------------------------------------------
    carboxyl_roles: List[CarboxylRole] = []
    for i, cg in enumerate(carboxyls):
        if i in carb_dual_partners:
            role = "dual"
        elif i in carb_to_amides:
            role = "single"
        else:
            role = "free"
        dual_partners = {
            (carboxyls[j].mol_index, carb_idx_in_mol[j])
            for j in carb_dual_partners.get(i, set())
        }
        single_acceptors = {
            (amides[k].mol_index, amide_idx_in_mol[k])
            for k in carb_to_amides.get(i, set())
        }
        carboxyl_roles.append(CarboxylRole(
            mol_index=cg.mol_index,
            carboxyl_index=carb_idx_in_mol[i],
            role=role,
            dual_partners=dual_partners,
            single_acceptors=single_acceptors,
        ))

    # Per-amide role ----------------------------------------------------------
    amide_roles: List[AmideRole] = []
    for k, ag in enumerate(amides):
        role = "accept" if k in amide_from_carbs else "free"
        donor_carbs = {
            (carboxyls[d].mol_index, carb_idx_in_mol[d])
            for d in amide_from_carbs.get(k, set())
        }
        amide_roles.append(AmideRole(
            mol_index=ag.mol_index,
            amide_index=amide_idx_in_mol[k],
            role=role,
            donor_carboxyls=donor_carbs,
        ))

    # Aggregate counts --------------------------------------------------------
    n_carboxyls = len(carboxyls)
    n_carb_dual = sum(1 for r in carboxyl_roles if r.role == "dual")
    n_carb_single = sum(1 for r in carboxyl_roles if r.role == "single")
    n_carb_free = sum(1 for r in carboxyl_roles if r.role == "free")
    n_amides_n = len(amides)
    n_amide_acc = sum(1 for r in amide_roles if r.role == "accept")
    n_amide_free = sum(1 for r in amide_roles if r.role == "free")

    # Per-mol aggregation (representative role) -------------------------------
    mol_carb_counts: Dict[int, Dict[str, int]] = defaultdict(
        lambda: {"dual": 0, "single": 0, "free": 0}
    )
    for r in carboxyl_roles:
        mol_carb_counts[r.mol_index][r.role] += 1
    mol_amide_counts: Dict[int, Dict[str, int]] = defaultdict(
        lambda: {"accept": 0, "free": 0}
    )
    for r in amide_roles:
        mol_amide_counts[r.mol_index][r.role] += 1

    dual_mol_partners: Dict[int, Set[int]] = defaultdict(set)
    dual_mol_pairs: Set[Tuple[int, int]] = set()
    for (i, j) in dual_carb_pairs:
        mi, mj = carboxyls[i].mol_index, carboxyls[j].mol_index
        if mi != mj:
            dual_mol_partners[mi].add(mj)
            dual_mol_partners[mj].add(mi)
            dual_mol_pairs.add((min(mi, mj), max(mi, mj)))

    single_mol_partners: Dict[int, Set[int]] = defaultdict(set)
    for r in carboxyl_roles:
        if r.role == "single":
            for (acc_mol, _) in r.single_acceptors:
                if acc_mol != r.mol_index:
                    single_mol_partners[r.mol_index].add(acc_mol)
                    single_mol_partners[acc_mol].add(r.mol_index)

    mol_roles: List[MolRole] = []
    n_dual = n_single = n_free = 0
    for m in range(n_molecules):
        ccnt = mol_carb_counts.get(m, {"dual": 0, "single": 0, "free": 0})
        acnt = mol_amide_counts.get(m, {"accept": 0, "free": 0})
        if ccnt["dual"] > 0:
            role, partners = "dual", dual_mol_partners.get(m, set())
            n_dual += 1
        elif ccnt["single"] > 0:
            role, partners = "single", single_mol_partners.get(m, set())
            n_single += 1
        else:
            role, partners = "free", set()
            n_free += 1
        mol_roles.append(MolRole(
            mol_index=m,
            role=role,
            partners=partners,
            n_carboxyls_dual=ccnt["dual"],
            n_carboxyls_single=ccnt["single"],
            n_carboxyls_free=ccnt["free"],
            n_amides_accept=acnt["accept"],
            n_amides_free=acnt["free"],
        ))

    single_pairs_mol = [
        (hb.donor_mol, hb.acceptor_mol) for hb in hbonds_carb_amide
    ]

    return FunctionalGroupClassification(
        n_molecules=n_molecules,
        roles=mol_roles,
        dual_pairs=sorted(dual_mol_pairs),
        single_pairs=single_pairs_mol,
        n_dual_mols=n_dual,
        n_single_mols=n_single,
        n_free_mols=n_free,
        n_carboxyls=n_carboxyls,
        n_carboxyls_dual=n_carb_dual,
        n_carboxyls_single=n_carb_single,
        n_carboxyls_free=n_carb_free,
        n_amides=n_amides_n,
        n_amides_accept=n_amide_acc,
        n_amides_free=n_amide_free,
        ratio_carboxyl_dual=_safe_ratio(n_carb_dual, n_carboxyls),
        ratio_carboxyl_single=_safe_ratio(n_carb_single, n_carboxyls),
        ratio_carboxyl_free=_safe_ratio(n_carb_free, n_carboxyls),
        ratio_amide_accept=_safe_ratio(n_amide_acc, n_amides_n),
        ratio_amide_free=_safe_ratio(n_amide_free, n_amides_n),
        carboxyl_roles=carboxyl_roles,
        amide_roles=amide_roles,
    )
