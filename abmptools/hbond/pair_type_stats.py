"""
pair_type_stats.py
------------------
Generic donor-type × acceptor-type H-bond statistics, an alternative to the
COOH-specific dual / chain / single / free classifier (``classifier.py``) for
non-IMC-like systems (PVA, peptides, alcohols, polymer blends, etc.).

Use this when:
- you set ``AnalyzerConfig.classify_mode = "generic"``
- the donor/acceptor selection does not center on COOH carbonyl species
- you want the report aggregated by (donor functional-group type,
  acceptor functional-group type) instead of the dual/chain/single/free axis.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from .hbond_detector import HBond


@dataclass
class PairTypeStat:
    """Aggregate of all H-bonds between one donor type and one acceptor type."""
    donor_type: str
    acceptor_type: str
    n_hbonds: int = 0
    unique_donors: Set[Tuple[int, int]] = field(default_factory=set)
    unique_acceptors: Set[Tuple[int, int]] = field(default_factory=set)

    @property
    def n_uniq_donors(self) -> int:
        return len(self.unique_donors)

    @property
    def n_uniq_acceptors(self) -> int:
        return len(self.unique_acceptors)


@dataclass
class GenericPairClassification:
    """Top-level result of ``classify_generic``.

    Attributes
    ----------
    pair_stats : per (donor_type, acceptor_type) → PairTypeStat
    atom_role : (mol_index, atom_local_idx) → role string in
                {"Donor", "Acceptor", "Both", "Candidate"}.
                ``Candidate`` = listed as a donor or acceptor candidate but
                not engaged in any H-bond in this frame.
    n_donor_candidates / n_acceptor_candidates : group-type totals (for the
                ratio_*_busy denominator).
    """
    pair_stats: Dict[Tuple[str, str], PairTypeStat]
    atom_role: Dict[Tuple[int, int], str]
    n_donor_candidates: Dict[str, int]
    n_acceptor_candidates: Dict[str, int]


def classify_generic(
    donor_groups_by_type: Dict[str, List[Tuple[int, int, int]]],
    acceptor_groups_by_type: Dict[str, List[Tuple[int, int]]],
    hbonds_by_pair_type: Dict[Tuple[str, str], List[HBond]],
) -> GenericPairClassification:
    """Aggregate generic donor × acceptor H-bond statistics.

    Parameters
    ----------
    donor_groups_by_type : ``{donor_type: [(mol_idx, donor_d, donor_h), ...]}``
        e.g. ``{'hydroxyl': [(0, 3, 40), (1, 3, 40), ...]}``
    acceptor_groups_by_type : ``{acceptor_type: [(mol_idx, acceptor_a), ...]}``
        e.g. ``{'hydroxyl_O': [(0, 3), ...], 'amide_O': [(0, 2), ...]}``
    hbonds_by_pair_type : ``{(donor_type, acceptor_type): [HBond, ...]}``
        Pre-detected H-bonds bucketed by the participating types. The donor
        atom of each HBond must match an entry in
        ``donor_groups_by_type[donor_type]``, and likewise for acceptors.
    """
    pair_stats: Dict[Tuple[str, str], PairTypeStat] = {}
    atom_role: Dict[Tuple[int, int], str] = {}

    # Initialize every candidate atom as "Candidate"
    for dt, sites in donor_groups_by_type.items():
        for (mi, d, h) in sites:
            atom_role.setdefault((mi, d), "Candidate")
            atom_role.setdefault((mi, h), "Candidate")
    for at, sites in acceptor_groups_by_type.items():
        for (mi, a) in sites:
            atom_role.setdefault((mi, a), "Candidate")

    def _bump_role(key: Tuple[int, int], new: str) -> None:
        cur = atom_role.get(key, "Candidate")
        if cur == "Both":
            return
        if cur in ("Candidate",) or cur == new:
            atom_role[key] = new
            return
        # cur and new differ and neither is Both → Both
        atom_role[key] = "Both"

    for (dt, at), hblist in hbonds_by_pair_type.items():
        st = PairTypeStat(donor_type=dt, acceptor_type=at)
        for hb in hblist:
            st.n_hbonds += 1
            st.unique_donors.add((hb.donor_mol, hb.donor_d))
            st.unique_acceptors.add((hb.acceptor_mol, hb.acceptor_a))
            _bump_role((hb.donor_mol, hb.donor_d), "Donor")
            _bump_role((hb.donor_mol, hb.donor_h), "Donor")
            _bump_role((hb.acceptor_mol, hb.acceptor_a), "Acceptor")
        pair_stats[(dt, at)] = st

    n_donor_candidates = {dt: len(s) for dt, s in donor_groups_by_type.items()}
    n_acceptor_candidates = {at: len(s) for at, s in acceptor_groups_by_type.items()}

    return GenericPairClassification(
        pair_stats=pair_stats,
        atom_role=atom_role,
        n_donor_candidates=n_donor_candidates,
        n_acceptor_candidates=n_acceptor_candidates,
    )


def summarize_pair_stats(
    result: GenericPairClassification,
) -> Dict[str, float]:
    """Quick per-type aggregate suitable for stdout log lines.

    Returns a flat dict like::

        {
          "n_hbonds_total": 132,
          "hydroxyl->hydroxyl_O": 25,
          "hydroxyl->amide_O": 7,
          "ratio_donor_busy_hydroxyl": 0.30,
          "ratio_acceptor_busy_hydroxyl_O": 0.34,
        }
    """
    out: Dict[str, float] = {"n_hbonds_total": 0}
    donor_busy: Dict[str, int] = defaultdict(int)
    acceptor_busy: Dict[str, int] = defaultdict(int)

    for (dt, at), st in result.pair_stats.items():
        out[f"{dt}->{at}"] = st.n_hbonds
        out["n_hbonds_total"] += st.n_hbonds
        donor_busy[dt] = max(donor_busy[dt], st.n_uniq_donors)
        acceptor_busy[at] = max(acceptor_busy[at], st.n_uniq_acceptors)

    for dt, n in donor_busy.items():
        denom = result.n_donor_candidates.get(dt, 0)
        out[f"ratio_donor_busy_{dt}"] = n / denom if denom else 0.0
    for at, n in acceptor_busy.items():
        denom = result.n_acceptor_candidates.get(at, 0)
        out[f"ratio_acceptor_busy_{at}"] = n / denom if denom else 0.0
    return out
