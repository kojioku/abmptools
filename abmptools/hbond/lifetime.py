"""
lifetime.py
-----------
H-bond lifetime statistics across multi-record trajectories.

Three quantities are computed:

1. **Continuous lifetime** — longest run of consecutive frames in which a
   given (donor_d, donor_h, acceptor_a, donor_mol, acceptor_mol) pair is
   present. Strict definition: a single missing frame ends the run.

2. **Intermittent lifetime** — runs allowed to bridge ``gap_tolerance``
   consecutive missing frames. Captures rapid librations that briefly
   break the H-bond geometry without "reorganizing" the network.

3. **Luzar-Chandler autocorrelation** — population autocorrelation
   ``C(t) = <h(0)h(t)> / <h(0)>`` where ``h(t) ∈ {0,1}`` indicates pair
   presence at frame t. The integral ``τ_HB = ∫ C(t) dt`` is the standard
   H-bond lifetime in the Luzar-Chandler theory of liquid water.

References:
- A. Luzar and D. Chandler, Nature 379, 55-57 (1996).
- D. C. Rapaport, Mol. Phys. 50, 1151-1162 (1983).
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .hbond_detector import HBond


# A pair is uniquely identified by its donor+acceptor atom IDs
PairKey = Tuple[int, int, int, int, int]
# = (donor_mol, donor_d_local, donor_h_local, acceptor_mol, acceptor_a_local)


def hbond_to_key(hb: HBond) -> PairKey:
    return (hb.donor_mol, hb.donor_d, hb.donor_h,
            hb.acceptor_mol, hb.acceptor_a)


@dataclass
class PairLifetime:
    """Lifetime statistics for one unique donor-acceptor pair."""
    pair: PairKey
    present_records: List[int]                # records where the pair was present
    continuous_intervals: List[Tuple[int, int]]  # (start, end) inclusive, strict
    continuous_max: int                       # longest strict continuous run
    continuous_mean: float                    # mean strict run length
    intermittent_intervals: List[Tuple[int, int]]  # bridged by gap_tolerance
    intermittent_max: int
    intermittent_mean: float
    total_present: int                        # total frames present
    occupancy: float                          # total_present / n_records


def _runs(present_set: set, n_records: int,
          gap_tolerance: int = 0) -> List[Tuple[int, int]]:
    """Compute (start, end) intervals from a presence set.

    With gap_tolerance=k, gaps of up to k missing frames inside a run are
    bridged (treated as part of the same interval).
    """
    if not present_set:
        return []
    intervals: List[Tuple[int, int]] = []
    sorted_recs = sorted(present_set)
    cur_start = sorted_recs[0]
    cur_end = sorted_recs[0]
    for r in sorted_recs[1:]:
        if r - cur_end <= 1 + gap_tolerance:
            cur_end = r
        else:
            intervals.append((cur_start, cur_end))
            cur_start = r
            cur_end = r
    intervals.append((cur_start, cur_end))
    return intervals


def compute_lifetimes(
    per_record_hbonds: Sequence[Sequence[HBond]],
    record_indices: Optional[Sequence[int]] = None,
    gap_tolerance: int = 0,
) -> List[PairLifetime]:
    """Compute lifetime stats per unique H-bond pair across a trajectory.

    Parameters
    ----------
    per_record_hbonds : list of (list of HBond), one per recorded frame
    record_indices : optional, absolute record numbers corresponding to
                     each entry of per_record_hbonds. If None, uses
                     ``list(range(len(per_record_hbonds)))``.
    gap_tolerance : max gap (in records) to bridge when computing
                    intermittent intervals. 0 = strict (= continuous).
    """
    n_records = len(per_record_hbonds)
    if record_indices is None:
        record_indices = list(range(n_records))

    # Bucket: pair_key -> set of records where pair appeared
    presence: Dict[PairKey, set] = defaultdict(set)
    for rec_idx, hbond_list in zip(record_indices, per_record_hbonds):
        for hb in hbond_list:
            presence[hbond_to_key(hb)].add(rec_idx)

    out: List[PairLifetime] = []
    for pair, recs in presence.items():
        strict_intervals = _runs(recs, n_records, gap_tolerance=0)
        strict_lengths = [e - s + 1 for s, e in strict_intervals]
        if gap_tolerance > 0:
            interm_intervals = _runs(recs, n_records, gap_tolerance=gap_tolerance)
            interm_lengths = [e - s + 1 for s, e in interm_intervals]
        else:
            interm_intervals = list(strict_intervals)
            interm_lengths = list(strict_lengths)
        out.append(PairLifetime(
            pair=pair,
            present_records=sorted(recs),
            continuous_intervals=strict_intervals,
            continuous_max=max(strict_lengths) if strict_lengths else 0,
            continuous_mean=(
                float(np.mean(strict_lengths)) if strict_lengths else 0.0
            ),
            intermittent_intervals=interm_intervals,
            intermittent_max=max(interm_lengths) if interm_lengths else 0,
            intermittent_mean=(
                float(np.mean(interm_lengths)) if interm_lengths else 0.0
            ),
            total_present=len(recs),
            occupancy=len(recs) / n_records if n_records else 0.0,
        ))
    return out


def compute_autocorrelation(
    per_record_hbonds: Sequence[Sequence[HBond]],
    max_lag: Optional[int] = None,
    normalize: bool = True,
) -> np.ndarray:
    """Luzar-Chandler population autocorrelation C(t).

    ``C(t) = <h(0)h(t)> / <h(0)>``

    where ``h(t) ∈ {0,1}`` indicates that the same pair (same key) is
    present at frame t. The average is over **pairs** and **time origins**.

    Parameters
    ----------
    per_record_hbonds : list of (list of HBond), one per frame
    max_lag : maximum lag (in frames) to compute. Default = N/2.
    normalize : if True, normalize so C(0) = 1.

    Returns
    -------
    C : np.ndarray of shape (max_lag + 1,)
    """
    n = len(per_record_hbonds)
    if max_lag is None:
        max_lag = n // 2

    # Build presence matrix: rows = pairs, cols = frames
    keys: set = set()
    for hbs in per_record_hbonds:
        for hb in hbs:
            keys.add(hbond_to_key(hb))
    key_to_idx = {k: i for i, k in enumerate(sorted(keys))}
    n_pairs = len(key_to_idx)
    if n_pairs == 0:
        return np.zeros(max_lag + 1)

    H = np.zeros((n_pairs, n), dtype=np.int8)
    for f, hbs in enumerate(per_record_hbonds):
        for hb in hbs:
            H[key_to_idx[hbond_to_key(hb)], f] = 1

    # C(t) = <h(0) h(t)> / <h(0)>  (Luzar-Chandler, unbiased estimator)
    # <h(0)h(t)> = (1/(N-t)) sum_{i=0}^{N-t-1} h(i) h(i+t)  averaged over pairs
    # <h(0)>    = (1/N) sum_i h(i)                          averaged over pairs
    c = np.zeros(max_lag + 1, dtype=float)
    for t in range(max_lag + 1):
        n_origins = n - t
        if n_origins <= 0:
            c[t] = 0.0
            continue
        # numerator averaged over n_origins time origins and n_pairs pairs
        num = float(np.sum(H[:, :n - t] * H[:, t:])) / (n_origins * n_pairs)
        c[t] = num
    if normalize and c[0] > 0:
        c = c / c[0]
    return c


def integrate_autocorrelation(
    c: np.ndarray, dt: float = 1.0,
) -> float:
    """Estimate τ_HB = ∫_0^∞ C(t) dt via simple trapezoidal rule.

    ``dt`` is the time step between consecutive records (in user units, e.g. ps).
    The result has the same time unit.
    """
    return float(np.trapz(c, dx=dt))


def summarize_lifetimes(lifetimes: List[PairLifetime]) -> dict:
    """Aggregate lifetime statistics across all detected pairs."""
    if not lifetimes:
        return {
            "n_unique_pairs": 0,
            "mean_occupancy": 0.0,
            "max_continuous": 0,
            "mean_continuous_max": 0.0,
            "mean_intermittent_max": 0.0,
        }
    occs = np.array([l.occupancy for l in lifetimes])
    cmaxes = np.array([l.continuous_max for l in lifetimes])
    imaxes = np.array([l.intermittent_max for l in lifetimes])
    return {
        "n_unique_pairs": len(lifetimes),
        "mean_occupancy": float(np.mean(occs)),
        "median_occupancy": float(np.median(occs)),
        "max_continuous": int(cmaxes.max()),
        "mean_continuous_max": float(np.mean(cmaxes)),
        "mean_intermittent_max": float(np.mean(imaxes)),
    }
