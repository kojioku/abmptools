"""Tests for H-bond lifetime + autocorrelation."""
import numpy as np

from abmptools.hbond.hbond_detector import HBond
from abmptools.hbond.lifetime import (
    compute_autocorrelation, compute_lifetimes,
    integrate_autocorrelation, summarize_lifetimes,
)


def _hb(donor_mol=0, acceptor_mol=1):
    return HBond(
        donor_mol=donor_mol, donor_d=3, donor_h=40,
        acceptor_mol=acceptor_mol, acceptor_a=4,
        d_da=2.8, d_ha=1.9, angle=170.0,
    )


def test_lifetime_continuous_simple():
    """Pair present in 3 consecutive frames → continuous_max=3."""
    per_rec = [[_hb()], [_hb()], [_hb()], [], []]
    lt = compute_lifetimes(per_rec)
    assert len(lt) == 1
    assert lt[0].total_present == 3
    assert lt[0].continuous_max == 3
    assert lt[0].continuous_mean == 3.0
    assert lt[0].occupancy == 3 / 5


def test_lifetime_two_intervals():
    """Frames 0-1 present, 2 absent, 3 present → two intervals strict."""
    per_rec = [[_hb()], [_hb()], [], [_hb()]]
    lt = compute_lifetimes(per_rec)
    assert lt[0].continuous_intervals == [(0, 1), (3, 3)]
    assert lt[0].continuous_max == 2
    assert lt[0].total_present == 3


def test_lifetime_gap_tolerance_bridges():
    """gap_tolerance=1 bridges single-frame gap."""
    per_rec = [[_hb()], [_hb()], [], [_hb()]]
    lt = compute_lifetimes(per_rec, gap_tolerance=1)
    assert lt[0].intermittent_intervals == [(0, 3)]
    assert lt[0].intermittent_max == 4
    # strict still reports two intervals
    assert lt[0].continuous_intervals == [(0, 1), (3, 3)]


def test_lifetime_multiple_pairs():
    pair_a = _hb(donor_mol=0, acceptor_mol=1)
    pair_b = _hb(donor_mol=2, acceptor_mol=3)
    per_rec = [[pair_a, pair_b], [pair_a], [pair_b]]
    lt = compute_lifetimes(per_rec)
    assert len(lt) == 2
    occs = {l.pair[0]: l.occupancy for l in lt}
    assert occs[0] == 2 / 3
    assert occs[2] == 2 / 3


def test_autocorrelation_normalized_at_zero():
    per_rec = [[_hb()], [_hb()], [_hb()]]
    c = compute_autocorrelation(per_rec, normalize=True)
    assert abs(c[0] - 1.0) < 1e-9
    # Constant pair: C(t) = 1 for all t (within max_lag)
    assert all(abs(ct - 1.0) < 1e-9 for ct in c)


def test_autocorrelation_decays_with_absence():
    """Pair present at t=0,1 but absent at t=2 → C(2) < C(1)."""
    per_rec = [[_hb()], [_hb()], []]
    c = compute_autocorrelation(per_rec, max_lag=2, normalize=True)
    assert c[0] == 1.0
    assert c[1] > 0
    assert c[2] == 0  # only h(0)*h(2) = 1*0 = 0


def test_integrate_autocorrelation():
    """∫_0^T 1 dt = T."""
    c = np.array([1.0, 1.0, 1.0, 1.0])  # constant
    assert abs(integrate_autocorrelation(c, dt=1.0) - 3.0) < 1e-9
    assert abs(integrate_autocorrelation(c, dt=0.5) - 1.5) < 1e-9


def test_summarize_empty():
    s = summarize_lifetimes([])
    assert s["n_unique_pairs"] == 0
    assert s["mean_occupancy"] == 0.0
