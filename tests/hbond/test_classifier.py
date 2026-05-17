"""Tests for the dual/single/free classifier."""
from abmptools.hbond.classifier import classify
from abmptools.hbond.hbond_detector import HBond


def _hb(donor_mol, acceptor_mol):
    return HBond(
        donor_mol=donor_mol, donor_d=0, donor_h=1,
        acceptor_mol=acceptor_mol, acceptor_a=0,
        d_da=3.0, d_ha=2.0, angle=170.0
    )


def test_no_hbonds():
    r = classify(5, [], [])
    assert r.n_dual_mols == 0
    assert r.n_single_mols == 0
    assert r.n_free_mols == 5


def test_single_only():
    # mol 0 donates to mol 1 (carb→amide)
    r = classify(3, [], [_hb(0, 1)])
    assert r.n_dual_mols == 0
    assert r.n_single_mols == 2     # 0 (donor) and 1 (acceptor)
    assert r.n_free_mols == 1
    assert r.roles[0].role == "single"
    assert r.roles[1].role == "single"
    assert r.roles[2].role == "free"


def test_dual_pair():
    # mol 0 <-> mol 1 dual carb-carb
    r = classify(3, [_hb(0, 1), _hb(1, 0)], [])
    assert r.n_dual_mols == 2
    assert r.n_single_mols == 0
    assert r.n_free_mols == 1
    assert r.dual_pairs == [(0, 1)]


def test_dual_precedence_over_single():
    # mol 0 forms dual with mol 1, AND mol 0 also single-bonds to mol 2.
    # Priority: dual > single, so 0 is "dual".
    r = classify(4, [_hb(0, 1), _hb(1, 0)], [_hb(0, 2)])
    assert r.roles[0].role == "dual"
    assert r.roles[1].role == "dual"
    # mol 2 only participated in single (acceptor) -> single
    assert r.roles[2].role == "single"
    assert r.roles[3].role == "free"


def test_one_way_carb_carb_not_dual():
    """Only one direction of carb-carb is NOT a dual pair."""
    r = classify(3, [_hb(0, 1)], [])
    assert r.n_dual_mols == 0
    # one-way carb-carb still leaves both mols as "free" (not classified as single
    # because classifier looks for carb→amide for single)
    assert r.n_free_mols == 3
