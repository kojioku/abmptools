"""Tests for the generic donor-type x acceptor-type pair statistics."""
from abmptools.hbond.hbond_detector import HBond
from abmptools.hbond.pair_type_stats import (
    GenericPairClassification, classify_generic, summarize_pair_stats,
)


def _hb(donor_mol, donor_d, donor_h, acceptor_mol, acceptor_a):
    return HBond(
        donor_mol=donor_mol, donor_d=donor_d, donor_h=donor_h,
        acceptor_mol=acceptor_mol, acceptor_a=acceptor_a,
        d_da=3.0, d_ha=2.0, angle=170.0,
    )


def test_no_hbonds_all_candidates():
    donors = {"hydroxyl": [(0, 3, 4), (1, 3, 4)]}
    acceptors = {"hydroxyl_O": [(0, 3), (1, 3)]}
    r = classify_generic(donors, acceptors, {})
    assert r.n_donor_candidates == {"hydroxyl": 2}
    assert r.n_acceptor_candidates == {"hydroxyl_O": 2}
    # all atom-roles should be Candidate
    assert all(v == "Candidate" for v in r.atom_role.values())
    assert r.pair_stats == {}


def test_single_pair_type():
    donors = {"hydroxyl": [(0, 3, 4)]}
    acceptors = {"hydroxyl_O": [(1, 3)]}
    hbs = {("hydroxyl", "hydroxyl_O"): [_hb(0, 3, 4, 1, 3)]}
    r = classify_generic(donors, acceptors, hbs)
    st = r.pair_stats[("hydroxyl", "hydroxyl_O")]
    assert st.n_hbonds == 1
    assert st.n_uniq_donors == 1
    assert st.n_uniq_acceptors == 1
    assert r.atom_role[(0, 3)] == "Donor"     # donor D atom
    assert r.atom_role[(0, 4)] == "Donor"     # donor H atom
    assert r.atom_role[(1, 3)] == "Acceptor"  # acceptor A atom


def test_both_role_when_donor_and_acceptor():
    """Same atom acting as donor AND acceptor (e.g. hydroxyl O) => Both."""
    # hydroxyl on mol 0 (atom 3 is the O = donor heavy, atom 4 = H,
    # but atom 3 also appears as acceptor in a separate hbond)
    donors = {"hydroxyl": [(0, 3, 4)]}
    acceptors = {"hydroxyl_O": [(0, 3), (1, 3)]}
    hbs = {
        ("hydroxyl", "hydroxyl_O"): [
            _hb(0, 3, 4, 1, 3),  # mol 0 donates to mol 1
            _hb(1, 3, 4, 0, 3) if False else _hb(2, 3, 4, 0, 3),
            # ↑ mol 2 donates to mol 0's hydroxyl O (so mol 0 atom 3 is acceptor)
        ],
    }
    # We did not list mol 2 as a donor in donors dict, so add it for completeness
    donors["hydroxyl"].append((2, 3, 4))
    r = classify_generic(donors, acceptors, hbs)
    # mol 0 atom 3 should be Both (donor in hb1, acceptor in hb2)
    assert r.atom_role[(0, 3)] == "Both"


def test_multiple_pair_types_summarize():
    donors = {
        "hydroxyl":   [(0, 3, 4), (1, 3, 4)],
        "amide_donor": [(2, 1, 5)],
    }
    acceptors = {
        "hydroxyl_O": [(0, 3), (1, 3)],
        "amide_O":    [(2, 2)],
    }
    hbs = {
        ("hydroxyl", "hydroxyl_O"): [_hb(0, 3, 4, 1, 3)],
        ("hydroxyl", "amide_O"):    [_hb(1, 3, 4, 2, 2)],
        ("amide_donor", "hydroxyl_O"): [_hb(2, 1, 5, 0, 3)],
    }
    r = classify_generic(donors, acceptors, hbs)
    summary = summarize_pair_stats(r)
    assert summary["n_hbonds_total"] == 3
    assert summary["hydroxyl->hydroxyl_O"] == 1
    assert summary["hydroxyl->amide_O"] == 1
    assert summary["amide_donor->hydroxyl_O"] == 1


def test_unique_counting_dedup():
    """If the same donor donates to two acceptors, counts unique donors only once."""
    donors = {"hydroxyl": [(0, 3, 4)]}
    acceptors = {"hydroxyl_O": [(1, 3), (2, 3)]}
    hbs = {
        ("hydroxyl", "hydroxyl_O"): [
            _hb(0, 3, 4, 1, 3),
            _hb(0, 3, 4, 2, 3),
        ],
    }
    r = classify_generic(donors, acceptors, hbs)
    st = r.pair_stats[("hydroxyl", "hydroxyl_O")]
    assert st.n_hbonds == 2
    assert st.n_uniq_donors == 1      # mol 0 only
    assert st.n_uniq_acceptors == 2   # mol 1, mol 2
