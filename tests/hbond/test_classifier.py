"""Tests for the per-functional-group H-bond classifier."""
from abmptools.hbond.classifier import (
    ClassificationResult, FunctionalGroupClassification, classify,
)
from abmptools.hbond.functional_groups import AmideGroup, CarboxylGroup
from abmptools.hbond.hbond_detector import HBond


def _hb(donor_mol, donor_d, donor_h, acceptor_mol, acceptor_a):
    return HBond(
        donor_mol=donor_mol, donor_d=donor_d, donor_h=donor_h,
        acceptor_mol=acceptor_mol, acceptor_a=acceptor_a,
        d_da=3.0, d_ha=2.0, angle=170.0,
    )


def _carb(mol, c=10, o=11, oh=0, ho=1):
    return CarboxylGroup(
        mol_index=mol, c_atom=c, o_atom=o, oh_atom=oh, ho_atom=ho,
    )


def _amide(mol, c=20, o=0, n=21):
    return AmideGroup(
        mol_index=mol, c_atom=c, o_atom=o, n_atom=n,
        tert=True, nh_atom=None,
    )


def test_no_hbonds():
    r = classify(5, [], [])
    assert r.n_dual_mols == 0
    assert r.n_single_mols == 0
    assert r.n_free_mols == 5
    assert r.n_carboxyls == 0
    assert r.n_amides == 0


def test_no_hbonds_with_carbs_all_free():
    carbs = [_carb(0), _carb(1), _carb(2)]
    r = classify(3, [], [], carboxyls=carbs)
    assert r.n_carboxyls == 3
    assert r.n_carboxyls_free == 3
    assert r.ratio_carboxyl_free == 1.0
    assert all(c.role == "free" for c in r.carboxyl_roles)


def test_single_only():
    """COOH on mol 0 donates to amide on mol 1."""
    carbs = [_carb(0)]
    amides = [_amide(1)]
    hb = _hb(0, 0, 1, 1, 0)
    r = classify(3, [], [hb], carbs, amides)
    assert r.carboxyl_roles[0].role == "single"
    assert r.n_carboxyls_single == 1
    assert r.ratio_carboxyl_single == 1.0
    assert r.amide_roles[0].role == "accept"
    assert r.n_amides_accept == 1
    assert r.ratio_amide_accept == 1.0
    # Mol-level: mol 0 has single COOH so mol-role=single; mol 1 has no COOH so free
    assert r.roles[0].role == "single"
    assert r.roles[1].role == "free"
    assert r.roles[2].role == "free"


def test_dual_pair():
    """Two COOHs from different mols form a cyclic dimer (both directions)."""
    carbs = [_carb(0), _carb(1)]
    cc = [_hb(0, 0, 1, 1, 11), _hb(1, 0, 1, 0, 11)]
    r = classify(3, cc, [], carbs)
    assert r.n_carboxyls == 2
    assert r.n_carboxyls_dual == 2
    assert r.ratio_carboxyl_dual == 1.0
    assert r.dual_pairs == [(0, 1)]
    assert r.n_dual_mols == 2
    assert r.n_free_mols == 1


def test_dual_precedence_over_single():
    """Dual takes precedence over single for a COOH appearing in both relations."""
    carbs = [_carb(0), _carb(1)]
    amides = [_amide(2)]
    cc = [_hb(0, 0, 1, 1, 11), _hb(1, 0, 1, 0, 11)]
    ca = [_hb(0, 0, 1, 2, 0)]
    r = classify(4, cc, ca, carbs, amides)
    assert r.carboxyl_roles[0].role == "dual"
    assert r.carboxyl_roles[1].role == "dual"
    assert r.roles[0].role == "dual"
    assert r.roles[1].role == "dual"
    # mol 2 has no COOH so its mol-role is free, but its amide is "accept"
    assert r.roles[2].role == "free"
    assert r.n_amides_accept == 1


def test_one_way_carb_carb_not_dual():
    """One-direction carb-carb is not dual; both COOHs are free."""
    carbs = [_carb(0), _carb(1)]
    cc = [_hb(0, 0, 1, 1, 11)]
    r = classify(3, cc, [], carbs)
    assert r.n_carboxyls_dual == 0
    assert all(c.role == "free" for c in r.carboxyl_roles)


def test_multi_cooh_per_mol():
    """One molecule with two independent COOH groups in different roles."""
    cooh1 = CarboxylGroup(mol_index=0, c_atom=10, o_atom=11, oh_atom=0, ho_atom=1)
    cooh2 = CarboxylGroup(mol_index=0, c_atom=20, o_atom=21, oh_atom=2, ho_atom=3)
    cooh3 = CarboxylGroup(mol_index=1, c_atom=10, o_atom=11, oh_atom=0, ho_atom=1)
    amide1 = _amide(2)
    carbs = [cooh1, cooh2, cooh3]
    amides = [amide1]
    # Dual: cooh1 (mol 0) ↔ cooh3 (mol 1)
    cc = [_hb(0, 0, 1, 1, 11), _hb(1, 0, 1, 0, 11)]
    # Single: cooh2 (mol 0) → amide1 (mol 2)
    ca = [_hb(0, 2, 3, 2, 0)]
    r = classify(3, cc, ca, carbs, amides)
    # Per-COOH roles - the key claim of the per-group design
    assert r.carboxyl_roles[0].role == "dual"      # cooh1
    assert r.carboxyl_roles[1].role == "single"    # cooh2
    assert r.carboxyl_roles[2].role == "dual"      # cooh3
    # Per-mol counts on mol 0: one COOH dual + one COOH single
    assert r.roles[0].n_carboxyls_dual == 1
    assert r.roles[0].n_carboxyls_single == 1
    assert r.roles[0].n_carboxyls_free == 0
    # Mol-level representative role: dual wins
    assert r.roles[0].role == "dual"
    assert r.roles[1].role == "dual"
    assert r.roles[2].role == "free"
    # Aggregate
    assert r.n_carboxyls == 3
    assert r.n_carboxyls_dual == 2
    assert r.n_carboxyls_single == 1
    assert r.n_carboxyls_free == 0
    assert r.n_amides == 1
    assert r.n_amides_accept == 1
    # Ratios
    assert r.ratio_carboxyl_dual == 2 / 3
    assert r.ratio_carboxyl_single == 1 / 3
    assert r.ratio_amide_accept == 1.0


def test_carboxyl_role_carries_partner_groups():
    """CarboxylRole.dual_partners records (partner_mol, partner_carb_idx)."""
    carbs = [_carb(0), _carb(1)]
    cc = [_hb(0, 0, 1, 1, 11), _hb(1, 0, 1, 0, 11)]
    r = classify(2, cc, [], carbs)
    assert r.carboxyl_roles[0].dual_partners == {(1, 0)}
    assert r.carboxyl_roles[1].dual_partners == {(0, 0)}


def test_amide_role_carries_donor_carboxyls():
    carbs = [_carb(0)]
    amides = [_amide(1)]
    ca = [_hb(0, 0, 1, 1, 0)]
    r = classify(2, [], ca, carbs, amides)
    assert r.amide_roles[0].donor_carboxyls == {(0, 0)}


def test_ratios_zero_when_no_groups():
    r = classify(2, [], [])
    assert r.ratio_carboxyl_dual == 0.0
    assert r.ratio_amide_accept == 0.0


def test_back_compat_classification_result_alias():
    """ClassificationResult is an alias for FunctionalGroupClassification."""
    assert ClassificationResult is FunctionalGroupClassification
