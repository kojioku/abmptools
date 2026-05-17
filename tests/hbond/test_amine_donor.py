"""Tests for amine/amide N-H donor detection (secondary amide support).

Builds a synthetic peptide-like topology in-memory (no IMC dependency).
"""
from abmptools.hbond.bdf_reader import AtomInfo, BondInfo, MoleculeTopology
from abmptools.hbond.func_tags import CHARMM36, GAFF2
from abmptools.hbond.functional_groups import (
    detect_amides, detect_amine_donors,
)


def _make_secondary_amide_mol():
    """Build a tiny -C(=O)-N(H)-CH3 molecule in GAFF2 style.

    Atom layout (0-based local indices):
      0: c   (carbonyl C)
      1: o   (carbonyl O)
      2: n   (amide N) -- 2° amide (has 1 H)
      3: hn  (N-H)
      4: c3  (methyl C on N)
      5,6,7: hc (methyl H's)
      8: c3 (CH3 on C side)
      9,10,11: hc
    """
    atoms = [
        AtomInfo(atom_id=0,  atom_name="C", atom_type="c"),
        AtomInfo(atom_id=1,  atom_name="O", atom_type="o"),
        AtomInfo(atom_id=2,  atom_name="N", atom_type="n"),
        AtomInfo(atom_id=3,  atom_name="H", atom_type="hn"),
        AtomInfo(atom_id=4,  atom_name="C", atom_type="c3"),
        AtomInfo(atom_id=5,  atom_name="H", atom_type="hc"),
        AtomInfo(atom_id=6,  atom_name="H", atom_type="hc"),
        AtomInfo(atom_id=7,  atom_name="H", atom_type="hc"),
        AtomInfo(atom_id=8,  atom_name="C", atom_type="c3"),
        AtomInfo(atom_id=9,  atom_name="H", atom_type="hc"),
        AtomInfo(atom_id=10, atom_name="H", atom_type="hc"),
        AtomInfo(atom_id=11, atom_name="H", atom_type="hc"),
    ]
    bonds = [
        BondInfo(type_name="c-o", atom1=0, atom2=1, order=2.0),
        BondInfo(type_name="c-n", atom1=0, atom2=2, order=1.0),
        BondInfo(type_name="n-hn", atom1=2, atom2=3, order=1.0),
        BondInfo(type_name="n-c3", atom1=2, atom2=4, order=1.0),
        BondInfo(type_name="c3-hc", atom1=4, atom2=5, order=1.0),
        BondInfo(type_name="c3-hc", atom1=4, atom2=6, order=1.0),
        BondInfo(type_name="c3-hc", atom1=4, atom2=7, order=1.0),
        BondInfo(type_name="c-c3", atom1=0, atom2=8, order=1.0),
        BondInfo(type_name="c3-hc", atom1=8, atom2=9, order=1.0),
        BondInfo(type_name="c3-hc", atom1=8, atom2=10, order=1.0),
        BondInfo(type_name="c3-hc", atom1=8, atom2=11, order=1.0),
    ]
    return MoleculeTopology(mol_name="N-methylacetamide", atoms=atoms, bonds=bonds)


def test_detect_amides_includes_sec_amide():
    mol = _make_secondary_amide_mol()
    amides = detect_amides([mol], mapping=GAFF2)
    assert len(amides) == 1
    a = amides[0]
    assert a.tert is False
    assert a.nh_atom == 3
    assert a.c_atom == 0
    assert a.n_atom == 2


def test_detect_amine_donors_sec_amide():
    mol = _make_secondary_amide_mol()
    donors = detect_amine_donors([mol], mapping=GAFF2)
    assert len(donors) == 1
    d = donors[0]
    assert d.from_amide is True
    assert d.n_atom == 2
    assert d.h_atom == 3


def test_amine_donor_filter_amide_only():
    """include_amide=False excludes secondary amide N-H donors."""
    mol = _make_secondary_amide_mol()
    no_amide = detect_amine_donors(
        [mol], include_amide=False, include_amine=True, mapping=GAFF2,
    )
    assert len(no_amide) == 0


def test_imc_tert_amide_has_no_nh_donor():
    """IMC is tertiary amide → no N-H donor."""
    import os
    bdf = "/home/okuwaki/llm-project/SI/IMC_result450.0_out_rec900.bdf"
    if not os.path.exists(bdf):
        return  # skip
    from abmptools.hbond.bdf_reader import BDFTrajectory
    traj = BDFTrajectory(bdf)
    traj.load_topology()
    donors = detect_amine_donors(traj.molecules)
    assert len(donors) == 0


def test_charmm_secondary_amide():
    """Same topology but with CHARMM36 atom types."""
    mol = _make_secondary_amide_mol()
    # Re-tag atoms for CHARMM36
    mol.atoms[0].atom_type = "C"     # carbonyl C
    mol.atoms[1].atom_type = "O"     # carbonyl O
    mol.atoms[2].atom_type = "NH1"   # secondary amide N
    mol.atoms[3].atom_type = "H"     # N-H
    for i in [4, 8]:
        mol.atoms[i].atom_type = "CT3"
    for i in [5, 6, 7, 9, 10, 11]:
        mol.atoms[i].atom_type = "HA3"
    amides = detect_amides([mol], mapping=CHARMM36)
    assert len(amides) == 1
    assert amides[0].tert is False
    donors = detect_amine_donors([mol], mapping=CHARMM36)
    assert len(donors) == 1
