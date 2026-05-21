"""Tests for the element + bond-graph fallback tagging.

Covers atoms whose ``atom_type`` is None or a per-atom unique string
(OpenFF SMIRNOFF case): the fallback must still tag O/H/N/C correctly
so that ``detect_carboxyls`` / ``detect_amides`` / ``detect_hydroxyls``
keep working.
"""
from abmptools.hbond.bdf_reader import AtomInfo, BondInfo, MoleculeTopology
from abmptools.hbond.func_tags import (
    GAFF2, TAG_AMIDE_H, TAG_AMIDE_N, TAG_CARBONYL_C, TAG_CARBONYL_O,
    TAG_HYDROXYL_H, TAG_HYDROXYL_O, fallback_tag_by_element,
)
from abmptools.hbond.functional_groups import (
    detect_amides, detect_carboxyls, detect_hydroxyls,
)
from abmptools.hbond import functional_groups as fg


def _bond(a, b):
    return BondInfo(type_name="", atom1=a, atom2=b, order=1.0)


def _atom(name, atom_id, atom_type=None):
    return AtomInfo(atom_id=atom_id, atom_name=name, atom_type=atom_type)


def test_fallback_tags_alcohol_oh():
    """Methanol (CH3-OH): O is hydroxyl_O, H is hydroxyl_H, C is sp3 (no tag)."""
    names = ["C", "O", "H1", "H2", "H3", "HO"]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)]
    initial = [None] * 6
    tags = fallback_tag_by_element(names, bonds, initial)
    assert tags[1] == TAG_HYDROXYL_O   # O
    assert tags[5] == TAG_HYDROXYL_H   # HO


def test_fallback_tags_carboxyl_pattern():
    """Acetic acid (CH3-COOH): C(=O)-O-H is the carboxyl carbon."""
    # atom indices: 0=CH3, 1=C(=O), 2=O(=C, carbonyl), 3=O(H), 4=H(O), 5..7=CH3-H
    names = ["C", "C", "O", "O", "HO", "H1", "H2", "H3"]
    bonds = [
        (0, 1),  # C(methyl) - C(carbonyl)
        (1, 2),  # C(carbonyl) = O
        (1, 3),  # C(carbonyl) - O(H)
        (3, 4),  # O - H
        (0, 5), (0, 6), (0, 7),
    ]
    tags = fallback_tag_by_element(names, bonds, [None] * 8)
    assert tags[2] == TAG_CARBONYL_O    # O without H
    assert tags[3] == TAG_HYDROXYL_O    # O with H
    assert tags[4] == TAG_HYDROXYL_H    # H bonded to O
    assert tags[1] == TAG_CARBONYL_C    # C bonded to carbonyl_O


def test_fallback_tags_amide_NH():
    """Acetamide (CH3-CO-NH2): N tagged amide_N, N-H tagged amide_H."""
    # atoms: 0=CH3 carbon, 1=C(=O), 2=O, 3=N, 4=H(N), 5=H(N), 6-8=CH3-H
    names = ["C", "C", "O", "N", "HN1", "HN2", "H1", "H2", "H3"]
    bonds = [
        (0, 1), (1, 2), (1, 3), (3, 4), (3, 5),
        (0, 6), (0, 7), (0, 8),
    ]
    tags = fallback_tag_by_element(names, bonds, [None] * 9)
    assert tags[2] == TAG_CARBONYL_O
    assert tags[1] == TAG_CARBONYL_C
    assert tags[3] == TAG_AMIDE_N
    assert tags[4] == TAG_AMIDE_H
    assert tags[5] == TAG_AMIDE_H


def test_fallback_preserves_existing_tags():
    """Atoms whose mapping already returned a tag are not overwritten."""
    names = ["O", "H"]
    bonds = [(0, 1)]
    initial = [TAG_CARBONYL_O, None]   # O pre-tagged as carbonyl (no H per ff)
    tags = fallback_tag_by_element(names, bonds, initial)
    # initial carbonyl_O should be preserved, H still tagged
    assert tags[0] == TAG_CARBONYL_O
    assert tags[1] == TAG_HYDROXYL_H


def test_detect_carboxyls_with_unknown_atom_type():
    """End-to-end: acetic acid topology with atom_type=None should yield 1 COOH."""
    fg.USE_ELEMENT_FALLBACK = True
    atoms = [
        _atom("C",  0),  # 0: CH3 C
        _atom("C",  1),  # 1: COOH C
        _atom("O",  2),  # 2: =O
        _atom("O",  3),  # 3: -O-H
        _atom("HO", 4),  # 4: -O-H
        _atom("H1", 5),
        _atom("H2", 6),
        _atom("H3", 7),
    ]
    bonds = [
        _bond(0, 1), _bond(1, 2), _bond(1, 3), _bond(3, 4),
        _bond(0, 5), _bond(0, 6), _bond(0, 7),
    ]
    topo = MoleculeTopology(mol_name="X", atoms=atoms, bonds=bonds)
    carboxyls = detect_carboxyls([topo], mapping=GAFF2)
    assert len(carboxyls) == 1
    cg = carboxyls[0]
    assert cg.c_atom == 1
    assert cg.o_atom == 2
    assert cg.oh_atom == 3
    assert cg.ho_atom == 4


def test_detect_hydroxyls_with_unknown_atom_type():
    """Methanol with atom_type=None should yield 1 alcohol OH."""
    fg.USE_ELEMENT_FALLBACK = True
    atoms = [
        _atom("C",  0),
        _atom("O",  1),
        _atom("HO", 2),
        _atom("H1", 3),
        _atom("H2", 4),
        _atom("H3", 5),
    ]
    bonds = [
        _bond(0, 1), _bond(1, 2),
        _bond(0, 3), _bond(0, 4), _bond(0, 5),
    ]
    topo = MoleculeTopology(mol_name="X", atoms=atoms, bonds=bonds)
    hydroxyls = detect_hydroxyls([topo], mapping=GAFF2)
    assert len(hydroxyls) == 1
    assert hydroxyls[0].oh_atom == 1
    assert hydroxyls[0].ho_atom == 2


def test_fallback_disabled_yields_no_groups():
    """With USE_ELEMENT_FALLBACK=False, unknown atom_type gives no detections."""
    saved = fg.USE_ELEMENT_FALLBACK
    try:
        fg.USE_ELEMENT_FALLBACK = False
        atoms = [_atom("C", 0), _atom("O", 1), _atom("HO", 2)]
        bonds = [_bond(0, 1), _bond(1, 2)]
        topo = MoleculeTopology(mol_name="X", atoms=atoms, bonds=bonds)
        # all atom_type are None, fallback off → no hydroxyls
        assert detect_hydroxyls([topo], mapping=GAFF2) == []
    finally:
        fg.USE_ELEMENT_FALLBACK = saved
