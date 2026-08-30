# -*- coding: utf-8 -*-
"""``_dedup_dihedrals`` must keep AMBER/GAFF multi-term torsions.

A GAFF torsion is a sum of Fourier terms that share the same four atoms and
differ only in (phase, k, multiplicity).  GROMACS writes them as one
``[ dihedrals ]`` line per term with function 9 ("proper dihedral, multiple").

The helper used to key on the atom quadruple alone, so every term after the
first was dropped as a duplicate -- silently flattening the torsion into
whichever term happened to be written first.  Measured on the aripiprazole /
lactic acid system, 39 of 268 terms were lost and the round-tripped topology
disagreed with the source for 35 of 229 quadruples.

Deduplication now also compares the function type and parameters, so only
genuinely identical entries collapse.
"""
import pytest

from abmptools.core.system_model import DihedralRecord
from abmptools.udf2gro.gromacs.writers.top_writer import _dedup_dihedrals


def d(atoms, funct, params):
    a1, a2, a3, a4 = atoms
    return DihedralRecord(atom1=a1, atom2=a2, atom3=a3, atom4=a4,
                          funct=funct, params=list(params))


class TestMultiTerm:

    def test_three_fourier_terms_survive(self):
        """The canonical GAFF case: one quadruple, three terms, function 9."""
        terms = [
            d((1, 2, 3, 4), "9", [180.0, 15.1670, 2]),
            d((1, 2, 3, 4), "9", [0.0, 1.2552, 1]),
            d((1, 2, 3, 4), "9", [180.0, 0.4184, 3]),
        ]
        out = _dedup_dihedrals(terms)
        assert len(out) == 3
        assert [x.params[2] for x in out] == [2, 1, 3]

    def test_terms_survive_when_written_reversed(self):
        """l-k-j-i is the same torsion, so its terms belong to the same set."""
        terms = [
            d((1, 2, 3, 4), "9", [180.0, 15.1670, 2]),
            d((4, 3, 2, 1), "9", [0.0, 1.2552, 1]),
        ]
        assert len(_dedup_dihedrals(terms)) == 2

    def test_identical_entries_still_collapse(self):
        """Same atoms, same funct, same params -- a real duplicate."""
        terms = [
            d((1, 2, 3, 4), "9", [180.0, 15.1670, 2]),
            d((4, 3, 2, 1), "9", [180.0, 15.1670, 2]),
            d((1, 2, 3, 4), "9", [0.0, 1.2552, 1]),
        ]
        out = _dedup_dihedrals(terms)
        assert len(out) == 2

    def test_self_dihedrals_are_dropped(self):
        """A repeated atom is a structural impossibility."""
        terms = [
            d((1, 2, 3, 1), "9", [180.0, 1.0, 2]),
            d((1, 2, 2, 4), "9", [180.0, 1.0, 2]),
            d((1, 2, 3, 4), "9", [180.0, 1.0, 2]),
        ]
        out = _dedup_dihedrals(terms)
        assert len(out) == 1
        assert (out[0].atom1, out[0].atom4) == (1, 4)

    def test_improper_and_proper_on_one_quadruple_both_kept(self):
        """Different function types are different interactions."""
        terms = [
            d((1, 2, 3, 4), "9", [180.0, 15.1670, 2]),
            d((1, 2, 3, 4), "4", [180.0, 4.6024, 2]),
        ]
        out = _dedup_dihedrals(terms)
        assert sorted(x.funct for x in out) == ["4", "9"]


class TestFunctPromotion:
    """GROMACS accepts only one function-1 entry per quadruple."""

    def test_multiple_proper_terms_are_promoted_to_9(self):
        terms = [
            d((1, 2, 3, 4), "1", [180.0, 15.1670, 2]),
            d((1, 2, 3, 4), "1", [0.0, 1.2552, 1]),
        ]
        out = _dedup_dihedrals(terms)
        assert len(out) == 2
        assert [x.funct for x in out] == ["9", "9"]
        # parameters must be carried over untouched
        assert out[0].params == [180.0, 15.1670, 2]
        assert out[1].params == [0.0, 1.2552, 1]

    def test_a_lone_function_1_entry_is_left_alone(self):
        terms = [
            d((1, 2, 3, 4), "1", [180.0, 15.1670, 2]),
            d((2, 3, 4, 5), "1", [0.0, 1.2552, 1]),
        ]
        out = _dedup_dihedrals(terms)
        assert [x.funct for x in out] == ["1", "1"]

    def test_promotion_does_not_mutate_the_input(self):
        terms = [
            d((1, 2, 3, 4), "1", [180.0, 15.1670, 2]),
            d((1, 2, 3, 4), "1", [0.0, 1.2552, 1]),
        ]
        _dedup_dihedrals(terms)
        assert [x.funct for x in terms] == ["1", "1"]


def test_order_is_preserved():
    """Terms keep their source order; GROMACS sums them either way, but a
    stable order keeps the written topology diffable."""
    terms = [d((1, 2, 3, 4), "9", [0.0, float(i), i]) for i in (1, 2, 3, 4)]
    out = _dedup_dihedrals(terms)
    assert [x.params[2] for x in out] == [1, 2, 3, 4]
