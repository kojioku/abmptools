# -*- coding: utf-8 -*-
"""`[ atomtypes ]` carries the atomic number, so a reader can tell the element.

GROMACS does not need it — it works out what it needs from the parameters —
but an importer does. The names in the column are force field types ("c3",
"os", "hc"), which mean nothing outside that force field, and the atom names
in `[ atoms ]` are derived from them. Without an atomic number the only clue
left is the mass, and a reader that does not guess from mass has no way to
distinguish an all-atom system from coarse-grained beads: J-OCTA's
import_gromacs read a converted all-atom system as CG (2026-09-15, J-OCTA
11.0).

A coarse-grained bead genuinely has no element, and 0 — GROMACS' own value
for "unknown" — is the right answer there, so this is not a case of always
filling something in.
"""
from __future__ import annotations

from abmptools.udf2gro.gromacs.writers.top_writer import _atomic_number


def test_the_elements_a_force_field_actually_ships():
    """Rounded masses, as force fields write them."""
    assert _atomic_number(1.008) == 1
    assert _atomic_number(12.01) == 6      # GAFF carbon
    assert _atomic_number(12.011) == 6
    assert _atomic_number(14.01) == 7
    assert _atomic_number(16.0) == 8       # GAFF oxygen, rounded
    assert _atomic_number(15.999) == 8
    assert _atomic_number(32.06) == 16
    assert _atomic_number(35.45) == 17


def test_a_bead_gets_zero_rather_than_a_guess():
    """A CG bead has no element; 0 says so.

    Martini beads sit around 72 amu, nowhere near an element, and inventing
    one would turn "unknown" into a wrong answer.
    """
    assert _atomic_number(72.0) == 0
    assert _atomic_number(0.0) == 0
    assert _atomic_number(45.0) == 0


def test_neighbouring_elements_are_not_confused():
    """The tolerance must not let one element answer for the next."""
    assert _atomic_number(14.007) == 7
    assert _atomic_number(15.999) == 8
    assert _atomic_number(18.998) == 9
    assert _atomic_number(19.5) == 0       # 間はどちらでもない
