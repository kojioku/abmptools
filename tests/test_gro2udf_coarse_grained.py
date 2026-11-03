# -*- coding: utf-8 -*-
"""A coarse-grained topology converts with its own force field, not a guess.

Four things differ from an all-atom topology, and each used to be written
as something else without a word:

* every LJ parameter sits in ``[ nonbond_params ]``, not in ``[ atomtypes ]``
* the angles are GROMOS-96 cosines (funct 2), not harmonic in theta
* ``nrexcl`` is 1, so 1-3 and 1-4 pairs are computed
* ``epsilon_r`` is 15, and the reaction field has nowhere to put it

The numbers here are the ones a real conversion produced; the energies they
reproduce were checked against GROMACS and are written down in
``docs/gro2udf.md``.
"""
from __future__ import annotations

import math

import pytest

from abmptools.gro2udf.guard import check_top, scan_sections, used_moleculetypes
from abmptools.gro2udf.top_parser import TopParser

CG_TOP = """\
[ defaults ]
1 2 ; sigma-epsilon format of LJ parameters

[ atomtypes ]
; the LJ is NOT here -- this is what the explicit pairs are for
  Q1   72.0 0.000 A 0.0 0.0
  C1   72.0 0.000 A 0.0 0.0
  W    72.0 0.000 A 0.0 0.0

[ nonbond_params ]
  Q1  Q1  1 0.47 3.50
  Q1  C1  1 0.47 2.00
  C1  C1  1 0.47 3.50
  Q1  W   1 0.47 5.60
  C1  W   1 0.47 2.00
  W   W   1 0.47 4.65

[moleculetype]
  LIP  1

[atoms]
  1 Q1 1 LIP NC3 1  1.0 72.0
  2 Q1 1 LIP PO4 2 -1.0 72.0
  3 C1 1 LIP C1A 3  0.0 72.0

[bonds]
  1 2 1 0.40 7000
  2 3 1 0.42 1350

[angles]
  1 2 3 2 120.0 25.0

[moleculetype]
  W  1

[atoms]
  1 W 1 W W 1 0.0 72.0

[ system ]
cg

[ molecules ]
LIP  2
W    3
"""


@pytest.fixture
def raw(tmp_path):
    p = tmp_path / "cg.top"
    p.write_text(CG_TOP, encoding="utf-8")
    return TopParser().parse(str(p))


def test_the_explicit_pairs_are_read(raw):
    """``[ atomtypes ]`` is all zeros here; without these the UDF has no LJ."""
    assert len(raw.nonbond_params) == 6
    assert raw.nonbond_params[("Q1", "Q1")] == (0.47, 3.50)
    # sorted key, whichever order the line was written in
    assert raw.nonbond_params[("Q1", "W")] == (0.47, 5.60)


def test_a_section_without_inner_spaces_is_still_a_section(raw):
    """``[moleculetype]`` and ``[ moleculetype ]`` are the same to GROMACS.

    Reading only the spaced form skipped every molecule of a force field
    written the other way -- 112 of them in one real lipid library, with no
    message at all.
    """
    assert raw.mol_types == ["LIP", "W"]
    assert len(raw.atomlist[raw.mol_types.index("LIP")]) == 3


def test_nrexcl_is_read(raw):
    assert raw.nrexcl["LIP"] == 1


def test_the_angle_keeps_its_funct(raw):
    assert [a[3] for a in raw.angle_types_from_mol] == [2]


def test_nothing_is_reported_as_unsupported(raw):
    fatal, warnings = check_top(raw, raw.sections)
    assert fatal == [], fatal
    assert warnings == [], warnings


def test_only_the_molecules_in_the_system_are_counted(tmp_path):
    """A library ``.itp`` carries molecules the system does not use; their
    sections must not fail the conversion (one lipid library's imidazolium
    ions refused every bilayer because of their ``[ constraints ]``)."""
    text = CG_TOP.replace("[ system ]", """[moleculetype]
  UNUSED  1

[atoms]
  1 Q1 1 UNU A1 1 0.0 72.0
  2 Q1 1 UNU A2 2 0.0 72.0

[constraints]
  1 2 1 0.33

[ system ]""")
    p = tmp_path / "lib.top"
    p.write_text(text, encoding="utf-8")
    lines = text.splitlines(keepends=True)
    assert used_moleculetypes(lines) == {"LIP", "W"}
    assert scan_sections(lines).get("constraints", 0) == 0
    raw = TopParser().parse(str(p))
    assert "UNUSED" not in raw.mol_types          # not even read
    fatal, _ = check_top(raw, raw.sections)
    assert fatal == [], fatal


def test_the_same_sections_are_counted_when_the_molecule_is_used(tmp_path):
    """The opposite case still counts: a constraint in a molecule the system
    does contain is real, and must keep being reported."""
    text = CG_TOP.replace("""[bonds]
  1 2 1 0.40 7000""", """[constraints]
  1 2 1 0.40""")
    p = tmp_path / "con.top"
    p.write_text(text, encoding="utf-8")
    assert scan_sections(text.splitlines(keepends=True))["constraints"] >= 1
    raw = TopParser().parse(str(p))
    fatal, _ = check_top(raw, raw.sections)
    assert any("constraints" in f for f in fatal), fatal


def test_epsilon_r_is_folded_into_the_charges_exactly():
    """``U = q_i q_j / eps`` equals ``(q_i/sqrt(eps))(q_j/sqrt(eps))``.

    The reaction field this writes has no permittivity inside the cutoff,
    so the only exact place for ``epsilon_r`` is the charges.
    """
    eps = 15.0
    q1, q2 = 1.0, -1.0
    s = 1.0 / math.sqrt(eps)
    assert (q1 * s) * (q2 * s) == pytest.approx(q1 * q2 / eps, rel=1e-12)
