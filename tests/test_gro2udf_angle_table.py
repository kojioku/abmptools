# -*- coding: utf-8 -*-
"""Restricted bending as a table, and explicit exclusions as bonds with no k.

Both were refused until now because the engine has no matching form. Neither
replacement changes the physics:

* the table **is** the function, sampled; the sampling error is measured here
  and the engine reproduced GROMACS to six figures on a real three-bead
  system (2.82845 against 2.82845 at 100 degrees, 9.40145 against 9.40144 at
  160);
* a bond with no force constant adds no energy and removes exactly the one
  pair, which was measured against GROMACS' own ``[ exclusions ]``.
"""
from __future__ import annotations

import math
import pathlib

import pytest

from abmptools.gro2udf.angle_table import (DEFAULT_DIVISIONS, reb_energy,
                                           table_values, write_angle_table)
from abmptools.gro2udf.guard import check_top
from abmptools.gro2udf.top_adapter import TopAdapter
from abmptools.gro2udf.top_model import EXCLUSION_BOND_PREFIX
from abmptools.gro2udf.top_parser import TopParser

K, THETA0 = 25.0, 130.0


def _gromacs_reb(theta_deg: float, k: float = K, t0: float = THETA0) -> float:
    """V = k/2 (cos t - cos t0)^2 / sin^2 t, straight from the definition."""
    t = math.radians(theta_deg)
    return 0.5 * k * (math.cos(t) - math.cos(math.radians(t0))) ** 2 / math.sin(t) ** 2


@pytest.mark.parametrize("theta", [60.0, 90.0, 100.0, 130.0, 150.0, 160.0, 175.0])
def test_the_table_variable_is_the_exterior_cosine(theta):
    """The engine indexes the table by ``cos`` of the **exterior** angle.

    Reading it as the interior one mirrors the whole potential, which is
    invisible for a symmetric t0 and badly wrong for any other.
    """
    var = -math.cos(math.radians(theta))     # exterior
    assert reb_energy(var, K, THETA0) == pytest.approx(_gromacs_reb(theta),
                                                       rel=1e-12, abs=1e-12)


def test_the_minimum_sits_at_theta0():
    assert reb_energy(-math.cos(math.radians(THETA0)), K, THETA0) == pytest.approx(0.0)


def test_how_much_the_sampling_costs():
    """The only error a table adds is between its points. Measure it.

    Measured the hard way -- halfway between two points, with straight-line
    interpolation, which is the worst case; the engine is asked for 4th
    order, so it does better than this. The error is quoted in kJ/mol, not
    as a fraction: the potential passes through zero at t0, where any
    fraction is meaningless.

    Over 20 to 170 degrees the worst is 0.054 kJ/mol, at 170, where the
    potential is already 48 kJ/mol; from 40 to 150 it is below 0.001. Past
    170 the potential diverges (1305 kJ/mol at 178 degrees) and no run goes
    there.
    """
    values = table_values(K, THETA0, DEFAULT_DIVISIONS)
    step = 2.0 / DEFAULT_DIVISIONS
    worst, worst_at = 0.0, None
    for theta_i in range(20, 171):
        var = -math.cos(math.radians(theta_i))
        i = min(max(int((var + 1.0) / step), 0), DEFAULT_DIVISIONS - 1)
        mid = -1.0 + (i + 0.5) * step
        err = abs(0.5 * (values[i] + values[i + 1]) - reb_energy(mid, K, THETA0))
        if err > worst:
            worst, worst_at = err, theta_i
    assert worst < 0.1, (worst, worst_at)        # kJ/mol, against kT = 2.5


def test_the_divergent_end_is_written_finite():
    """``sin t`` reaches zero at 0 and 180 degrees. A table cannot hold
    infinity, and a run that reached it would already be broken."""
    values = table_values(K, THETA0, 200)
    assert math.isfinite(values[0]) and math.isfinite(values[-1])
    assert values[0] > 1.0e5 and values[-1] > 1.0e5


def test_the_table_file_says_what_it_is(tmp_path):
    p = write_angle_table(tmp_path / "t.udf", "reb", K, THETA0, divisions=200)
    text = p.read_text(encoding="utf-8")
    assert "Mesh:{" in text and "FieldValue:{" in text
    assert "[ -1.0 1.0 200 ]" in text
    assert text.count("\n          ") == 201          # divisions + 1 values
    assert "Unit_Parameter" in text


EXCL_TOP = """\
[ defaults ]
1 2

[ atomtypes ]
  C1   72.0 0.000 A 0.0 0.0

[ nonbond_params ]
  C1  C1  1 0.47 3.50

[moleculetype]
  REB  1

[atoms]
  1 C1 1 REB A1 1 0.0 72.0
  2 C1 1 REB B1 2 0.0 72.0
  3 C1 1 REB C1 3 0.0 72.0

[bonds]
  1 2 1 0.47 1250
  2 3 1 0.47 1250

[angles]
  1 2 3 10 130.0 25.0

[exclusions]
  1 3

[ system ]
reb

[ molecules ]
REB 1
"""


@pytest.fixture
def raw(tmp_path):
    p = tmp_path / "reb.top"
    p.write_text(EXCL_TOP, encoding="utf-8")
    return TopParser().parse(str(p))


def test_an_exclusion_becomes_a_bond_with_no_force_constant(raw):
    assert raw.exclusions["REB"] == [(1, 3)]
    i = raw.mol_types.index("REB")
    assert [b[:2] for b in raw.bondlist[i]] == [[1, 2], [2, 3], [1, 3]]
    k = [bt[3][1] for bt in raw.bond_types_from_mol]
    assert 0.0 in k, k


def test_that_bond_is_named_so_it_is_never_constrained(raw, tmp_path):
    """RATTLE takes bond potentials by name; one of these in the list would
    fix two atoms at a length that is not a bond length."""
    gro = tmp_path / "x.gro"
    gro.write_text("x\n    3\n"
                   "    1REB     A1    1   2.000   2.000   2.000\n"
                   "    1REB     B1    2   2.470   2.000   2.000\n"
                   "    1REB     C1    3   2.800   2.340   2.000\n"
                   "   6.00000   6.00000   6.00000\n", encoding="utf-8")
    model = TopAdapter().build(raw, str(gro))
    excl = [b for b in model.bond_type_specs
            if b.name.startswith(EXCLUSION_BOND_PREFIX)]
    assert len(excl) == 1 and excl[0].kb == 0.0 and excl[0].r0 == 0.0
    # and the real bond is not marked
    assert any(not b.name.startswith(EXCLUSION_BOND_PREFIX)
               and b.kb > 0 for b in model.bond_type_specs)


def test_restricted_bending_and_an_exclusion_both_pass_the_guard(raw):
    fatal, warnings = check_top(raw, raw.sections)
    assert fatal == [], fatal


def test_an_exclusion_with_a_larger_nrexcl_is_refused(tmp_path):
    """The trick only removes the one pair while nrexcl is 1."""
    p = tmp_path / "n3.top"
    p.write_text(EXCL_TOP.replace("  REB  1\n", "  REB  3\n"), encoding="utf-8")
    raw = TopParser().parse(str(p))
    fatal, _ = check_top(raw, raw.sections)
    assert any("nrexcl >= 2" in f for f in fatal), fatal


# --- a file written to be looked at ------------------------------------------

VIEW_TOP = EXCL_TOP.replace("""[angles]
  1 2 3 10 130.0 25.0""", """[angles]
  1 2 3 10 130.0 25.0

[dihedrals]
  1 2 3 1 1 -120.0 400.0 1

[virtual_sitesn]
  3 2 1 2""")


def test_a_viewing_file_is_written_and_says_what_it_lacks(tmp_path):
    """The point of such a file is to be opened, so it must not be refused.

    What it cannot reproduce is written down instead: the caller gets it in
    ``viewing_notes``, and a file beside the UDF carries the same words so
    the two cannot be separated.
    """
    from abmptools.gro2udf.top_exporter import TopExporter

    top = tmp_path / "v.top"
    top.write_text(VIEW_TOP, encoding="utf-8")
    gro = tmp_path / "v.gro"
    gro.write_text("v\n    3\n"
                   "    1REB     A1    1   2.000   2.000   2.000\n"
                   "    1REB     B1    2   2.470   2.000   2.000\n"
                   "    1REB     C1    3   2.800   2.340   2.000\n"
                   "   6.00000   6.00000   6.00000\n", encoding="utf-8")
    out = tmp_path / "v.udf"
    tpl = pathlib.Path(__file__).resolve().parents[1] / \
        "abmptools" / "gro2udf" / "default_template.udf"
    x = TopExporter()
    x.export(str(top), str(gro), str(tpl), str(out), for_viewing=True)
    assert out.is_file()

    joined = "\n".join(x.viewing_notes)
    assert "phase" in joined                      # the dihedral it left out
    assert "virtual site" in joined               # and what it did with those
    notes = pathlib.Path(str(out) + ".notes.txt")
    assert notes.is_file()
    text = notes.read_text(encoding="utf-8")
    assert "TO BE LOOKED AT" in text
    assert "phase" in text


def test_the_same_topology_is_still_refused_when_it_has_to_run(tmp_path):
    """Without that flag nothing changes: a file meant to be run is still
    refused rather than written wrong."""
    from abmptools.gro2udf.guard import UnsupportedTopFeatureError
    from abmptools.gro2udf.top_exporter import TopExporter

    top = tmp_path / "v.top"
    top.write_text(VIEW_TOP, encoding="utf-8")
    gro = tmp_path / "v.gro"
    gro.write_text("v\n    3\n"
                   "    1REB     A1    1   2.000   2.000   2.000\n"
                   "    1REB     B1    2   2.470   2.000   2.000\n"
                   "    1REB     C1    3   2.800   2.340   2.000\n"
                   "   6.00000   6.00000   6.00000\n", encoding="utf-8")
    tpl = pathlib.Path(__file__).resolve().parents[1] / \
        "abmptools" / "gro2udf" / "default_template.udf"
    with pytest.raises(UnsupportedTopFeatureError):
        TopExporter().export(str(top), str(gro), str(tpl),
                             str(tmp_path / "r.udf"))


# --- where the exact route ends ----------------------------------------------

BILAYER_TOP = """\
[ defaults ]
1 2

[ atomtypes ]
  Q5   72.0 0.000 A 0.0 0.0
  C1   72.0 0.000 A 0.0 0.0

[ nonbond_params ]
  Q5  Q5  1 0.47 3.50
  Q5  C1  1 0.47 2.00
  C1  C1  1 0.47 3.50

[moleculetype]
  LIP  1

[atoms]
  1 Q5 1 LIP H1 1 -1.0 72.0
  2 C1 1 LIP T1 2  0.0 72.0
  3 C1 1 LIP T2 3  0.0 72.0

[bonds]
  1 2 1 0.47 1250

[constraints]
  2 3 1 0.47

[angles]
  1 2 3 2 180.0 35.0

[ system ]
bilayer

[ molecules ]
LIP 2
"""


def test_a_lipid_carries_over_with_nothing_left_out(tmp_path):
    """The shape a bilayer actually has -- a bond, a constraint, a GROMOS-96
    angle and per-pair LJ -- goes across whole.

    This is the line to hold: a lipid-and-water box needs no viewing mode,
    and the real one was checked the same way (a POPC bilayer with water and
    ions: 4 molecule types, 5692 molecules, nothing left out). What breaks it
    is a peptide or a protein in the same box, whose fitted dihedrals have
    phases no torsion here can hold -- not the bilayer.
    """
    top = tmp_path / "bi.top"
    top.write_text(BILAYER_TOP, encoding="utf-8")
    raw = TopParser(constraints_as_bonds=50000.0).parse(str(top))
    fatal, _ = check_top(raw, raw.sections)
    assert fatal == [], fatal


def test_the_same_box_with_a_fitted_dihedral_is_not_exact(tmp_path):
    """And the cause is named, so the boundary cannot quietly move."""
    top = tmp_path / "bi2.top"
    top.write_text(BILAYER_TOP.replace("""[ system ]""", """[dihedrals]
  1 2 3 1 1 -33.3 75.0 1

[ system ]"""), encoding="utf-8")
    raw = TopParser(constraints_as_bonds=50000.0).parse(str(top))
    fatal, _ = check_top(raw, raw.sections)
    assert any("phase" in f for f in fatal), fatal
