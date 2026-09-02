# -*- coding: utf-8 -*-
"""Display type names must stay unique across the whole topology.

An interchange per-atom type is named ``<moleculetype>_<index>``, and the
index counts atoms **within its molecule**, restarting at zero for every
component. ``_display_type_name`` rewrites it to ``<element><index>`` so
the OCTA viewer can read the element off the front — and that rewrite
throws away the only thing separating one component from another. In a
mixture, two types collide as soon as they share an element and a
position: aripiprazole's ``APZ_2`` and lactic acid's ``LAC_2`` are both
oxygens at index 2, and both became ``O2``.

The UDF then held two ``Atom_Type`` entries under one name with different
Lennard-Jones parameters, and everything that resolves a type by name
took the first: **the second component silently ran on the first's
parameters.** Measured on APZ + lactic acid, 69 types collapsed to 68
names and LJ-14 came back 0.9 % wrong after a round trip.

A single-component system never hits this, which is why it went
unnoticed.
"""
import pytest

from abmptools.gro2udf.top_exporter import (
    _display_type_name,
    build_display_type_map,
)
from abmptools.gro2udf.top_model import MolAtomSpec, MolSpec, TopModel


def atom(i, element, type_name):
    return MolAtomSpec(index_1based=i, atom_name="%s%d" % (element, i),
                       element=element, type_name=type_name, charge=0.0,
                       global_atom_id=i - 1)


def model(*mols):
    """A TopModel carrying only what the name map reads."""
    specs = [
        MolSpec(name=name,
                atoms=[atom(i + 1, el, tn) for i, (el, tn) in enumerate(atoms)])
        for name, atoms in mols
    ]
    return TopModel(
        comb_rule=2, fudge_lj=0.5, fudge_qq=0.8333,
        atom_type_specs=[], bond_type_specs=[], angle_type_specs=[],
        torsion_type_specs=[], mass_dict={},
        mol_type_names=[s.name for s in specs], mol_specs=specs,
        mol_instance_list=[s.name for s in specs], frames=[],
        n_atoms_total=sum(len(s.atoms) for s in specs),
    )


class TestSingleComponent:
    """The common case must come out exactly as it did before."""

    def test_interchange_names_become_element_plus_index(self):
        m = model(("MOL0", [("C", "MOL0_0"), ("H", "MOL0_1"),
                            ("O", "MOL0_2")]))
        assert build_display_type_map(m) == {
            "MOL0_0": "C0", "MOL0_1": "H1", "MOL0_2": "O2"}

    def test_a_named_molecule_works_the_same_way(self):
        """interchange names the type after the molecule once it has a name."""
        m = model(("methane", [("C", "methane_0"), ("H", "methane_1")]))
        assert build_display_type_map(m) == {
            "methane_0": "C0", "methane_1": "H1"}

    def test_conventional_force_field_types_pass_through(self):
        m = model(("MOL", [("C", "c3"), ("H", "hc"), ("C", "ca")]))
        assert build_display_type_map(m) == {
            "c3": "c3", "hc": "hc", "ca": "ca"}


class TestCollision:

    def test_two_components_do_not_share_a_name(self):
        """The regression: APZ_2 and LAC_2 are both O at index 2."""
        m = model(
            ("APZ", [("Cl", "APZ_0"), ("Cl", "APZ_1"), ("O", "APZ_2"),
                     ("O", "APZ_3")]),
            ("LAC", [("C", "LAC_0"), ("C", "LAC_1"), ("O", "LAC_2")]),
        )
        display = build_display_type_map(m)
        assert len(set(display.values())) == len(display)
        # the first component keeps the names it had
        assert display["APZ_2"] == "O2"
        assert display["APZ_3"] == "O3"
        # the collided one continues above the highest O already used
        assert display["LAC_2"] == "O4"

    def test_the_first_component_is_never_displaced(self):
        """Whoever appears first keeps its preferred name, so a
        single-component result does not change when a second component is
        added alongside it."""
        alone = build_display_type_map(
            model(("APZ", [("O", "APZ_0"), ("O", "APZ_1")])))
        together = build_display_type_map(model(
            ("APZ", [("O", "APZ_0"), ("O", "APZ_1")]),
            ("LAC", [("O", "LAC_0"), ("O", "LAC_1")]),
        ))
        for raw, name in alone.items():
            assert together[raw] == name

    def test_three_components_all_stay_distinct(self):
        m = model(
            ("A", [("O", "A_0")]),
            ("B", [("O", "B_0")]),
            ("C", [("O", "C_0")]),
        )
        display = build_display_type_map(m)
        assert sorted(display.values()) == ["O0", "O1", "O2"]

    def test_different_elements_at_the_same_index_do_not_collide(self):
        """Only same element *and* same index is a collision."""
        m = model(
            ("A", [("C", "A_0"), ("O", "A_1")]),
            ("B", [("N", "B_0"), ("S", "B_1")]),
        )
        display = build_display_type_map(m)
        assert display == {"A_0": "C0", "A_1": "O1",
                           "B_0": "N0", "B_1": "S1"}

    def test_a_generated_name_never_lands_on_a_force_field_type(self):
        """Pass-through names are the force field's own and are reserved
        before anything is generated."""
        m = model(
            ("A", [("O", "O1")]),                 # a literal FF type named O1
            ("B", [("O", "B_0"), ("O", "B_1")]),  # B_1 would prefer O1
        )
        display = build_display_type_map(m)
        assert display["O1"] == "O1"
        assert display["B_1"] != "O1"
        assert len(set(display.values())) == 3

    def test_the_map_is_deterministic(self):
        m = model(
            ("APZ", [("O", "APZ_2")]),
            ("LAC", [("O", "LAC_2")]),
        )
        assert build_display_type_map(m) == build_display_type_map(m)

    def test_a_collision_is_logged(self, caplog):
        """Silently renaming would be its own trap: the name no longer
        matches the interchange index the user sees in the .top."""
        m = model(
            ("APZ", [("O", "APZ_2")]),
            ("LAC", [("O", "LAC_2")]),
        )
        with caplog.at_level("WARNING"):
            build_display_type_map(m)
        messages = [r.getMessage() for r in caplog.records]
        assert any("LAC_2" in msg and "already taken" in msg
                   for msg in messages), messages


class TestDisplayTypeNameItself:
    """The per-name helper is unchanged; pin its contract."""

    def test_rewrites_when_the_prefix_is_a_moleculetype(self):
        assert _display_type_name("MOL0_4", "C", {"MOL0"}) == "C4"

    def test_leaves_an_unrelated_prefix_alone(self):
        assert _display_type_name("opls_267", "C", {"MOL0"}) == "opls_267"

    def test_falls_back_to_the_MOL_literal_without_mol_names(self):
        assert _display_type_name("MOL0_4", "C") == "C4"
