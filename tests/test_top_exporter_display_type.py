# -*- coding: utf-8 -*-
"""Atom-type display names for OCTA viewer (`_display_type_name`).

OCTA viewer reads the leading character of a type name as the element
symbol, so OpenFF interchange's per-atom type names have to be rewritten
to `<element><index>` or CPK colouring breaks.

The trap these tests pin down: interchange names the type after the
*molecule*. It produced ``MOL0_4`` when the molecule carried no name and
produces ``methane_4`` once one is set (openff-interchange 0.4.2 reading
``Molecule.name``). A rule hard-coded to ``MOL\\d+`` silently stopped
converting anything the moment callers began naming their molecules —
silently, because an unconverted name is still a valid string and only
the rendering downstream goes wrong.
"""
from __future__ import annotations

import pytest

from abmptools.gro2udf.top_exporter import _display_type_name


class TestOpenFFTypeNames:
    def test_named_molecule_is_converted(self):
        # The regression: was passed through unchanged, so OCTA viewer read
        # "m" as Mg/Mn.
        assert _display_type_name("methane_0", "C", {"methane"}) == "C0"
        assert _display_type_name("methane_11", "H", {"methane"}) == "H11"

    def test_unnamed_molecule_still_converted(self):
        # interchange's older output, and what callers get with no name set.
        assert _display_type_name("MOL0_4", "C", {"MOL0"}) == "C4"

    def test_legacy_literal_without_mol_names(self):
        # Callers that cannot supply the topology keep the old behaviour.
        assert _display_type_name("MOL0_4", "C") == "C4"
        assert _display_type_name("MOL1_2", "N") == "N2"

    def test_underscore_name_is_handled(self):
        # A molecule name may itself contain underscores; only the final
        # `_<digits>` group is the atom index.
        assert _display_type_name("ethyl_acetate_7", "O",
                                  {"ethyl_acetate"}) == "O7"


class TestConventionalTypeNamesPassThrough:
    @pytest.mark.parametrize("type_name", ["c3", "ca", "hc", "n4"])
    def test_gaff(self, type_name):
        assert _display_type_name(type_name, "C", {"methane"}) == type_name

    def test_opls_is_not_mistaken_for_an_interchange_name(self):
        # `opls_267` has the same <prefix>_<digits> shape as an interchange
        # type. What separates them is that the prefix of an interchange
        # name is a moleculetype in the same topology; "opls" is not.
        assert _display_type_name("opls_267", "C", {"methane"}) == "opls_267"

    def test_prefix_not_in_topology_is_left_alone(self):
        assert _display_type_name("other_3", "C", {"methane"}) == "other_3"

    def test_name_without_index_is_left_alone(self):
        assert _display_type_name("methane", "C", {"methane"}) == "methane"

    def test_non_numeric_suffix_is_left_alone(self):
        assert _display_type_name("methane_x", "C", {"methane"}) == "methane_x"


class TestUniquenessIsPreserved:
    def test_per_atom_types_stay_distinct(self):
        # LJ Pair_Interaction needs one entry per atom type; collapsing two
        # atoms onto one display name would corrupt the UDF cross-references.
        names = {"methane"}
        got = [
            _display_type_name(f"methane_{i}", e, names)
            for i, e in enumerate(["C", "H", "H", "H", "H"])
        ]
        assert got == ["C0", "H1", "H2", "H3", "H4"]
        assert len(set(got)) == len(got)
