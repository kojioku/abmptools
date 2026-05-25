# -*- coding: utf-8 -*-
"""Tests for abmptools.gro2udf.top_model."""
import pytest

from abmptools.gro2udf.top_model import (
    mass_to_element,
    KB_AMU_A2_PS2_K,
    TopModel,
    MolSpec,
    MolAtomSpec,
    AtomTypeSpec,
    BondTypeSpec,
    AngleTypeSpec,
    TorsionTypeSpec,
)


# ---------------------------------------------------------------------------
# mass_to_element tests
# ---------------------------------------------------------------------------

class TestMassToElement:

    @pytest.mark.parametrize("mass, expected", [
        (1.0, "H"),
        (12.0, "C"),
        (14.0, "N"),
        (16.0, "O"),
        (32.0, "S"),
    ])
    def test_exact_integer_masses(self, mass, expected):
        assert mass_to_element(mass) == expected

    def test_fractional_mass_rounds(self):
        """12.011 rounds to 12 -> 'C'."""
        assert mass_to_element(12.011) == "C"

    def test_unknown_mass(self):
        assert mass_to_element(999.0) == "X"

    @pytest.mark.parametrize("mass, expected", [
        (1.008, "H"),       # GAFF h1/hc/ho/hn 等
        (4.0026, "He"),     # He-4 (rare)
        (6.94, "Li"),       # Li (rounds to 7)
        (9.012, "Be"),      # Be-9
        (10.81, "B"),       # B (rounds to 11)
        (12.011, "C"),      # GAFF c/c3/ca 等
        (14.007, "N"),      # GAFF n/n3/n4 等
        (15.999, "O"),      # GAFF o/oh/os 等
        (18.998, "F"),      # GAFF f
        (32.065, "S"),      # GAFF s/sh 等
        (35.453, "Cl"),     # GAFF cl
        (79.904, "Br"),     # GAFF br
        (126.904, "I"),     # GAFF i
    ])
    def test_force_field_typical_masses(self, mass, expected):
        """Common GROMACS / AMBER FF masses (the ones GAFF / GAFF2 / OpenFF
        actually emit) must map to the correct element."""
        assert mass_to_element(mass) == expected

    @pytest.mark.parametrize("mass", [
        3.0,    # was incorrectly mapped to 'Li' in pre-fix table
        5.0,    # was incorrectly mapped to 'B'
        6.0,    # was incorrectly mapped to 'Li' (overlapped 7 amu entry)
        8.0,    # no element at this mass
    ])
    def test_unphysical_masses_yield_X(self, mass):
        """3/5/6/8 amu have no naturally-occurring element. A bogus mass
        from a broken topology should yield the explicit fallback "X",
        not silently mislabel atoms as Li/B (pre-2026-05-25 bug)."""
        assert mass_to_element(mass) == "X"


# ---------------------------------------------------------------------------
# Constant test
# ---------------------------------------------------------------------------

class TestConstants:

    def test_kb_value(self):
        assert KB_AMU_A2_PS2_K == pytest.approx(0.83144626)


# ---------------------------------------------------------------------------
# Dataclass creation tests
# ---------------------------------------------------------------------------

class TestDataclassCreation:

    def test_mol_spec(self):
        mol = MolSpec(name="MOL")
        assert mol.name == "MOL"
        assert mol.atoms == []
        assert mol.bonds == []

    def test_mol_atom_spec(self):
        atom = MolAtomSpec(
            index_1based=1,
            atom_name="ca1",
            element="C",
            type_name="ca",
            charge=-0.1,
            global_atom_id=0,
        )
        assert atom.element == "C"
        assert atom.charge == pytest.approx(-0.1)

    def test_top_model_creation(self):
        model = TopModel(
            comb_rule=2,
            fudge_lj=0.5,
            fudge_qq=0.8333,
            atom_type_specs=[],
            bond_type_specs=[],
            angle_type_specs=[],
            torsion_type_specs=[],
            mass_dict={},
            mol_type_names=[],
            mol_specs=[],
            mol_instance_list=[],
        )
        assert model.comb_rule == 2
        assert model.n_atoms_total == 0
        assert model.frames == []
