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
