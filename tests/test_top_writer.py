# -*- coding: utf-8 -*-
"""Tests for abmptools.udf2gro.gromacs.writers.top_writer."""
import pytest

from abmptools.core.system_model import (
    AngleRecord,
    AtomRecord,
    AtomType,
    BondRecord,
    CellGeometry,
    DihedralRecord,
    MoleculeTopology,
    PairRecord,
    SystemModel,
)
from abmptools.udf2gro.gromacs.writers.top_writer import (
    TopWriter,
    _f2s,
    _strl,
    _strr,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _minimal_model(include_topology: bool = True) -> SystemModel:
    """Build a SystemModel with one molecule topology for testing."""
    atom_types = [
        AtomType(name="ca", mass=12.011, sigma=0.339967, epsilon=0.359824),
        AtomType(name="ha", mass=1.008,  sigma=0.259964, epsilon=0.062760),
    ]

    mol_topos = []
    if include_topology:
        atoms = [
            AtomRecord(index=1, type_name="ca", gro_name="C001", charge=-0.130),
            AtomRecord(index=2, type_name="ha", gro_name="H002", charge=0.130),
        ]
        bonds = [
            BondRecord(atom1=1, atom2=2, funct="1", r0=0.108, kb=307105.6),
        ]
        pairs = [
            PairRecord(atom1=1, atom2=2),
        ]
        angles = [
            AngleRecord(atom1=1, atom2=2, atom3=1, theta0=120.0, k=292.88),
        ]
        dihedrals = [
            DihedralRecord(
                atom1=1, atom2=2, atom3=1, atom4=2,
                funct="9", params=[180.0, 15.167, 2],
            ),
        ]
        topo = MoleculeTopology(
            udf_name="benzene",
            gro_name="M0000",
            nrexcl=3,
            atoms=atoms,
            bonds=bonds,
            pairs=pairs,
            angles=angles,
            dihedrals=dihedrals,
        )
        mol_topos.append(topo)

    return SystemModel(
        title="test topology",
        udf_path="/tmp/test.udf",
        comb_rule=2,
        flags14=0,
        fudgeLJ=0.5,
        fudgeQQ=0.8333,
        calcQQ=0,
        atom_types=atom_types,
        mol_topologies=mol_topos,
        mol_sequence=[("M0000", 10)] if include_topology else [],
        cell=CellGeometry(a=3.0, b=3.0, c=3.0),
    )


# ---------------------------------------------------------------------------
# Tests for helper functions
# ---------------------------------------------------------------------------

class TestHelpers:

    def test_strl_left_justifies_and_pads(self):
        result = _strl("", "abc", 8)
        # "abc" left-justified in 8 chars + " "
        assert result == "abc      "

    def test_strr_right_justifies_and_pads(self):
        result = _strr("", "abc", 8)
        # "abc" right-justified in 8 chars + " "
        assert result == "     abc "

    def test_f2s_truncates_float(self):
        result = _f2s(3.14159265, 6)
        assert len(result) == 6
        assert result.startswith("3.1415")


# ---------------------------------------------------------------------------
# Tests for TopWriter
# ---------------------------------------------------------------------------

class TestTopWriter:

    def test_write_contains_defaults_section(self, tmp_path):
        model = _minimal_model()
        filepath = str(tmp_path / "out.top")
        TopWriter().write(model, filepath)

        content = open(filepath).read()
        assert "[ defaults ]" in content

    def test_write_contains_atomtypes_section(self, tmp_path):
        model = _minimal_model()
        filepath = str(tmp_path / "out.top")
        TopWriter().write(model, filepath)

        content = open(filepath).read()
        assert "[ atomtypes ]" in content

    def test_write_contains_moleculetype_per_topology(self, tmp_path):
        model = _minimal_model()
        filepath = str(tmp_path / "out.top")
        TopWriter().write(model, filepath)

        content = open(filepath).read()
        assert "[ moleculetype ]" in content
        assert "M0000" in content

    def test_write_contains_molecules_section(self, tmp_path):
        model = _minimal_model()
        filepath = str(tmp_path / "out.top")
        TopWriter().write(model, filepath)

        content = open(filepath).read()
        assert "[ molecules ]" in content

    def test_full_write_with_all_sections(self, tmp_path):
        """Verify a full write with one topology includes all expected sections."""
        model = _minimal_model(include_topology=True)
        filepath = str(tmp_path / "out.top")
        TopWriter().write(model, filepath)

        content = open(filepath).read()
        for section in [
            "[ defaults ]",
            "[ atomtypes ]",
            "[ moleculetype ]",
            "[ atoms ]",
            "[ bonds ]",
            "[ pairs ]",
            "[ angles ]",
            "[ dihedrals ]",
            "[ system ]",
            "[ molecules ]",
        ]:
            assert section in content, f"Missing section: {section}"
