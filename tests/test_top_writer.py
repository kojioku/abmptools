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
    _dedup_angles,
    _dedup_dihedrals,
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


# ---------------------------------------------------------------------------
# Angle / dihedral dedup
# ---------------------------------------------------------------------------

class TestDedupAngles:
    def test_drops_self_angle(self):
        """ai == ak は物理的に不可能なので落とす (acetone の '3 2 3' 等)。"""
        angles = [
            AngleRecord(atom1=3, atom2=2, atom3=3, theta0=120.0, k=400.0),
            AngleRecord(atom1=1, atom2=2, atom3=3, theta0=120.0, k=400.0),
        ]
        out = _dedup_angles(angles)
        assert len(out) == 1
        assert (out[0].atom1, out[0].atom2, out[0].atom3) == (1, 2, 3)

    def test_drops_duplicate_with_same_endpoints(self):
        """vertex が同じで end-pair が同じなら重複。"""
        angles = [
            AngleRecord(atom1=1, atom2=2, atom3=3, theta0=120.0, k=400.0),
            AngleRecord(atom1=1, atom2=2, atom3=3, theta0=120.0, k=400.0),
        ]
        assert len(_dedup_angles(angles)) == 1

    def test_treats_reversed_endpoints_as_duplicate(self):
        """1-2-3 と 3-2-1 は同じ angle。"""
        angles = [
            AngleRecord(atom1=1, atom2=2, atom3=3, theta0=120.0, k=400.0),
            AngleRecord(atom1=3, atom2=2, atom3=1, theta0=120.0, k=400.0),
        ]
        assert len(_dedup_angles(angles)) == 1

    def test_keeps_different_vertices(self):
        """vertex が違うなら別 angle。"""
        angles = [
            AngleRecord(atom1=1, atom2=2, atom3=3, theta0=120.0, k=400.0),
            AngleRecord(atom1=1, atom2=3, atom3=2, theta0=120.0, k=400.0),
        ]
        assert len(_dedup_angles(angles)) == 2


class TestDedupDihedrals:
    def test_drops_self_dihedral(self):
        """4 atom 中に重複があれば落とす。"""
        dihs = [
            DihedralRecord(atom1=1, atom2=2, atom3=3, atom4=2,  # atom2==atom4
                           funct="1", params=[0.0, 1.0, 1.0]),
            DihedralRecord(atom1=1, atom2=2, atom3=3, atom4=4,
                           funct="1", params=[0.0, 1.0, 1.0]),
        ]
        out = _dedup_dihedrals(dihs)
        assert len(out) == 1
        assert (out[0].atom1, out[0].atom4) == (1, 4)

    def test_treats_reversed_as_duplicate(self):
        """1-2-3-4 と 4-3-2-1 は同じ dihedral。"""
        dihs = [
            DihedralRecord(atom1=1, atom2=2, atom3=3, atom4=4,
                           funct="1", params=[0.0, 1.0, 1.0]),
            DihedralRecord(atom1=4, atom2=3, atom3=2, atom4=1,
                           funct="1", params=[0.0, 1.0, 1.0]),
        ]
        assert len(_dedup_dihedrals(dihs)) == 1

    def test_keeps_distinct_dihedrals(self):
        dihs = [
            DihedralRecord(atom1=1, atom2=2, atom3=3, atom4=4,
                           funct="1", params=[0.0, 1.0, 1.0]),
            DihedralRecord(atom1=2, atom2=3, atom3=4, atom4=5,
                           funct="1", params=[0.0, 1.0, 1.0]),
        ]
        assert len(_dedup_dihedrals(dihs)) == 2


class TestWriteWithDuplicateAngles:
    """End-to-end: spurious self-angles in MoleculeTopology don't show in .top."""

    def test_acetone_like_angles_get_cleaned(self, tmp_path):
        """C=O 系で UDF にあるような '3 2 3' (self) と duplicate を入れて、
        書き出された .top に self/duplicate が残らないことを確認。"""
        atom_types = [AtomType(name="c", mass=12.011, sigma=0.34, epsilon=0.36)]
        atoms = [
            AtomRecord(index=i, type_name="c", gro_name=f"C{i:03d}", charge=0.0)
            for i in (1, 2, 3, 4)
        ]
        bonds = [
            BondRecord(atom1=1, atom2=2, funct="1", r0=0.15, kb=2e5),
            BondRecord(atom1=2, atom2=3, funct="1", r0=0.13, kb=4e5),
            BondRecord(atom1=2, atom2=4, funct="1", r0=0.15, kb=2e5),
        ]
        # spurious entries an upstream gen_udf bug produces:
        angles = [
            AngleRecord(atom1=3, atom2=2, atom3=3, theta0=120.0, k=400.0),  # self
            AngleRecord(atom1=1, atom2=2, atom3=3, theta0=120.0, k=400.0),
            AngleRecord(atom1=1, atom2=2, atom3=3, theta0=120.0, k=400.0),  # dup
            AngleRecord(atom1=4, atom2=2, atom3=3, theta0=120.0, k=400.0),
            AngleRecord(atom1=4, atom2=2, atom3=3, theta0=120.0, k=400.0),  # dup
            AngleRecord(atom1=1, atom2=2, atom3=4, theta0=109.5, k=400.0),
            AngleRecord(atom1=2, atom2=3, atom3=2, theta0=0.0, k=0.0),      # self
        ]
        topo = MoleculeTopology(
            udf_name="acetone",
            gro_name="M0000",
            nrexcl=3,
            atoms=atoms,
            bonds=bonds,
            pairs=[],
            angles=angles,
            dihedrals=[],
        )
        model = SystemModel(
            title="dedup test",
            udf_path="/tmp/test.udf",
            comb_rule=2,
            flags14=0,
            fudgeLJ=0.5,
            fudgeQQ=0.8333,
            calcQQ=0,
            atom_types=atom_types,
            mol_topologies=[topo],
            mol_sequence=[("M0000", 1)],
            cell=CellGeometry(a=2.0, b=2.0, c=2.0),
        )

        filepath = str(tmp_path / "out.top")
        TopWriter().write(model, filepath)

        # extract [ angles ] section
        text = open(filepath).read()
        ang_start = text.index("[ angles ]")
        ang_end = text.index("[ dihedrals ]")
        ang_block = text[ang_start:ang_end]

        # count non-comment, non-empty lines
        body_lines = [
            ln for ln in ang_block.splitlines()
            if ln.strip() and not ln.strip().startswith(";")
            and not ln.strip().startswith("[")
        ]
        # 7 input → 3 unique non-self: (1,2,3), (4,2,3), (1,2,4)
        assert len(body_lines) == 3, (
            f"expected 3 unique angles, got {len(body_lines)}: {body_lines}"
        )
