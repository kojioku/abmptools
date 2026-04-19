# -*- coding: utf-8 -*-
"""
Tests for abmptools.amorphous.system_model_adapter.from_interchange.

Two layers:

1. A mock-based test that feeds a hand-built :class:`TopModel` into the
   private ``_topmodel_to_systemmodel`` projection, checking the field
   mapping deterministically without OpenFF or Packmol.
2. An integration test (@pytest.mark.slow) that runs a real Packmol +
   OpenFF + Interchange build for methane × 4 and verifies the
   returned SystemModel is usable by GroWriter.
"""
from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest

from abmptools.amorphous.system_model_adapter import _topmodel_to_systemmodel
from abmptools.gro2udf.top_model import (
    AtomTypeSpec,
    BondTypeSpec,
    AngleTypeSpec,
    TorsionTypeSpec,
    GROFrameData,
    MolAtomSpec,
    MolSpec,
    TopModel,
)


# ---------------------------------------------------------------------------
# Unit: _topmodel_to_systemmodel
# ---------------------------------------------------------------------------

def _make_fake_topmodel():
    atom_type_specs = [
        AtomTypeSpec(name="c3", mass=12.011, sigma=0.34, epsilon=0.46),
        AtomTypeSpec(name="hc", mass=1.008,  sigma=0.26, epsilon=0.07),
    ]

    mol = MolSpec(
        name="MOL0",
        atoms=[
            MolAtomSpec(index_1based=1, atom_name="C1", element="C",
                        type_name="c3", charge=-0.1, global_atom_id=0),
            MolAtomSpec(index_1based=2, atom_name="H1", element="H",
                        type_name="hc", charge=0.025, global_atom_id=1),
            MolAtomSpec(index_1based=3, atom_name="H2", element="H",
                        type_name="hc", charge=0.025, global_atom_id=2),
        ],
    )

    frame = GROFrameData(
        step=0, time=0.0,
        coord_list=[
            [0.10, 0.10, 0.10],
            [0.20, 0.10, 0.10],
            [0.10, 0.20, 0.10],
            [1.10, 1.10, 1.10],
            [1.20, 1.10, 1.10],
            [1.10, 1.20, 1.10],
        ],
        cell=[2.0, 2.0, 2.0],
    )

    tm = TopModel(
        comb_rule=2,
        fudge_lj=0.5,
        fudge_qq=0.5,
        atom_type_specs=atom_type_specs,
        bond_type_specs=[],
        angle_type_specs=[],
        torsion_type_specs=[],
        mass_dict={"c3": 12.011, "hc": 1.008},
        mol_type_names=["MOL0"],
        mol_specs=[mol],
        mol_instance_list=["MOL0", "MOL0"],  # two molecules
        frames=[frame],
        n_atoms_total=6,
        ref_t=300.0,
        tau_t=0.1,
        ewald_r_cutoff=10.0,
    )
    return tm


def test_projection_copies_force_field_defaults():
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="x")
    assert model.comb_rule == 2
    assert model.fudgeLJ == 0.5
    assert model.fudgeQQ == 0.5
    assert model.calcQQ == 1   # OpenFF path → always ON
    assert model.flags14 == 0


def test_projection_copies_atom_types():
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="x")
    names = [a.name for a in model.atom_types]
    assert names == ["c3", "hc"]
    assert model.atom_types[0].mass == pytest.approx(12.011)
    assert model.atom_types[1].epsilon == pytest.approx(0.07)


def test_projection_populates_positions_and_cell():
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="x")
    assert model.cell is not None
    assert model.cell.a == pytest.approx(2.0)
    # 2 molecules × 3 atoms = 6 positions
    assert len(model.atom_positions) == 6
    # Second molecule starts from coord index 3 onwards
    assert model.atom_positions[3].x == pytest.approx(1.10)


def test_projection_produces_collapsed_mol_sequence():
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="x")
    # Two consecutive MOL0 → collapsed to one tuple
    assert model.mol_sequence == [("MOL0", 2)]


def test_projection_leaves_mol_topologies_empty_by_design():
    """Phase 1 keeps bond/angle/dihedral empty; documented in module docstring."""
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="x")
    assert model.mol_topologies == []


def test_projection_defaults_to_gromacs_ok():
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="x")
    assert model.ensemble_family == "gromacs_ok"


def test_projection_title_propagates():
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="my_run_42")
    assert model.title == "my_run_42"


def test_projection_has_no_cluster_or_fixed_labels():
    tm = _make_fake_topmodel()
    model = _topmodel_to_systemmodel(tm, title="x")
    assert model.cluster_data is None
    assert model.fixed_labels == []


# ---------------------------------------------------------------------------
# Integration (slow): real Packmol + OpenFF + Interchange -> SystemModel
# ---------------------------------------------------------------------------

_openff = pytest.importorskip("openff.toolkit", reason="openff-toolkit not installed")
pytest.importorskip("openff.interchange", reason="openff-interchange not installed")
pytest.importorskip("rdkit", reason="rdkit not installed")
_packmol = shutil.which("packmol")
_sqm = shutil.which("sqm")

integration_mark = [
    pytest.mark.skipif(_packmol is None, reason="packmol not on PATH"),
    pytest.mark.skipif(_sqm is None, reason="sqm (AmberTools) not on PATH"),
    pytest.mark.slow,
]


@pytest.fixture(scope="module")
def _methane_system(tmp_path_factory):
    from abmptools.amorphous.molecule_prep import prepare_molecule, write_single_mol_pdb
    from abmptools.amorphous.packing import run_packmol
    from abmptools.amorphous.parameterizer import create_interchange
    from abmptools.amorphous.system_model_adapter import from_interchange

    work = tmp_path_factory.mktemp("ad_smoke")
    mol = prepare_molecule(smiles="C", name="methane")
    pdb = str(work / "methane.pdb")
    write_single_mol_pdb(mol, pdb)

    mixture = str(work / "mixture.pdb")
    run_packmol(pdb_paths=[pdb], counts=[4], box_size_nm=2.0,
                output_pdb=mixture, build_dir=str(work),
                tolerance=2.0, seed=1)

    inter = create_interchange(molecules=[mol], counts=[4],
                               box_size_nm=2.0, mixture_pdb=mixture)
    return from_interchange(inter, title="methane_integration")


@pytest.mark.parametrize("_p", [None], ids=[""])  # parametrize lets us tag slow
@pytest.mark.slow
def test_integration_methane_systemmodel_shape(_methane_system, _p):
    m = _methane_system
    assert m.title == "methane_integration"
    assert m.ensemble_family == "gromacs_ok"
    assert len(m.atom_types) == 5           # C + 4H unique per methane
    assert m.mol_sequence == [("MOL0", 4)]
    assert len(m.atom_positions) == 20      # 5 atoms × 4 methanes
    assert m.cell is not None
    assert m.cell.a == pytest.approx(2.0)


@pytest.mark.slow
def test_integration_methane_gro_round_trip(_methane_system, tmp_path):
    from abmptools.udf2gro.gromacs.writers.gro_writer import GroWriter
    out = tmp_path / "m.gro"
    GroWriter().write(_methane_system, str(out))
    lines = out.read_text().splitlines()
    assert lines[0] == "methane_integration"
    assert int(lines[1].strip()) == 20
    # box line matches
    box = [float(x) for x in lines[-1].split()[:3]]
    assert all(b == pytest.approx(2.0) for b in box)
