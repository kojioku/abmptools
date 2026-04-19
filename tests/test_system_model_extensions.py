# -*- coding: utf-8 -*-
"""
Tests for the SystemModel extensions added in the
feature/system-model-for-amorphous branch:

- ClusterData / FixedLabel dataclasses
- SystemModel.cluster_data / fixed_labels / ensemble_family defaults
- classify_ensemble() and COGNAC_ONLY_ALGOS
- Writer-side rejection of ensemble_family == "cognac_only"
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from abmptools.core.system_model import (
    AtomPosition,
    CellGeometry,
    ClusterData,
    FixedLabel,
    SimulationParams,
    SystemModel,
    COGNAC_ONLY_ALGOS,
    classify_ensemble,
)


# ---------------------------------------------------------------------------
# classify_ensemble / COGNAC_ONLY_ALGOS
# ---------------------------------------------------------------------------

def test_classify_cognac_algorithms():
    for algo in ("NPT_Andersen_Kremer_Grest", "NPT_Andersen_Nose_Hoover"):
        assert classify_ensemble(algo) == "cognac_only"


def test_classify_gromacs_native_defaults_ok():
    for algo in ("md", "md-vv", "sd", "", "NVT_Nose_Hoover", "NVT_Berendsen"):
        assert classify_ensemble(algo) == "gromacs_ok"


def test_cognac_only_algos_is_frozenset():
    assert isinstance(COGNAC_ONLY_ALGOS, frozenset)
    # Core membership claims
    assert "NPT_Andersen_Kremer_Grest" in COGNAC_ONLY_ALGOS


# ---------------------------------------------------------------------------
# ClusterData / FixedLabel
# ---------------------------------------------------------------------------

def test_cluster_data_defaults():
    c = ClusterData()
    assert c.xyz == []
    assert c.n_per_cluster == 1
    assert c.cluster_file == ""


def test_fixed_label_defaults():
    f = FixedLabel()
    assert f.atom_indices == []
    assert f.label == "fixed"


def test_cluster_data_custom():
    c = ClusterData(
        xyz=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
        n_per_cluster=4,
        cluster_file="water_tetramer.xyz",
    )
    assert len(c.xyz) == 2
    assert c.n_per_cluster == 4
    assert c.cluster_file == "water_tetramer.xyz"


# ---------------------------------------------------------------------------
# SystemModel new-field defaults
# ---------------------------------------------------------------------------

def _minimal_system_model(**overrides) -> SystemModel:
    base = dict(
        title="t", udf_path="",
        comb_rule=2, flags14=0,
        fudgeLJ=0.5, fudgeQQ=0.5, calcQQ=1,
    )
    base.update(overrides)
    return SystemModel(**base)


def test_system_model_defaults_for_cognac_extensions():
    m = _minimal_system_model()
    assert m.cluster_data is None
    assert m.fixed_labels == []
    assert m.ensemble_family == "gromacs_ok"


def test_system_model_accepts_cluster_and_fixed_labels():
    m = _minimal_system_model()
    m.cluster_data = ClusterData(
        xyz=[(0.0, 0.0, 0.0)],
        n_per_cluster=4,
        cluster_file="water.xyz",
    )
    m.fixed_labels = [FixedLabel(atom_indices=[0, 1], label="fixed")]
    assert m.cluster_data.n_per_cluster == 4
    assert m.fixed_labels[0].label == "fixed"


# ---------------------------------------------------------------------------
# Writer validator: reject cognac_only
# ---------------------------------------------------------------------------

def _cognac_only_model() -> SystemModel:
    m = _minimal_system_model()
    m.ensemble_family = "cognac_only"
    m.sim_params = SimulationParams(
        title="t",
        algorithm="NPT_Andersen_Kremer_Grest",
        nsteps=1, dt=0.001,
        outputinterval=1, outputinterval2=1,
        integrator="md",
    )
    m.atom_positions = [
        AtomPosition(mol_id=1, mol_name_short="A", atom_gro_name="C1",
                     atom_id=1, x=0.0, y=0.0, z=0.0,
                     vx=0.0, vy=0.0, vz=0.0)
    ]
    m.cell = CellGeometry(a=1.0, b=1.0, c=1.0)
    return m


def test_gro_writer_rejects_cognac_only(tmp_path):
    from abmptools.udf2gro.gromacs.writers.gro_writer import GroWriter

    m = _cognac_only_model()
    with pytest.raises(ValueError, match="cognac_only"):
        GroWriter().write(m, str(tmp_path / "x.gro"))


def test_top_writer_rejects_cognac_only(tmp_path):
    from abmptools.udf2gro.gromacs.writers.top_writer import TopWriter

    m = _cognac_only_model()
    with pytest.raises(ValueError, match="cognac_only"):
        TopWriter().write(m, str(tmp_path / "x.top"))


def test_mdp_writer_rejects_cognac_only(tmp_path):
    from abmptools.udf2gro.gromacs.writers.mdp_writer import MdpWriter

    m = _cognac_only_model()
    with pytest.raises(ValueError, match="cognac_only"):
        MdpWriter().write(m, str(tmp_path / "x.mdp"))


def test_itp_writer_rejects_cognac_only(tmp_path):
    from abmptools.udf2gro.gromacs.writers.itp_writer import ItpWriter

    m = _cognac_only_model()
    with pytest.raises(ValueError, match="cognac_only"):
        ItpWriter().write(m, str(tmp_path / "x"))


def test_writers_include_algorithm_in_error_message():
    from abmptools.udf2gro.gromacs.writers.gro_writer import GroWriter

    m = _cognac_only_model()
    with pytest.raises(ValueError) as exc:
        GroWriter().write(m, "/tmp/whatever.gro")
    assert "NPT_Andersen_Kremer_Grest" in str(exc.value)
