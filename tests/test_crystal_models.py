# -*- coding: utf-8 -*-
"""Unit tests for ``abmptools.crystal.models``.

Round-trip tests (JSON + YAML) plus per-dataclass validation. No
external scientific stacks required — pyyaml is the only optional
dependency, and the YAML tests skip if it is unavailable.
"""
from __future__ import annotations

from dataclasses import asdict

import pytest

from abmptools.crystal.models import (
    CIFEngineConfig,
    CIFInputSpec,
    CrystalBuildConfig,
    FMOMethod,
    FragmentTemplate,
    HPCJobSpec,
    PostProcessSpec,
)


def _example_config() -> CrystalBuildConfig:
    return CrystalBuildConfig(
        inputs=[
            CIFInputSpec(
                cif="XXXI-MMFF-R00001.cif",
                layer=5,
                atoms_in_mol=[32],
            ),
            CIFInputSpec(
                cif="XXXI-MMFF-R00002.cif",
                layer=5,
                atoms_in_mol=[32, 32],
            ),
        ],
        cif_engine=CIFEngineConfig(engine="legacy"),
        fragment=FragmentTemplate(
            cutmode="around", solutes=[0], criteria=6.0,
            molname=["UNK"], pieda=True, cmm=True,
        ),
        fmo=FMOMethod(
            method="MP2", basis_set="6-31Gdag", memory=6000,
            is_xyz=True,
        ),
        hpc=HPCJobSpec(
            scheduler="PJM", queue="small", group="hp190133",
            nodes=12, proc_per_node=2, elapse="24:00:00",
            abinit_dir="/data/hp190133/programs/ABINIT-MP",
        ),
        postproc=PostProcessSpec(
            enable=True, frag_target="1-10",
            annotate_nearest_atoms=True,
        ),
        output_dir="./out",
        project_name="csp7_layer5_around6",
    )


# ---------------------------------------------------------------------------
# CIFInputSpec
# ---------------------------------------------------------------------------

def test_cif_input_spec_defaults():
    spec = CIFInputSpec(cif="x.cif")
    assert spec.layer == 5
    assert spec.atoms_in_mol == [32]
    assert spec.asymmetric_only is False
    assert spec.name is None


def test_cif_input_spec_requires_cif():
    with pytest.raises(ValueError, match="cif is required"):
        CIFInputSpec(cif="")


def test_cif_input_spec_layer_validation():
    with pytest.raises(ValueError, match="layer must be"):
        CIFInputSpec(cif="x.cif", layer=0)


def test_cif_input_spec_atoms_in_mol_validation():
    with pytest.raises(ValueError, match="atoms_in_mol must be non-empty"):
        CIFInputSpec(cif="x.cif", atoms_in_mol=[])
    with pytest.raises(ValueError, match="atoms_in_mol entries must be > 0"):
        CIFInputSpec(cif="x.cif", atoms_in_mol=[32, 0])


# ---------------------------------------------------------------------------
# CIFEngineConfig
# ---------------------------------------------------------------------------

def test_cif_engine_config_defaults():
    cfg = CIFEngineConfig()
    assert cfg.engine == "legacy"
    assert cfg.detect_z_prime is True
    assert cfg.use_neighbor_bonds is True
    assert cfg.bond_tolerance == pytest.approx(0.4)


def test_cif_engine_config_engine_validation():
    with pytest.raises(ValueError, match="engine must be"):
        CIFEngineConfig(engine="gemmi")


def test_cif_engine_config_bond_tolerance_validation():
    with pytest.raises(ValueError, match="bond_tolerance"):
        CIFEngineConfig(bond_tolerance=-0.1)


# ---------------------------------------------------------------------------
# FragmentTemplate
# ---------------------------------------------------------------------------

def test_fragment_template_defaults():
    ft = FragmentTemplate()
    assert ft.cutmode == "around"
    assert ft.solutes == [0]
    assert ft.criteria == pytest.approx(6.0)
    assert ft.molname == ["UNK"]


def test_fragment_template_cutmode_validation():
    with pytest.raises(ValueError, match="cutmode must be"):
        FragmentTemplate(cutmode="bogus")


def test_fragment_template_cube_criteria():
    # cube cutmode requires [x, y, z]
    with pytest.raises(ValueError, match="\\[x, y, z\\]"):
        FragmentTemplate(cutmode="cube", criteria=6.0)
    # ok: list of 3
    FragmentTemplate(cutmode="cube", criteria=[6.0, 6.0, 6.0])


def test_fragment_template_molname_required():
    with pytest.raises(ValueError, match="molname must be non-empty"):
        FragmentTemplate(molname=[])


# ---------------------------------------------------------------------------
# FMOMethod
# ---------------------------------------------------------------------------

def test_fmo_method_defaults():
    m = FMOMethod()
    assert m.method == "MP2"
    assert m.basis_set == "6-31Gdag"
    assert m.is_xyz is True  # crystal pipeline default = direct coords


def test_fmo_method_validations():
    with pytest.raises(ValueError, match="memory"):
        FMOMethod(memory=0)
    with pytest.raises(ValueError, match="npro"):
        FMOMethod(npro=0)


# ---------------------------------------------------------------------------
# HPCJobSpec
# ---------------------------------------------------------------------------

def test_hpc_jobspec_defaults():
    hpc = HPCJobSpec()
    assert hpc.scheduler == "PJM"
    assert hpc.queue == "small"


def test_hpc_jobspec_scheduler_validation():
    with pytest.raises(ValueError, match="scheduler must be"):
        HPCJobSpec(scheduler="LSF")


def test_hpc_jobspec_local():
    HPCJobSpec(scheduler="local")  # valid


def test_hpc_jobspec_count_validations():
    with pytest.raises(ValueError, match="nodes"):
        HPCJobSpec(nodes=0)
    with pytest.raises(ValueError, match="proc_per_node"):
        HPCJobSpec(proc_per_node=0)


# ---------------------------------------------------------------------------
# PostProcessSpec
# ---------------------------------------------------------------------------

def test_postproc_spec_defaults():
    pp = PostProcessSpec()
    assert pp.enable is False
    assert pp.nearest_atom_count == 3


def test_postproc_spec_validations():
    with pytest.raises(ValueError, match="nearest_atom_count"):
        PostProcessSpec(nearest_atom_count=0)


# ---------------------------------------------------------------------------
# CrystalBuildConfig
# ---------------------------------------------------------------------------

def test_crystal_build_config_requires_inputs():
    with pytest.raises(ValueError, match="inputs must list"):
        CrystalBuildConfig(inputs=[])


def test_crystal_build_config_requires_project_name():
    with pytest.raises(ValueError, match="project_name"):
        CrystalBuildConfig(
            inputs=[CIFInputSpec(cif="x.cif")],
            project_name="",
        )


def test_crystal_build_config_json_round_trip(tmp_path):
    cfg = _example_config()
    path = tmp_path / "crystal.json"
    cfg.to_json(str(path))
    loaded = CrystalBuildConfig.from_json(str(path))
    assert asdict(loaded) == asdict(cfg)


def test_crystal_build_config_yaml_round_trip(tmp_path):
    pytest.importorskip("yaml")
    cfg = _example_config()
    path = tmp_path / "crystal.yaml"
    cfg.to_yaml(str(path))
    loaded = CrystalBuildConfig.from_yaml(str(path))
    assert asdict(loaded) == asdict(cfg)


def test_crystal_build_config_yaml_empty_file_raises(tmp_path):
    pytest.importorskip("yaml")
    path = tmp_path / "empty.yaml"
    path.write_text("")
    with pytest.raises(ValueError, match="empty"):
        CrystalBuildConfig.from_yaml(str(path))


def test_crystal_build_config_round_trip_preserves_zprime_2(tmp_path):
    """Z'>=2 atoms_in_mol = [32, 32] survives the JSON round-trip."""
    cfg = _example_config()
    path = tmp_path / "crystal.json"
    cfg.to_json(str(path))
    loaded = CrystalBuildConfig.from_json(str(path))
    assert loaded.inputs[1].atoms_in_mol == [32, 32]
