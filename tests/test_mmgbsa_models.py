# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.mmgbsa.models dataclasses."""
from __future__ import annotations

import pytest

from abmptools.genesis.mmgbsa.models import (
    DEFAULT_ADD_PDB_RES_MAP,
    EnergyProtocol,
    ForceFieldSet,
    LigandParameterization,
    MMGBSABuildConfig,
    MinimizationProtocol,
    TargetSpec,
)


# ---------------------------------------------------------------------------
# TargetSpec
# ---------------------------------------------------------------------------

class TestTargetSpec:
    def test_minimal_valid(self):
        spec = TargetSpec(pdb="3772L.pdb", ligand_resno=201)
        assert spec.pdb == "3772L.pdb"
        assert spec.ligand_resno == 201
        assert spec.chain is None
        assert spec.name is None

    def test_with_chain(self):
        spec = TargetSpec(pdb="x.pdb", ligand_resno=42, chain="A")
        assert spec.chain == "A"

    def test_empty_pdb_raises(self):
        with pytest.raises(ValueError, match="pdb is required"):
            TargetSpec(pdb="", ligand_resno=1)

    def test_zero_resno_raises(self):
        with pytest.raises(ValueError, match="ligand_resno"):
            TargetSpec(pdb="x.pdb", ligand_resno=0)

    def test_negative_resno_raises(self):
        with pytest.raises(ValueError, match="ligand_resno"):
            TargetSpec(pdb="x.pdb", ligand_resno=-1)

    def test_multi_char_chain_raises(self):
        with pytest.raises(ValueError, match="1-character"):
            TargetSpec(pdb="x.pdb", ligand_resno=1, chain="AB")


# ---------------------------------------------------------------------------
# ForceFieldSet
# ---------------------------------------------------------------------------

class TestForceFieldSet:
    def test_defaults_match_poc(self):
        ff = ForceFieldSet()
        assert ff.leaprc_protein == "leaprc.protein.ff14SB"
        assert ff.leaprc_dna == "leaprc.DNA.OL15"
        assert ff.leaprc_rna == "leaprc.RNA.OL3"
        assert ff.leaprc_water == "leaprc.water.tip3p"
        assert ff.leaprc_gaff2 == "leaprc.gaff2"
        assert ff.leaprc_gaff == "leaprc.gaff"
        assert ff.add_pdb_res_map == DEFAULT_ADD_PDB_RES_MAP
        assert ff.add_pdb_atom_map == DEFAULT_ADD_PDB_RES_MAP

    def test_custom_overrides(self):
        ff = ForceFieldSet(
            leaprc_protein="leaprc.protein.ff14iSB",
            extra_lines=["loadAmberParams custom.frcmod"],
        )
        assert ff.leaprc_protein == "leaprc.protein.ff14iSB"
        assert ff.extra_lines == ["loadAmberParams custom.frcmod"]

    def test_json_list_normalised_to_tuples(self):
        # When loaded from JSON, list entries become inner lists; the
        # __post_init__ should normalise them to tuples for hashing/etc.
        ff = ForceFieldSet(
            add_pdb_res_map=[["LI+", "LI"], ["NA+", "NA"]],
        )
        assert all(isinstance(p, tuple) for p in ff.add_pdb_res_map)

    def test_invalid_pair_length_raises(self):
        with pytest.raises(ValueError, match="2-tuples"):
            ForceFieldSet(add_pdb_res_map=[("LI+", "LI", "extra")])


# ---------------------------------------------------------------------------
# LigandParameterization
# ---------------------------------------------------------------------------

class TestLigandParameterization:
    def test_defaults_match_poc(self):
        lp = LigandParameterization()
        assert lp.charge_method == "bcc"
        assert lp.extra_keys == "maxcyc=0"
        assert lp.atom_type == "gaff2"
        assert lp.skip_if_cached is True

    def test_invalid_charge_method_raises(self):
        with pytest.raises(ValueError, match="charge_method"):
            LigandParameterization(charge_method="nonexistent")

    def test_invalid_atom_type_raises(self):
        with pytest.raises(ValueError, match="atom_type"):
            LigandParameterization(atom_type="opls")

    def test_user_charge_method(self):
        lp = LigandParameterization(charge_method="user", net_charge=-1)
        assert lp.charge_method == "user"
        assert lp.net_charge == -1


# ---------------------------------------------------------------------------
# EnergyProtocol
# ---------------------------------------------------------------------------

class TestEnergyProtocol:
    def test_defaults_match_poc(self):
        e = EnergyProtocol()
        assert e.electrostatic == "CUTOFF"
        assert e.cutoffdist_A == 99.9
        assert e.pairlistdist_A == 100.0
        assert e.implicit_solvent == "GBSA"
        assert e.gbsa_salt_cons == 0.15
        assert e.gbsa_surf_tens == 0.0072
        assert e.gbsa_vdw_offset == 0.09

    def test_invalid_implicit_solvent_raises(self):
        with pytest.raises(ValueError, match="implicit_solvent"):
            EnergyProtocol(implicit_solvent="POISSON")

    def test_invalid_electrostatic_raises(self):
        with pytest.raises(ValueError, match="electrostatic"):
            EnergyProtocol(electrostatic="EWALD")

    def test_zero_cutoff_raises(self):
        with pytest.raises(ValueError, match="cutoffdist_A"):
            EnergyProtocol(cutoffdist_A=0.0)

    def test_pairlist_lt_cutoff_raises(self):
        with pytest.raises(ValueError, match="pairlistdist_A"):
            EnergyProtocol(cutoffdist_A=12.0, pairlistdist_A=8.0)

    def test_implicit_none_allowed(self):
        e = EnergyProtocol(implicit_solvent="NONE")
        assert e.implicit_solvent == "NONE"


# ---------------------------------------------------------------------------
# MinimizationProtocol
# ---------------------------------------------------------------------------

class TestMinimizationProtocol:
    def test_defaults_match_poc_singlepoint(self):
        m = MinimizationProtocol()
        assert m.method == "SD"
        assert m.nsteps == 1  # POC: 1-step single-point
        assert m.eneout_period == 1
        assert m.crdout_period == 1
        assert m.rstout_period == 1
        assert m.nbupdate_period == 1

    def test_lbfgs_allowed(self):
        m = MinimizationProtocol(method="LBFGS", nsteps=100)
        assert m.method == "LBFGS"
        assert m.nsteps == 100

    def test_invalid_method_raises(self):
        with pytest.raises(ValueError, match="method must be SD or LBFGS"):
            MinimizationProtocol(method="CG")

    def test_zero_nsteps_raises(self):
        with pytest.raises(ValueError, match="nsteps must be >= 1"):
            MinimizationProtocol(nsteps=0)

    def test_zero_period_raises(self):
        with pytest.raises(ValueError, match="eneout_period"):
            MinimizationProtocol(eneout_period=0)


# ---------------------------------------------------------------------------
# MMGBSABuildConfig
# ---------------------------------------------------------------------------

def _valid_targets() -> list:
    return [TargetSpec(pdb="3772L.pdb", ligand_resno=201)]


class TestMMGBSABuildConfig:
    def test_defaults_with_targets(self):
        cfg = MMGBSABuildConfig(targets=_valid_targets())
        assert cfg.project_name == "mmgbsa_run"
        assert cfg.force_field.leaprc_protein == "leaprc.protein.ff14SB"
        assert cfg.energy.implicit_solvent == "GBSA"
        assert cfg.minimize.nsteps == 1
        assert cfg.mpi_processes == 1
        assert cfg.fail_fast is True

    def test_folder_mode_only(self):
        cfg = MMGBSABuildConfig(input_dir="./input", targets=[])
        assert cfg.input_dir == "./input"
        assert cfg.targets == []

    def test_neither_targets_nor_input_dir_raises(self):
        with pytest.raises(ValueError, match="targets.*input_dir"):
            MMGBSABuildConfig(targets=[], input_dir="")

    def test_zero_mpi_raises(self):
        with pytest.raises(ValueError, match="mpi_processes"):
            MMGBSABuildConfig(targets=_valid_targets(), mpi_processes=0)

    def test_empty_project_name_raises(self):
        with pytest.raises(ValueError, match="project_name"):
            MMGBSABuildConfig(targets=_valid_targets(), project_name="")

    def test_to_json_roundtrip(self, tmp_path):
        cfg = MMGBSABuildConfig(
            targets=[
                TargetSpec(pdb="3772L.pdb", ligand_resno=201),
                TargetSpec(pdb="9MM52.pdb", ligand_resno=201, chain="A"),
            ],
            project_name="poc_reproduction",
        )
        out = tmp_path / "cfg.json"
        cfg.to_json(str(out))
        loaded = MMGBSABuildConfig.from_json(str(out))
        assert loaded.project_name == "poc_reproduction"
        assert len(loaded.targets) == 2
        assert loaded.targets[0].pdb == "3772L.pdb"
        assert loaded.targets[0].ligand_resno == 201
        assert loaded.targets[1].chain == "A"
        assert loaded.force_field.leaprc_protein == cfg.force_field.leaprc_protein
        assert loaded.energy.gbsa_salt_cons == cfg.energy.gbsa_salt_cons
        assert loaded.minimize.nsteps == cfg.minimize.nsteps

    def test_to_json_roundtrip_folder_mode(self, tmp_path):
        cfg = MMGBSABuildConfig(input_dir="./input", project_name="folder_mode")
        out = tmp_path / "cfg.json"
        cfg.to_json(str(out))
        loaded = MMGBSABuildConfig.from_json(str(out))
        assert loaded.input_dir == "./input"
        assert loaded.targets == []

    def test_pdb_res_map_roundtrip_normalises_tuples(self, tmp_path):
        # JSON stores tuples as lists; from_json must normalise back.
        cfg = MMGBSABuildConfig(targets=_valid_targets())
        out = tmp_path / "cfg.json"
        cfg.to_json(str(out))
        loaded = MMGBSABuildConfig.from_json(str(out))
        assert all(isinstance(p, tuple) for p in loaded.force_field.add_pdb_res_map)
