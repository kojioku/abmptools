# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.models."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from abmptools.formulation.models import (
    BileSaltSpec,
    EnhancerSpec,
    EquilibrationProtocol,
    FormulationBuildConfig,
    PeptideSpec,
    ProductionProtocol,
    SystemSpec,
    USProtocol,
)


# ---------------------------------------------------------------------------
# PeptideSpec
# ---------------------------------------------------------------------------


def test_peptide_spec_from_sequence():
    p = PeptideSpec(name="ala5", sequence="AAAAA", n_copies=2)
    assert p.sequence == "AAAAA"
    assert p.n_copies == 2
    assert p.pdb_path == ""


def test_peptide_spec_from_pdb():
    p = PeptideSpec(name="insulin", pdb_path="insulin.pdb",
                    disulfide_bonds=[("A:6:SG", "A:11:SG")])
    assert p.pdb_path == "insulin.pdb"
    assert p.disulfide_bonds == [("A:6:SG", "A:11:SG")]


def test_peptide_spec_requires_one_of_sequence_or_pdb():
    with pytest.raises(ValueError, match="exactly one"):
        PeptideSpec(name="x")


def test_peptide_spec_rejects_both_sequence_and_pdb():
    with pytest.raises(ValueError, match="exactly one"):
        PeptideSpec(name="x", sequence="AAA", pdb_path="x.pdb")


def test_peptide_spec_rejects_zero_copies():
    with pytest.raises(ValueError, match="n_copies must be >= 1"):
        PeptideSpec(name="x", sequence="A", n_copies=0)


def test_peptide_spec_disulfide_pair_format():
    with pytest.raises(ValueError, match="2-tuples"):
        PeptideSpec(name="x", sequence="CC", disulfide_bonds=[("A:1:SG",)])


# ---------------------------------------------------------------------------
# EnhancerSpec
# ---------------------------------------------------------------------------


def test_enhancer_spec_neutral_only():
    e = EnhancerSpec(name="caprate", resname="CPR",
                     smiles_neutral="CCCCCCCCCC(=O)O", n_neutral=16)
    assert e.n_neutral == 16 and e.n_charged == 0


def test_enhancer_spec_charged_only():
    e = EnhancerSpec(name="caprate", resname="CPC",
                     smiles_charged="CCCCCCCCCC(=O)[O-]", n_charged=16)
    assert e.n_charged == 16


def test_enhancer_spec_both_forms_50_50():
    e = EnhancerSpec(
        name="caprate", resname="CPR",
        smiles_neutral="CCCCCCCCCC(=O)O",
        smiles_charged="CCCCCCCCCC(=O)[O-]",
        n_neutral=16, n_charged=16,
    )
    assert e.n_neutral + e.n_charged == 32


def test_enhancer_spec_resname_required():
    with pytest.raises(ValueError, match="resname required"):
        EnhancerSpec(name="x", smiles_neutral="C", n_neutral=1)


def test_enhancer_spec_resname_too_long():
    with pytest.raises(ValueError, match="resname must be <= 4"):
        EnhancerSpec(name="x", resname="LONG5", smiles_neutral="C", n_neutral=1)


def test_enhancer_spec_zero_copies_rejected():
    with pytest.raises(ValueError, match="must be > 0"):
        EnhancerSpec(name="x", resname="XXX")


def test_enhancer_spec_missing_neutral_smiles():
    with pytest.raises(ValueError, match="no\n? *smiles_neutral"):
        EnhancerSpec(name="x", resname="XXX", n_neutral=5)


# ---------------------------------------------------------------------------
# BileSaltSpec
# ---------------------------------------------------------------------------


def test_bile_salt_spec_basic():
    b = BileSaltSpec(name="taurocholate", resname="TCH",
                     smiles="O", n_copies=2, net_charge=-1)
    assert b.net_charge == -1


def test_bile_salt_spec_zero_copies_allowed():
    """n_copies=0 means 'skip this species' — allowed but trivial."""
    b = BileSaltSpec(name="tch", resname="TCH", n_copies=0)
    assert b.n_copies == 0


def test_bile_salt_spec_resname_required():
    with pytest.raises(ValueError, match="resname required"):
        BileSaltSpec(name="x", smiles="O", n_copies=2)


# ---------------------------------------------------------------------------
# SystemSpec
# ---------------------------------------------------------------------------


def _make_minimal_system() -> SystemSpec:
    return SystemSpec(
        peptides=[PeptideSpec(name="p", sequence="AAAAA", n_copies=2)],
        enhancers=[EnhancerSpec(name="caprate", resname="CPR",
                                smiles_neutral="C", n_neutral=4)],
        bile_salts=[BileSaltSpec(name="tch", resname="TCH",
                                 smiles="O", n_copies=2)],
    )


def test_system_spec_minimal():
    s = _make_minimal_system()
    assert s.total_peptide_copies == 2
    assert s.box_size_nm == 10.0


def test_system_spec_requires_peptide():
    with pytest.raises(ValueError, match="at least one PeptideSpec"):
        SystemSpec(peptides=[])


def test_system_spec_rejects_negative_box():
    with pytest.raises(ValueError, match="box_size_nm"):
        SystemSpec(
            peptides=[PeptideSpec(name="p", sequence="A")],
            box_size_nm=-1.0,
        )


def test_system_spec_rejects_resname_collision():
    with pytest.raises(ValueError, match="resnames must be unique"):
        SystemSpec(
            peptides=[PeptideSpec(name="p", sequence="A")],
            enhancers=[EnhancerSpec(name="x", resname="XXX",
                                    smiles_neutral="C", n_neutral=1)],
            bile_salts=[BileSaltSpec(name="y", resname="XXX",
                                     smiles="O", n_copies=1)],
        )


# ---------------------------------------------------------------------------
# JSON round-trip
# ---------------------------------------------------------------------------


def test_build_config_json_round_trip(tmp_path: Path):
    cfg = FormulationBuildConfig(
        system=_make_minimal_system(),
        equilibration=EquilibrationProtocol(em_steps=5000),
        production=ProductionProtocol(nsteps=10000),
        release_us=USProtocol(n_windows=8),
        output_dir=str(tmp_path),
    )
    path = tmp_path / "config.json"
    cfg.to_json(str(path))
    back = FormulationBuildConfig.from_json(str(path))
    assert back.equilibration.em_steps == 5000
    assert back.production.nsteps == 10000
    assert back.release_us is not None
    assert back.release_us.n_windows == 8
    assert back.system.enhancers[0].resname == "CPR"
    assert back.system.peptides[0].n_copies == 2


def test_build_config_json_without_release_us(tmp_path: Path):
    cfg = FormulationBuildConfig(
        system=_make_minimal_system(),
        output_dir=str(tmp_path),
    )
    path = tmp_path / "cfg.json"
    cfg.to_json(str(path))
    back = FormulationBuildConfig.from_json(str(path))
    assert back.release_us is None


def test_build_config_default_force_fields():
    cfg = FormulationBuildConfig(system=_make_minimal_system())
    assert cfg.amber_protein_ff == "leaprc.protein.ff14SB"
    assert cfg.amber_water_ff == "leaprc.water.tip3p"
    assert cfg.gaff_version == "gaff2"
