# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.models module."""
import json

import pytest

from abmptools.cg.peptide.models import PeptideBuildConfig, PeptideSpec


# ---------------------------------------------------------------------------
# PeptideSpec
# ---------------------------------------------------------------------------

class TestPeptideSpec:
    """Tests for the PeptideSpec dataclass."""

    def test_valid_peptide(self):
        spec = PeptideSpec(name="kgg", sequence="KGG", count=10)
        assert spec.name == "kgg"
        assert spec.sequence == "KGG"
        assert spec.count == 10

    def test_lowercase_sequence_normalized_to_upper(self):
        spec = PeptideSpec(name="rgg", sequence="rgg", count=1)
        assert spec.sequence == "RGG"

    def test_default_count_is_one(self):
        spec = PeptideSpec(name="x", sequence="A")
        assert spec.count == 1

    def test_invalid_amino_acid_raises(self):
        with pytest.raises(ValueError, match="Invalid amino acid"):
            PeptideSpec(name="bad", sequence="KGZ")

    def test_invalid_amino_acid_lowercase_raises(self):
        # 大小区別なく Z は invalid
        with pytest.raises(ValueError, match="Invalid amino acid"):
            PeptideSpec(name="bad", sequence="kgz")

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError, match="sequence"):
            PeptideSpec(name="x", sequence="")

    def test_empty_name_raises(self):
        with pytest.raises(ValueError, match="name"):
            PeptideSpec(name="", sequence="K")

    def test_count_zero_raises(self):
        with pytest.raises(ValueError, match="count"):
            PeptideSpec(name="x", sequence="K", count=0)

    def test_count_negative_raises(self):
        with pytest.raises(ValueError, match="count"):
            PeptideSpec(name="x", sequence="K", count=-1)

    def test_all_20_standard_aa_accepted(self):
        spec = PeptideSpec(name="all", sequence="ACDEFGHIKLMNPQRSTVWY")
        assert len(spec.sequence) == 20


# ---------------------------------------------------------------------------
# PeptideBuildConfig
# ---------------------------------------------------------------------------

class TestPeptideBuildConfig:
    """Tests for the PeptideBuildConfig dataclass."""

    def test_default_values(self):
        cfg = PeptideBuildConfig()
        assert cfg.peptides == []
        assert cfg.box_size_nm == pytest.approx(0.0)
        assert cfg.box_lengths_nm is None
        assert cfg.solvent_enabled is True
        assert cfg.neutralize is True
        assert cfg.nacl_molar == pytest.approx(0.15)
        assert cfg.temperature == pytest.approx(310.0)
        assert cfg.martinize2_path == "martinize2"
        assert cfg.gmx_path == "gmx"
        assert cfg.tleap_path == "tleap"
        assert cfg.dt_fs == pytest.approx(20.0)
        assert cfg.martini_itp_dir == ""
        assert cfg.output_dir == "."
        assert cfg.seed is None

    def test_invalid_box_lengths_length(self):
        with pytest.raises(ValueError, match="3 values"):
            PeptideBuildConfig(box_lengths_nm=[10.0, 10.0])

    def test_invalid_box_lengths_negative(self):
        with pytest.raises(ValueError, match="positive"):
            PeptideBuildConfig(box_lengths_nm=[10.0, 10.0, -1.0])

    def test_negative_box_size_raises(self):
        with pytest.raises(ValueError, match="box_size_nm"):
            PeptideBuildConfig(box_size_nm=-5.0)

    def test_negative_nacl_raises(self):
        with pytest.raises(ValueError, match="nacl_molar"):
            PeptideBuildConfig(nacl_molar=-0.1)

    def test_zero_dt_raises(self):
        with pytest.raises(ValueError, match="dt_fs"):
            PeptideBuildConfig(dt_fs=0.0)

    def test_negative_em_steps_raises(self):
        with pytest.raises(ValueError, match="em_steps"):
            PeptideBuildConfig(em_steps=-1)

    def test_duplicate_peptide_names_raise(self):
        peps = [
            PeptideSpec(name="dup", sequence="K"),
            PeptideSpec(name="dup", sequence="G"),
        ]
        with pytest.raises(ValueError, match="unique"):
            PeptideBuildConfig(peptides=peps)

    def test_unique_peptide_names_ok(self):
        peps = [
            PeptideSpec(name="a", sequence="K"),
            PeptideSpec(name="b", sequence="G"),
        ]
        cfg = PeptideBuildConfig(peptides=peps)
        assert len(cfg.peptides) == 2

    def test_to_json_from_json_roundtrip(self, tmp_path):
        peps = [
            PeptideSpec(name="kgg", sequence="KGG", count=5),
            PeptideSpec(name="rgg", sequence="RGGRGG", count=3),
        ]
        original = PeptideBuildConfig(
            peptides=peps,
            box_size_nm=12.0,
            temperature=300.0,
            nacl_molar=0.20,
            seed=42,
        )
        json_path = str(tmp_path / "cfg.json")
        original.to_json(json_path)
        loaded = PeptideBuildConfig.from_json(json_path)

        assert len(loaded.peptides) == 2
        assert loaded.peptides[0].name == "kgg"
        assert loaded.peptides[0].sequence == "KGG"
        assert loaded.peptides[0].count == 5
        assert loaded.peptides[1].sequence == "RGGRGG"
        assert loaded.box_size_nm == pytest.approx(12.0)
        assert loaded.temperature == pytest.approx(300.0)
        assert loaded.nacl_molar == pytest.approx(0.20)
        assert loaded.seed == 42

    def test_json_file_is_valid_json(self, tmp_path):
        cfg = PeptideBuildConfig()
        json_path = str(tmp_path / "out.json")
        cfg.to_json(json_path)
        with open(json_path) as f:
            data = json.load(f)
        assert isinstance(data, dict)
        assert "box_size_nm" in data
        assert "peptides" in data
        assert data["peptides"] == []

    def test_mdp_default_steps(self):
        cfg = PeptideBuildConfig()
        assert cfg.em_steps == 50000
        assert cfg.nvt_nsteps == 250000
        assert cfg.npt_nsteps == 250000
        assert cfg.md_nsteps == 5000000

    def test_rectangular_box(self):
        cfg = PeptideBuildConfig(box_lengths_nm=[10.0, 12.0, 15.0])
        assert cfg.box_lengths_nm == [10.0, 12.0, 15.0]

    def test_box_lengths_round_trip(self, tmp_path):
        cfg = PeptideBuildConfig(box_lengths_nm=[10.0, 12.0, 15.0])
        json_path = str(tmp_path / "cfg.json")
        cfg.to_json(json_path)
        loaded = PeptideBuildConfig.from_json(json_path)
        assert loaded.box_lengths_nm == [10.0, 12.0, 15.0]
