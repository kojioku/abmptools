# -*- coding: utf-8 -*-
"""Tests for abmptools.amorphous.models module."""
import json

import pytest

from abmptools.amorphous.models import BuildConfig, ComponentSpec


# ---------------------------------------------------------------------------
# ComponentSpec
# ---------------------------------------------------------------------------

class TestComponentSpec:
    """Tests for the ComponentSpec dataclass."""

    def test_valid_with_smiles(self):
        """ComponentSpec can be created with a SMILES string."""
        spec = ComponentSpec(name="water", smiles="O")
        assert spec.smiles == "O"
        assert spec.sdf_path == ""

    def test_valid_with_sdf_path(self):
        """ComponentSpec can be created with an SDF file path."""
        spec = ComponentSpec(name="mol", sdf_path="/tmp/mol.sdf")
        assert spec.sdf_path == "/tmp/mol.sdf"
        assert spec.smiles == ""

    def test_raises_when_both_empty(self):
        """Raises ValueError if neither smiles nor sdf_path is provided."""
        with pytest.raises(ValueError, match="smiles.*sdf_path"):
            ComponentSpec(name="empty")


# ---------------------------------------------------------------------------
# BuildConfig
# ---------------------------------------------------------------------------

class TestBuildConfig:
    """Tests for the BuildConfig dataclass."""

    def test_default_values(self):
        """BuildConfig has expected defaults."""
        cfg = BuildConfig()
        assert cfg.density_g_cm3 == pytest.approx(0.8)
        assert cfg.temperature == pytest.approx(300.0)
        assert cfg.T_high == pytest.approx(600.0)
        assert cfg.pressure == pytest.approx(1.0)
        assert cfg.forcefield == "openff_unconstrained-2.1.0.offxml"
        assert cfg.packmol_tolerance == pytest.approx(2.0)
        assert cfg.dt == pytest.approx(0.001)
        assert cfg.components == []

    def test_to_json_from_json_roundtrip(self, tmp_path):
        """Saving to JSON and loading back produces an equivalent object."""
        comp = ComponentSpec(name="ethanol", smiles="CCO", weight_fraction=1.0)
        original = BuildConfig(
            components=[comp],
            total_molecules=100,
            density_g_cm3=0.9,
            temperature=350.0,
            seed=42,
        )
        json_path = str(tmp_path / "config.json")
        original.to_json(json_path)
        loaded = BuildConfig.from_json(json_path)

        assert loaded.total_molecules == original.total_molecules
        assert loaded.density_g_cm3 == pytest.approx(original.density_g_cm3)
        assert loaded.temperature == pytest.approx(original.temperature)
        assert loaded.seed == original.seed

    def test_from_json_preserves_component_specs(self, tmp_path):
        """ComponentSpec list survives the JSON roundtrip."""
        comps = [
            ComponentSpec(name="A", smiles="C"),
            ComponentSpec(name="B", sdf_path="/some/path.sdf", n_mol=10),
        ]
        cfg = BuildConfig(components=comps)
        json_path = str(tmp_path / "cfg.json")
        cfg.to_json(json_path)
        loaded = BuildConfig.from_json(json_path)

        assert len(loaded.components) == 2
        assert loaded.components[0].name == "A"
        assert loaded.components[0].smiles == "C"
        assert loaded.components[1].name == "B"
        assert loaded.components[1].sdf_path == "/some/path.sdf"
        assert loaded.components[1].n_mol == 10

    def test_json_file_is_valid_json(self, tmp_path):
        """The written file is parseable JSON."""
        cfg = BuildConfig()
        json_path = str(tmp_path / "out.json")
        cfg.to_json(json_path)
        with open(json_path) as f:
            data = json.load(f)
        assert isinstance(data, dict)
        assert "density_g_cm3" in data

    def test_multiple_mdp_override_defaults(self):
        """MDP override fields have expected default values."""
        cfg = BuildConfig()
        assert cfg.em_steps == 50000
        assert cfg.em_tol == pytest.approx(1000.0)
        assert cfg.nvt_high_nsteps == 100000
        assert cfg.npt_high_nsteps == 200000
        assert cfg.anneal_nsteps == 500000
        assert cfg.npt_low_nsteps == 500000
