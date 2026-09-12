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

    def test_valid_with_pdb_path(self):
        """ComponentSpec accepts a PDB path (Phase 9-a oligomer)."""
        spec = ComponentSpec(name="oligo", pdb_path="/tmp/poly.pdb")
        assert spec.pdb_path == "/tmp/poly.pdb"
        assert spec.smiles == ""
        assert spec.sdf_path == ""

    def test_raises_when_all_three_empty(self):
        """Raises ValueError when none of smiles / sdf_path / pdb_path set."""
        with pytest.raises(ValueError, match="smiles.*sdf_path.*pdb_path"):
            ComponentSpec(name="empty")

    def test_raises_when_two_sources_set(self):
        """Setting both smiles and pdb_path is not allowed (must be exclusive)."""
        with pytest.raises(ValueError, match="exactly one"):
            ComponentSpec(name="ambig", smiles="O", pdb_path="/tmp/x.pdb")


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

    def test_forcefield_default_is_openff_organic(self):
        """forcefield default = openff_unconstrained-2.1.0.offxml (organic only)。"""
        cfg = BuildConfig()
        assert cfg.forcefield == "openff_unconstrained-2.1.0.offxml"

    def test_forcefield_accepts_list_for_water_override(self):
        """forcefield に list[str] を渡すと stacked FF として扱われる。

        典型的用途: water に TIP3P を後ろから上書き。
        """
        cfg = BuildConfig(
            forcefield=[
                "openff_unconstrained-2.1.0.offxml",
                "tip3p.offxml",
            ],
        )
        assert isinstance(cfg.forcefield, list)
        assert len(cfg.forcefield) == 2
        assert cfg.forcefield[0] == "openff_unconstrained-2.1.0.offxml"
        assert cfg.forcefield[1] == "tip3p.offxml"

    def test_forcefield_list_round_trips_via_json(self, tmp_path):
        """list[str] forcefield も to_json / from_json で round-trip する。"""
        cfg = BuildConfig(
            forcefield=[
                "openff_unconstrained-2.1.0.offxml",
                "tip3p.offxml",
            ],
        )
        json_path = str(tmp_path / "cfg.json")
        cfg.to_json(json_path)
        with open(json_path) as f:
            data = json.load(f)
        assert data["forcefield"] == [
            "openff_unconstrained-2.1.0.offxml",
            "tip3p.offxml",
        ]


# ---------------------------------------------------------------------------
# 熱浴・圧力浴が BuildConfig から .mdp まで届くか
#
# ここが繋がっていないと、 --thermostat / --barostat を指定しても既定の
# .mdp が黙って出る。 型も名前も合っているので気付けない類の抜け。
# ---------------------------------------------------------------------------

class TestThermostatBarostatPassthrough:

    def test_defaults_match_the_protocol_defaults(self):
        from abmptools.core.system_model import AnnealProtocol

        cfg = BuildConfig()
        proto = AnnealProtocol()
        assert cfg.thermostat == proto.thermostat == "V-rescale"
        assert cfg.barostat == proto.barostat == "C-rescale"

    def test_the_builder_hands_them_to_the_protocol(self):
        from abmptools.amorphous.builder import AmorphousBuilder

        cfg = BuildConfig(thermostat="Nose-Hoover",
                          barostat="Parrinello-Rahman")
        builder = AmorphousBuilder.__new__(AmorphousBuilder)
        builder.config = cfg
        proto = builder._make_protocol()

        assert proto.thermostat == "Nose-Hoover"
        assert proto.barostat == "Parrinello-Rahman"

    def test_the_choice_reaches_the_mdp_text(self):
        from abmptools.amorphous.builder import AmorphousBuilder
        from abmptools.amorphous.mdp_protocol import generate_npt_high_mdp

        builder = AmorphousBuilder.__new__(AmorphousBuilder)
        builder.config = BuildConfig(thermostat="Nose-Hoover",
                                     barostat="Parrinello-Rahman")
        text = generate_npt_high_mdp(builder._make_protocol())

        assert "Nose-Hoover" in text
        assert "Parrinello-Rahman" in text
        assert "V-rescale" not in text and "C-rescale" not in text
