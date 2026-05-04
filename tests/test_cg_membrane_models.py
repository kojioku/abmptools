# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.models dataclasses."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from abmptools.cg.membrane.models import (
    EquilibrationCGProtocol,
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
    PullingCGProtocol,
    UmbrellaCGProtocol,
)


# ---------------------------------------------------------------------------
# LipidMix
# ---------------------------------------------------------------------------

class TestLipidMix:
    def test_defaults(self):
        lp = LipidMix()
        assert lp.resname == "POPC"
        assert lp.n_per_leaflet == 64

    def test_custom(self):
        lp = LipidMix(resname="DOPC", n_per_leaflet=32)
        assert lp.resname == "DOPC"
        assert lp.n_per_leaflet == 32

    def test_empty_resname_raises(self):
        with pytest.raises(ValueError, match="resname"):
            LipidMix(resname="")

    def test_zero_n_raises(self):
        with pytest.raises(ValueError, match="n_per_leaflet"):
            LipidMix(n_per_leaflet=0)

    def test_negative_n_raises(self):
        with pytest.raises(ValueError, match="n_per_leaflet"):
            LipidMix(n_per_leaflet=-1)


# ---------------------------------------------------------------------------
# PeptideMembraneSpec
# ---------------------------------------------------------------------------

class TestPeptideMembraneSpec:
    def test_default_kgg(self):
        spec = PeptideMembraneSpec()
        assert spec.name == "kgg"
        assert spec.sequence == "KGG"
        assert spec.initial_z_offset_nm == 3.0

    def test_lowercase_sequence_normalized(self):
        spec = PeptideMembraneSpec(sequence="kgg")
        assert spec.sequence == "KGG"

    def test_invalid_aa_raises(self):
        with pytest.raises(ValueError, match="Invalid amino acid"):
            PeptideMembraneSpec(sequence="KGZ")

    def test_empty_name_raises(self):
        with pytest.raises(ValueError, match="name"):
            PeptideMembraneSpec(name="", sequence="KGG")

    def test_no_seq_no_pdb_raises(self):
        with pytest.raises(ValueError, match="sequence.*cg_pdb_path"):
            PeptideMembraneSpec(sequence="")

    def test_seq_and_pdb_raises(self):
        with pytest.raises(ValueError, match="exactly one"):
            PeptideMembraneSpec(
                sequence="KGG",
                cg_pdb_path="/tmp/a.pdb",
                cg_itp_path="/tmp/a.itp",
            )

    def test_pdb_without_itp_raises(self):
        with pytest.raises(ValueError, match="both"):
            PeptideMembraneSpec(
                sequence="",
                cg_pdb_path="/tmp/a.pdb",
            )

    def test_prebuilt_cg_files_ok(self):
        spec = PeptideMembraneSpec(
            sequence="",
            cg_pdb_path="/tmp/peptide_cg.pdb",
            cg_itp_path="/tmp/peptide.itp",
        )
        assert spec.cg_pdb_path == "/tmp/peptide_cg.pdb"
        assert spec.cg_itp_path == "/tmp/peptide.itp"


# ---------------------------------------------------------------------------
# Protocols
# ---------------------------------------------------------------------------

class TestProtocols:
    def test_equil_defaults(self):
        eq = EquilibrationCGProtocol()
        assert eq.dt_ps == 0.020
        assert eq.temperature_K == 310.0
        assert eq.tau_p_ps == 12.0  # Martini standard

    def test_equil_negative_dt_raises(self):
        with pytest.raises(ValueError, match="dt_ps"):
            EquilibrationCGProtocol(dt_ps=0.0)

    def test_equil_negative_temperature_raises(self):
        with pytest.raises(ValueError, match="temperature"):
            EquilibrationCGProtocol(temperature_K=-1.0)

    def test_equil_negative_steps_raises(self):
        with pytest.raises(ValueError, match="em_steps"):
            EquilibrationCGProtocol(em_steps=-1)

    def test_pulling_defaults(self):
        pl = PullingCGProtocol()
        assert pl.pull_force_constant == 1000.0
        assert pl.pull_rate_nm_per_ps == 0.001
        assert pl.nsteps == 250_000  # 5 ns at dt=20 fs

    def test_pulling_negative_k_raises(self):
        with pytest.raises(ValueError, match="pull_force_constant"):
            PullingCGProtocol(pull_force_constant=-1.0)

    def test_umbrella_defaults_yield_13_windows(self):
        umb = UmbrellaCGProtocol()
        assert umb.n_windows == 13   # -1.5 .. +1.5 step 0.25 inclusive
        assert umb.force_constant_kj_mol_nm2 == 1000.0
        assert umb.window_nsteps == 50_000  # 1 ns at dt=20 fs

    def test_umbrella_zmin_geq_zmax_raises(self):
        with pytest.raises(ValueError, match="z_min_nm"):
            UmbrellaCGProtocol(z_min_nm=2.0, z_max_nm=1.0)

    def test_umbrella_zero_spacing_raises(self):
        with pytest.raises(ValueError, match="window_spacing_nm"):
            UmbrellaCGProtocol(window_spacing_nm=0.0)

    def test_umbrella_negative_k_raises(self):
        with pytest.raises(ValueError, match="force_constant"):
            UmbrellaCGProtocol(force_constant_kj_mol_nm2=-1.0)

    def test_umbrella_pull_defaults(self):
        umb = UmbrellaCGProtocol()
        assert umb.pull_geometry == "direction"
        assert umb.pull_group1 == "Bilayer"
        assert umb.pull_group2 == "Peptide"
        assert umb.pull_dim == "N N Y"


# ---------------------------------------------------------------------------
# MembraneCGBuildConfig
# ---------------------------------------------------------------------------

class TestMembraneCGBuildConfig:
    def _minimal_cfg(self):
        return MembraneCGBuildConfig(
            lipids=[LipidMix()],
            peptide=PeptideMembraneSpec(),
            martini_itp_dir="./ff",
        )

    def test_defaults(self):
        cfg = self._minimal_cfg()
        assert cfg.insane_d_nm == 8.0
        assert cfg.box_z_nm == 14.0
        assert cfg.solvent_type == "W"
        assert cfg.nacl_molar == 0.15
        assert cfg.insane_path == "insane"
        assert cfg.grompp_maxwarn == 5

    def test_empty_lipids_raises(self):
        with pytest.raises(ValueError, match="lipids"):
            MembraneCGBuildConfig(
                lipids=[],
                peptide=PeptideMembraneSpec(),
            )

    def test_no_peptide_raises(self):
        with pytest.raises(ValueError, match="peptide"):
            MembraneCGBuildConfig(
                lipids=[LipidMix()],
                peptide=None,
            )

    def test_multiple_lipids_v1_raises(self):
        with pytest.raises(ValueError, match="single lipid"):
            MembraneCGBuildConfig(
                lipids=[LipidMix(resname="POPC"), LipidMix(resname="POPE")],
                peptide=PeptideMembraneSpec(),
            )

    def test_negative_d_raises(self):
        with pytest.raises(ValueError, match="insane_d_nm"):
            MembraneCGBuildConfig(
                lipids=[LipidMix()],
                peptide=PeptideMembraneSpec(),
                insane_d_nm=-1.0,
            )

    def test_negative_box_z_raises(self):
        with pytest.raises(ValueError, match="box_z_nm"):
            MembraneCGBuildConfig(
                lipids=[LipidMix()],
                peptide=PeptideMembraneSpec(),
                box_z_nm=0.0,
            )

    def test_negative_nacl_raises(self):
        with pytest.raises(ValueError, match="nacl_molar"):
            MembraneCGBuildConfig(
                lipids=[LipidMix()],
                peptide=PeptideMembraneSpec(),
                nacl_molar=-0.1,
            )

    def test_negative_maxwarn_raises(self):
        with pytest.raises(ValueError, match="grompp_maxwarn"):
            MembraneCGBuildConfig(
                lipids=[LipidMix()],
                peptide=PeptideMembraneSpec(),
                grompp_maxwarn=-1,
            )

    def test_to_json_from_json_roundtrip(self, tmp_path):
        cfg = self._minimal_cfg()
        cfg.umbrella.window_nsteps = 100  # tweak to verify round-trip
        path = tmp_path / "cfg.json"
        cfg.to_json(str(path))
        loaded = MembraneCGBuildConfig.from_json(str(path))
        assert loaded.peptide.sequence == "KGG"
        assert loaded.peptide.name == "kgg"
        assert loaded.lipids[0].resname == "POPC"
        assert loaded.lipids[0].n_per_leaflet == 64
        assert loaded.umbrella.window_nsteps == 100
        assert loaded.equilibration.temperature_K == 310.0
        assert loaded.pulling.pull_force_constant == 1000.0

    def test_json_file_is_valid_json(self, tmp_path):
        cfg = self._minimal_cfg()
        path = tmp_path / "cfg.json"
        cfg.to_json(str(path))
        # parse manually to confirm valid JSON
        data = json.loads(path.read_text())
        assert "lipids" in data
        assert "peptide" in data
        assert "umbrella" in data

    def test_n_windows_propagates(self):
        cfg = self._minimal_cfg()
        cfg.umbrella.z_min_nm = -2.0
        cfg.umbrella.z_max_nm = +2.0
        cfg.umbrella.window_spacing_nm = 0.5
        # 4.0 / 0.5 = 8 + 1 = 9
        assert cfg.umbrella.n_windows == 9
