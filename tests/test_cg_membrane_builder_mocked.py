# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.builder.MembraneCGBuilder (mocked stages).

All external subprocess calls (insane, gmx via cg.peptide sub-call) are
patched with file-creating fakes. Verifies pipeline orchestration,
directory layout, and result dict.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from abmptools.cg.membrane.builder import MembraneCGBuilder
from abmptools.cg.membrane.models import (
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
)


def _build_config(tmp_path):
    return MembraneCGBuildConfig(
        lipids=[LipidMix(resname="POPC", n_per_leaflet=32)],
        peptide=PeptideMembraneSpec(name="kgg", sequence="KGG"),
        martini_itp_dir=str(tmp_path / "ff"),
        output_dir=str(tmp_path),
    )


def _stage_mocks(monkeypatch, tmp_path):
    """Patch insane / cg.peptide / gmx with file-creating fakes."""

    # Materialize the 4 ITPs that _copy_ff_files expects.
    ff_src = tmp_path / "ff_src"
    ff_src.mkdir(parents=True, exist_ok=True)
    for fname in [
        "martini_v3.0.0.itp",
        "martini_v3.0.0_solvents_v1.itp",
        "martini_v3.0.0_ions_v1.itp",
        "martini_v3.0.0_phospholipids_v1.itp",
    ]:
        (ff_src / fname).write_text("; stub itp\n")

    # 1. Mock cg.peptide PeptideCGBuilder so stage1 produces stub CG files.
    def fake_peptide_build(self):
        # Mirror real cg.peptide output layout (sub-call writes into a
        # dedicated scratch dir under the membrane builder's output_dir).
        out = Path(self.config.output_dir).resolve()
        spec = self.config.peptides[0]
        cg_dir = out / "molecules" / spec.name
        cg_dir.mkdir(parents=True, exist_ok=True)
        cg_pdb = cg_dir / f"{spec.name}_cg.pdb"
        cg_pdb.write_text("REMARK stub CG PDB\n")
        itp = cg_dir / f"{spec.name}.itp"
        itp.write_text(
            "[ moleculetype ]\nmolecule_0 1\n[ atoms ]\n"
        )
        atomistic = cg_dir / f"{spec.name}_atomistic.pdb"
        atomistic.write_text("REMARK stub atomistic\n")
        return {"output_dir": out, "gro": cg_pdb, "top": cg_pdb,
                "ndx": cg_pdb}

    monkeypatch.setattr(
        "abmptools.cg.membrane.builder.PeptideCGBuilder.build",
        fake_peptide_build,
    )

    # 2. Mock insane to materialize bilayer.gro + insane_topol.top.
    def fake_run_insane(output_dir, **kwargs):
        out_dir = Path(output_dir).resolve()
        gro = out_dir / "bilayer.gro"
        gro.write_text(
            "stub insane gro\n"
            "    3\n"
            f"{1:5d}{'LYS':<5s}{'BB':>5s}{1:5d}"
            f"{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}\n"
            f"{2:5d}{'POPC':<5s}{'NC3':>5s}{2:5d}"
            f"{4.0:8.3f}{5.0:8.3f}{6.0:8.3f}\n"
            f"{3:5d}{'NA+':<5s}{'NA+':>5s}{3:5d}"
            f"{7.0:8.3f}{8.0:8.3f}{9.0:8.3f}\n"
            "  10.000  10.000  10.000\n"
        )
        top = out_dir / "insane_topol.top"
        top.write_text(
            '#include "martini.itp"\n'
            '\n'
            '[ system ]\n'
            'INSANE!\n'
            '\n'
            '[ molecules ]\n'
            '; name  number\n'
            'Protein 1\n'
            'POPC    32\n'
            'POPC    32\n'
            'W       100\n'
            'NA+     2\n'
            'CL-     3\n'
        )
        return gro, top

    monkeypatch.setattr(
        "abmptools.cg.membrane.builder.insane_runner.run_insane",
        fake_run_insane,
    )

    return ff_src


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestMembraneCGBuilder:

    def test_build_returns_expected_keys(self, tmp_path, monkeypatch):
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)

        result = MembraneCGBuilder(cfg).build()

        for key in (
            "output_dir", "gro", "top", "ndx", "equil_mdps", "pull_mdp",
            "window_mdps", "run_script", "config_json", "n_windows",
            "lipid_resname", "n_lipids_total",
        ):
            assert key in result, f"missing key: {key}"
        assert result["n_windows"] == 13
        assert result["lipid_resname"] == "POPC"
        assert result["n_lipids_total"] == 64  # 2 * 32

    def test_build_creates_directory_structure(
        self, tmp_path, monkeypatch,
    ):
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        MembraneCGBuilder(cfg).build()

        assert (tmp_path / "input" / "config.json").exists()
        assert (tmp_path / "molecules" / "kgg" / "kgg.itp").exists()
        assert (tmp_path / "molecules" / "kgg" / "kgg_cg.pdb").exists()
        assert (tmp_path / "topol.top").exists()
        assert (tmp_path / "system_ions.gro").exists()
        assert (tmp_path / "index.ndx").exists()
        assert (tmp_path / "mdp" / "em.mdp").exists()
        assert (tmp_path / "mdp" / "nvt.mdp").exists()
        assert (tmp_path / "mdp" / "npt.mdp").exists()
        assert (tmp_path / "pull" / "pull.mdp").exists()
        assert (tmp_path / "windows" / "win_000" / "window.mdp").exists()
        assert (tmp_path / "windows" / "win_012" / "window.mdp").exists()
        assert (tmp_path / "run.sh").exists()

    def test_run_script_executable(self, tmp_path, monkeypatch):
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        result = MembraneCGBuilder(cfg).build()
        run_sh = result["run_script"]
        assert run_sh.stat().st_mode & 0o100

    def test_topology_includes_4_itps_and_peptide(
        self, tmp_path, monkeypatch,
    ):
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        MembraneCGBuilder(cfg).build()
        topol = (tmp_path / "topol.top").read_text()
        assert '#include "martini_v3.0.0.itp"' in topol
        assert '#include "martini_v3.0.0_phospholipids_v1.itp"' in topol
        assert '#include "molecules/kgg/kgg.itp"' in topol
        assert "molecule_0" in topol
        # Ion names normalized in molecules block
        mol_block = topol.split("[ molecules ]")[1]
        assert "NA+" not in mol_block
        assert "CL-" not in mol_block

    def test_ions_normalized_in_gro(self, tmp_path, monkeypatch):
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        result = MembraneCGBuilder(cfg).build()
        gro_text = result["gro"].read_text()
        assert "NA+" not in gro_text
        # NA still present (normalized), POPC + LYS preserved
        assert "POPC" in gro_text

    def test_index_has_named_groups(self, tmp_path, monkeypatch):
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        result = MembraneCGBuilder(cfg).build()
        ndx_text = result["ndx"].read_text()
        assert "[ Bilayer ]" in ndx_text
        assert "[ Peptide ]" in ndx_text
        assert "[ NA ]" in ndx_text
        assert "[ Non_Bilayer ]" in ndx_text

    def test_no_peptide_raises(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch, tmp_path)
        # Construct config bypassing dataclass validation by toggling after
        # construction is impossible (peptide=None raises). Use a custom
        # builder that nullifies peptide post-construction.
        cfg = _build_config(tmp_path)
        builder = MembraneCGBuilder(cfg)
        builder.config.peptide = None
        with pytest.raises(ValueError, match="peptide"):
            builder.build()

    def test_config_json_round_trip(self, tmp_path, monkeypatch):
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        MembraneCGBuilder(cfg).build()
        loaded = MembraneCGBuildConfig.from_json(
            str(tmp_path / "input" / "config.json")
        )
        assert loaded.peptide.sequence == "KGG"
        assert loaded.peptide.name == "kgg"
        assert loaded.lipids[0].resname == "POPC"
        assert loaded.lipids[0].n_per_leaflet == 32

    def test_prebuilt_cg_files_skip_subcall(self, tmp_path, monkeypatch):
        """When user supplies cg_pdb_path + cg_itp_path, peptide subcall is skipped."""
        ff_src = _stage_mocks(monkeypatch, tmp_path)

        # Counter to detect cg.peptide.PeptideCGBuilder.build invocation.
        counter = {"n_calls": 0}
        original = monkeypatch
        # Re-mock to count calls
        from abmptools.cg.membrane.builder import PeptideCGBuilder

        def counting_build(self):
            counter["n_calls"] += 1
            return {}

        monkeypatch.setattr(PeptideCGBuilder, "build", counting_build)

        # Pre-build CG files.
        prebuilt = tmp_path / "prebuilt"
        prebuilt.mkdir()
        cg_pdb = prebuilt / "kgg_cg.pdb"
        cg_pdb.write_text("REMARK stub\n")
        cg_itp = prebuilt / "kgg.itp"
        cg_itp.write_text(
            "[ moleculetype ]\nmolecule_0 1\n[ atoms ]\n"
        )

        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        cfg.peptide = PeptideMembraneSpec(
            name="kgg",
            sequence="",
            cg_pdb_path=str(cg_pdb),
            cg_itp_path=str(cg_itp),
        )
        MembraneCGBuilder(cfg).build()

        # cg.peptide sub-call should NOT have been invoked
        assert counter["n_calls"] == 0

    def test_input_dir_isolated_from_molecules(self, tmp_path, monkeypatch):
        """input/ and molecules/ are separate directories."""
        ff_src = _stage_mocks(monkeypatch, tmp_path)
        cfg = _build_config(tmp_path)
        cfg.martini_itp_dir = str(ff_src)
        MembraneCGBuilder(cfg).build()
        assert (tmp_path / "input").is_dir()
        assert (tmp_path / "molecules").is_dir()
        # config.json lives in input/, peptide ITP in molecules/
        assert (tmp_path / "input" / "config.json").exists()
        assert (tmp_path / "molecules" / "kgg" / "kgg.itp").exists()
