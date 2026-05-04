# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.builder.PeptideCGBuilder (mocked stages)."""
from __future__ import annotations

from pathlib import Path

import pytest

from abmptools.cg.peptide.builder import PeptideCGBuilder
from abmptools.cg.peptide.models import PeptideBuildConfig, PeptideSpec
from abmptools.cg.peptide.system_packer import PackedSystem


# ---------------------------------------------------------------------------
# Stage mocks
# ---------------------------------------------------------------------------

def _build_config(tmp_path):
    return PeptideBuildConfig(
        peptides=[PeptideSpec(name="kgg", sequence="KGG", count=5)],
        box_size_nm=8.0,
        output_dir=str(tmp_path),
        martini_itp_dir=str(tmp_path / "ff"),
    )


def _stage_mocks(monkeypatch):
    """Patch all 6 stage helpers with file-creating fakes."""

    def fake_atomistic(seq, name, mol_dir, *, tleap_path="tleap"):
        out = mol_dir / f"{name}_atomistic.pdb"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("ATOM 1 N ALA 1\n")
        return out

    def fake_martinize2(input_pdb, output_dir, name,
                        *, martinize2_path="martinize2"):
        itp = output_dir / f"{name}.itp"
        cg = output_dir / f"{name}_cg.pdb"
        itp.write_text("[ moleculetype ]\nProtein_A 1\n")
        cg.write_text("ATOM 1 BB ALA 1\n")
        return itp, cg

    def fake_insert(specs, out_dir, **kw):
        out = out_dir / "packed.gro"
        out.write_text("packed\n")
        return out

    def fake_solvate(packed, topol, out_dir, water, **kw):
        out = out_dir / "system_solv.gro"
        out.write_text("solvated\n")
        return out

    def fake_add_ions(solv, topol, out_dir, **kw):
        ions = out_dir / "system_ions.gro"
        ions.write_text("ions\n")
        ions_mdp = out_dir / "ions.mdp"
        ions_mdp.write_text("integrator = steep\n")
        ions_tpr = out_dir / "ions.tpr"
        ions_tpr.write_text("tpr stub\n")
        return PackedSystem(
            packed_gro=out_dir / "packed.gro",
            solvated_gro=solv,
            ions_gro=ions,
            ions_tpr=ions_tpr,
            ions_mdp=ions_mdp,
        )

    def fake_make_index(gro, out_dir, **kw):
        out = out_dir / "index.ndx"
        out.write_text("[ System ]\n")
        return out

    def fake_water_box(out, **kw):
        out = Path(out)
        out.write_text("water box stub\n")
        return out

    monkeypatch.setattr(
        "abmptools.cg.peptide.builder."
        "peptide_atomistic.build_atomistic_pdb",
        fake_atomistic,
    )
    monkeypatch.setattr(
        "abmptools.cg.peptide.builder.martinize_runner.run_martinize2",
        fake_martinize2,
    )
    monkeypatch.setattr(
        "abmptools.cg.peptide.builder.system_packer.insert_peptides",
        fake_insert,
    )
    monkeypatch.setattr(
        "abmptools.cg.peptide.builder.system_packer.solvate",
        fake_solvate,
    )
    monkeypatch.setattr(
        "abmptools.cg.peptide.builder.system_packer.add_ions",
        fake_add_ions,
    )
    monkeypatch.setattr(
        "abmptools.cg.peptide.builder.system_packer.make_index",
        fake_make_index,
    )
    monkeypatch.setattr(
        "abmptools.cg.peptide.builder.water_box.make_martini_water_box",
        fake_water_box,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestPeptideCGBuilder:

    def test_build_returns_expected_keys(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        result = PeptideCGBuilder(cfg).build()
        for key in (
            "output_dir", "gro", "top", "ndx", "mdp_files",
            "run_script", "config_json", "box_nm", "n_peptides_total",
        ):
            assert key in result, f"missing key: {key}"
        assert result["n_peptides_total"] == 5
        assert result["box_nm"] == (8.0, 8.0, 8.0)

    def test_build_creates_directory_structure(
        self, tmp_path, monkeypatch,
    ):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        PeptideCGBuilder(cfg).build()
        assert (tmp_path / "molecules" / "kgg").exists()
        assert (tmp_path / "molecules" / "kgg" / "kgg_atomistic.pdb").exists()
        assert (tmp_path / "molecules" / "kgg" / "kgg_cg.pdb").exists()
        assert (tmp_path / "molecules" / "kgg" / "kgg.itp").exists()
        assert (tmp_path / "mdp" / "em.mdp").exists()
        assert (tmp_path / "mdp" / "nvt.mdp").exists()
        assert (tmp_path / "mdp" / "npt.mdp").exists()
        assert (tmp_path / "mdp" / "md.mdp").exists()
        assert (tmp_path / "run" / "run.sh").exists()
        assert (tmp_path / "system_ions.gro").exists()
        assert (tmp_path / "topol.top").exists()
        assert (tmp_path / "config.json").exists()

    def test_run_script_is_executable(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        result = PeptideCGBuilder(cfg).build()
        run_sh = result["run_script"]
        assert run_sh.stat().st_mode & 0o100

    def test_run_script_uses_configured_gmx(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        cfg.gmx_path = "/opt/gmx-2024/bin/gmx"
        result = PeptideCGBuilder(cfg).build()
        text = result["run_script"].read_text()
        assert "/opt/gmx-2024/bin/gmx" in text

    def test_solvent_disabled_skips_stages(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        cfg.solvent_enabled = False
        result = PeptideCGBuilder(cfg).build()
        assert result["gro"].name == "packed.gro"

    def test_empty_peptides_raises(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch)
        cfg = PeptideBuildConfig(output_dir=str(tmp_path))
        with pytest.raises(ValueError, match="empty"):
            PeptideCGBuilder(cfg).build()

    def test_config_json_round_trip(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        PeptideCGBuilder(cfg).build()
        loaded = PeptideBuildConfig.from_json(
            str(tmp_path / "config.json")
        )
        assert loaded.peptides[0].sequence == "KGG"
        assert loaded.box_size_nm == 8.0

    def test_topology_includes_m3_files_and_molecule(
        self, tmp_path, monkeypatch,
    ):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        PeptideCGBuilder(cfg).build()
        topol = (tmp_path / "topol.top").read_text()
        assert "martini_v3.0.0.itp" in topol
        assert '#include "molecules/kgg/kgg.itp"' in topol
        # molecule_count
        assert "Protein_A  5" in topol

    def test_box_lengths_takes_precedence(self, tmp_path, monkeypatch):
        _stage_mocks(monkeypatch)
        cfg = _build_config(tmp_path)
        cfg.box_lengths_nm = [10.0, 12.0, 15.0]
        result = PeptideCGBuilder(cfg).build()
        assert result["box_nm"] == (10.0, 12.0, 15.0)

    def test_water_gro_auto_generated_when_missing(
        self, tmp_path, monkeypatch,
    ):
        """ff_dir に water.gro が無いとき auto-gen 経路に入ることを検証。"""
        from abmptools.cg.peptide import water_box

        called = {}
        original_fake = lambda *a, **kw: None  # placeholder

        def tracking_water_box(out, **kw):
            called["out"] = out
            called["gmx_path"] = kw.get("gmx_path", "gmx")
            from pathlib import Path as _P
            _P(out).write_text("auto-gen W box\n")
            return _P(out)

        _stage_mocks(monkeypatch)
        # Override the water_box mock with our tracking version
        monkeypatch.setattr(
            "abmptools.cg.peptide.builder."
            "water_box.make_martini_water_box",
            tracking_water_box,
        )

        cfg = _build_config(tmp_path)
        # ff_dir does not exist, so water.gro is missing -> auto-gen
        PeptideCGBuilder(cfg).build()
        assert "out" in called
        assert called["gmx_path"] == "gmx"

    def test_water_gro_not_auto_generated_when_present(
        self, tmp_path, monkeypatch,
    ):
        """ff_dir に water.gro が既にあれば auto-gen は呼ばれない。"""
        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "martini_v3.0.0_water.gro").write_text("user-provided\n")

        called = {"count": 0}

        def tracking_water_box(out, **kw):
            called["count"] += 1
            from pathlib import Path as _P
            _P(out).write_text("X")
            return _P(out)

        _stage_mocks(monkeypatch)
        monkeypatch.setattr(
            "abmptools.cg.peptide.builder."
            "water_box.make_martini_water_box",
            tracking_water_box,
        )

        cfg = _build_config(tmp_path)
        PeptideCGBuilder(cfg).build()
        assert called["count"] == 0
