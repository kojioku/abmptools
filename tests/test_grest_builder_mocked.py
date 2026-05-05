# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.builder (mocked subprocess)."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.genesis.grest.builder import GrestBuilder
from abmptools.genesis.grest.models import (
    GrestBuildConfig,
    ReplicaTemperatureSpec,
    RESTSelectionSpec,
)


@pytest.fixture
def fake_prmtop_text() -> str:
    """Synthetic prmtop with BOX_DIMENSIONS + RESIDUE_LABEL."""
    return (
        "%VERSION  VERSION_STAMP = V0001.000\n"
        "%FLAG BOX_DIMENSIONS\n"
        "%FORMAT(5E16.8)\n"
        "  9.0E+01  8.10568E+01  8.35179E+01  9.29056E+01\n"
        "%FLAG RESIDUE_LABEL\n"
        "%FORMAT(20a4)\n"
        "LYS GLY GLY WAT WAT WAT WAT \n"
        "%FLAG NEXT\n"
    )


@pytest.fixture
def explicit_cfg(tmp_path) -> GrestBuildConfig:
    """Explicit-mode config that needs only tleap (no cpptraj)."""
    return GrestBuildConfig(
        input_pdb=str(tmp_path / "protein.pdb"),
        rest_selection=RESTSelectionSpec(
            mode="explicit", residues=["1-3"]
        ),
        replica_temperatures=ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.0, 318.11, 337.11, 357.10],
        ),
        output_dir=str(tmp_path / "out"),
    )


def _make_fake_run(workdir: Path, prmtop_text: str):
    """Return a side_effect callable mocking tleap + cpptraj."""
    def fake_run(cmd, cwd=None, capture=True, **kwargs):
        cwd = Path(cwd) if cwd else workdir
        # tleap path: write prmtop + coor + ref pdb.
        if cmd[0].endswith("tleap") or "tleap" in cmd[0]:
            (cwd / "system.prmtop").write_text(prmtop_text)
            (cwd / "system.coor").write_text("# fake coor\n")
            (cwd / "system_ref.pdb").write_text("# fake ref pdb\n")
        elif cmd[0].endswith("cpptraj") or "cpptraj" in cmd[0]:
            # cpptraj path: write rest_residues.txt.
            (cwd / "rest_residues.txt").write_text(
                "  21 ARG 100 123 23 21 1\n"
                "  96 LYS 1418 1439 22 96 1\n"
            )
        from subprocess import CompletedProcess
        return CompletedProcess(cmd, 0, "ok", "")
    return fake_run


# ---------------------------------------------------------------------------
# Full-pipeline test (explicit REST mode -- only tleap is shelled out)
# ---------------------------------------------------------------------------

class TestGrestBuilderExplicitMode:
    def test_full_pipeline_runs(
        self, explicit_cfg, fake_prmtop_text, tmp_path,
    ):
        # Build dir is auto-created under output_dir; the fake_run
        # writes into whatever cwd it's invoked with.
        build_dir = Path(explicit_cfg.output_dir) / "build"
        # Ensure the dir exists for the fake side_effect to write into.
        build_dir.mkdir(parents=True, exist_ok=True)
        fake_run = _make_fake_run(build_dir, fake_prmtop_text)

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run,
        ):
            builder = GrestBuilder(explicit_cfg)
            result = builder.build()

        # Result dict has the expected keys.
        for key in [
            "output_dir",
            "build_dir",
            "inp_dir",
            "config_json",
            "prmtop",
            "coor",
            "ref_pdb",
            "leap_log",
            "box_size_A",
            "n_protein_residues",
            "rest_residues",
            "rest_selection_string",
            "n_rest_residues",
            "temperature_ladder",
            "n_replicas",
            "inp_files",
            "run_script",
            "run_pjsub",
            "run_sbatch",
        ]:
            assert key in result, f"missing key: {key}"

    def test_output_layout(self, explicit_cfg, fake_prmtop_text, tmp_path):
        build_dir = Path(explicit_cfg.output_dir) / "build"
        build_dir.mkdir(parents=True, exist_ok=True)
        fake_run = _make_fake_run(build_dir, fake_prmtop_text)

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run,
        ):
            builder = GrestBuilder(explicit_cfg)
            result = builder.build()

        out = Path(result["output_dir"])
        # Expected files exist.
        assert (out / "input" / "config.json").exists()
        assert (out / "build" / "system.prmtop").exists()
        assert (out / "build" / "system.coor").exists()
        assert (out / "build" / "system_ref.pdb").exists()
        assert (out / "build" / "system_tleap.log").exists()
        assert (out / "build" / "system.tleap").exists()
        assert (out / "inp" / "step1_minimize.inp").exists()
        assert (out / "inp" / "step2_equilibrate.inp").exists()
        assert (out / "inp" / "step3_grest.inp").exists()
        assert (out / "inp" / "step5_remd_convert.inp").exists()
        assert (out / "run.sh").exists()
        assert (out / "run_pjsub.sh").exists()
        assert (out / "run_sbatch.sh").exists()

    def test_inp_files_have_correct_box_size(
        self, explicit_cfg, fake_prmtop_text, tmp_path,
    ):
        build_dir = Path(explicit_cfg.output_dir) / "build"
        build_dir.mkdir(parents=True, exist_ok=True)
        fake_run = _make_fake_run(build_dir, fake_prmtop_text)

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run,
        ):
            builder = GrestBuilder(explicit_cfg)
            result = builder.build()

        # Box sizes from prmtop are baked into the minimize.inp.
        text = result["inp_files"]["minimize"].read_text()
        assert "box_size_x = 81.056800" in text or "box_size_x = 81.056806" in text or "box_size_x = 81.05680" in text
        # Replica count flows through.
        grest_text = result["inp_files"]["grest"].read_text()
        assert "nreplica1       = 4" in grest_text

    def test_rest_selection_string_formatted(
        self, explicit_cfg, fake_prmtop_text,
    ):
        build_dir = Path(explicit_cfg.output_dir) / "build"
        build_dir.mkdir(parents=True, exist_ok=True)
        fake_run = _make_fake_run(build_dir, fake_prmtop_text)

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run,
        ):
            builder = GrestBuilder(explicit_cfg)
            result = builder.build()

        # explicit mode "1-3" -> "rno:1-3"
        assert result["rest_selection_string"] == "rno:1-3"
        assert result["rest_residues"] == [1, 2, 3]
        assert result["n_rest_residues"] == 3

    def test_run_sh_executable(self, explicit_cfg, fake_prmtop_text):
        build_dir = Path(explicit_cfg.output_dir) / "build"
        build_dir.mkdir(parents=True, exist_ok=True)
        fake_run = _make_fake_run(build_dir, fake_prmtop_text)

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run,
        ):
            builder = GrestBuilder(explicit_cfg)
            result = builder.build()

        run_sh = Path(result["run_script"])
        assert run_sh.exists()
        # 0o755 (rwxr-xr-x) was set by write_run_script.
        assert run_sh.stat().st_mode & 0o111  # at least one exec bit set

    def test_temperature_ladder_passthrough(
        self, explicit_cfg, fake_prmtop_text,
    ):
        build_dir = Path(explicit_cfg.output_dir) / "build"
        build_dir.mkdir(parents=True, exist_ok=True)
        fake_run = _make_fake_run(build_dir, fake_prmtop_text)

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run,
        ):
            builder = GrestBuilder(explicit_cfg)
            result = builder.build()

        assert result["temperature_ladder"] == [300.0, 318.11, 337.11, 357.10]
        assert result["n_replicas"] == 4


# ---------------------------------------------------------------------------
# around mode (mocks both tleap + cpptraj)
# ---------------------------------------------------------------------------

class TestGrestBuilderAroundMode:
    def test_around_mode_runs_cpptraj(self, tmp_path, fake_prmtop_text):
        cfg = GrestBuildConfig(
            input_pdb=str(tmp_path / "protein.pdb"),
            rest_selection=RESTSelectionSpec(
                mode="around", center="rno:96", radius_A=5.0
            ),
            replica_temperatures=ReplicaTemperatureSpec(
                mode="manual",
                temperatures=[300.0, 318.11, 337.11, 357.10],
            ),
            output_dir=str(tmp_path / "out"),
        )
        build_dir = Path(cfg.output_dir) / "build"
        build_dir.mkdir(parents=True, exist_ok=True)

        # Both tleap and cpptraj invoked here -- we patch run_command in
        # both modules.
        fake_run_sb = _make_fake_run(build_dir, fake_prmtop_text)
        fake_run_rs = _make_fake_run(build_dir, fake_prmtop_text)

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run_sb,
        ), patch(
            "abmptools.genesis.grest.rest_selection.run_command",
            side_effect=fake_run_rs,
        ):
            builder = GrestBuilder(cfg)
            result = builder.build()

        # cpptraj produced rest_residues.txt with [21, 96] -> compacted.
        assert result["rest_residues"] == [21, 96]
        assert result["rest_selection_string"] == "rno:21,96"
