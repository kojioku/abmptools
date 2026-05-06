# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.mmgbsa.{inp_writer,gbsa_runner}."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.genesis.mmgbsa.gbsa_runner import (
    GBSARunResult,
    run_gbsa_one_target,
)
from abmptools.genesis.mmgbsa.inp_writer import (
    SYSTEM_NAMES,
    render_gbsa_inp,
    render_three_inps,
)
from abmptools.genesis.mmgbsa.models import (
    EnergyProtocol,
    MinimizationProtocol,
)


# ---------------------------------------------------------------------------
# render_gbsa_inp
# ---------------------------------------------------------------------------

class TestRenderGbsaInp:
    def test_complex_default_blocks(self):
        text = render_gbsa_inp(
            "complex",
            EnergyProtocol(),
            MinimizationProtocol(),
        )
        for section in [
            "[INPUT]", "[OUTPUT]", "[ENERGY]", "[MINIMIZE]", "[BOUNDARY]",
        ]:
            assert section in text

    def test_amber_force_field_and_gbsa(self):
        text = render_gbsa_inp("ligand", EnergyProtocol(), MinimizationProtocol())
        assert "forcefield       = AMBER" in text
        assert "implicit_solvent = GBSA" in text
        # POC values.
        assert "gbsa_salt_cons   = 0.1500" in text
        assert "gbsa_surf_tens   = 0.007200" in text
        assert "gbsa_vdw_offset  = 0.0900" in text

    def test_nobc_boundary(self):
        text = render_gbsa_inp("receptor", EnergyProtocol(), MinimizationProtocol())
        assert "type             = NOBC" in text

    def test_cutoff_99_9_default(self):
        text = render_gbsa_inp("complex", EnergyProtocol(), MinimizationProtocol())
        assert "cutoffdist       = 99.9000" in text
        assert "switchdist       = 99.9000" in text
        assert "pairlistdist     = 100.0000" in text

    def test_minimize_singlepoint(self):
        text = render_gbsa_inp("complex", EnergyProtocol(), MinimizationProtocol())
        assert "method           = SD" in text
        assert "nsteps           = 1" in text  # POC: single-point

    def test_implicit_none_skips_gbsa_keys(self):
        text = render_gbsa_inp(
            "complex",
            EnergyProtocol(implicit_solvent="NONE"),
            MinimizationProtocol(),
        )
        assert "implicit_solvent = NONE" in text
        # GBSA-specific keys should NOT appear.
        assert "gbsa_salt_cons" not in text
        assert "gbsa_surf_tens" not in text

    def test_pme_alternate(self):
        text = render_gbsa_inp(
            "complex",
            EnergyProtocol(electrostatic="PME"),
            MinimizationProtocol(),
        )
        assert "electrostatic    = PME" in text

    def test_filenames_match_target_name(self):
        text = render_gbsa_inp("ligand", EnergyProtocol(), MinimizationProtocol())
        assert "prmtopfile = ligand.prmtop" in text
        assert "ambcrdfile = ligand.inpcrd" in text
        assert "dcdfile = ligand.dcd" in text
        assert "rstfile = ligand.rst" in text

    def test_invalid_name_raises(self):
        with pytest.raises(ValueError, match="name must be"):
            render_gbsa_inp("solvent", EnergyProtocol(), MinimizationProtocol())


# ---------------------------------------------------------------------------
# render_three_inps
# ---------------------------------------------------------------------------

class TestRenderThreeInps:
    def test_keys(self):
        d = render_three_inps(EnergyProtocol(), MinimizationProtocol())
        assert set(d.keys()) == set(SYSTEM_NAMES)

    def test_distinct_filenames(self):
        d = render_three_inps(EnergyProtocol(), MinimizationProtocol())
        assert "complex.prmtop" in d["complex"]
        assert "ligand.prmtop" in d["ligand"]
        assert "receptor.prmtop" in d["receptor"]


# ---------------------------------------------------------------------------
# run_gbsa_one_target (mocked atdyn)
# ---------------------------------------------------------------------------

class TestRunGbsaOneTarget:
    def test_full_flow_mocked(self, tmp_path):
        target_dir = tmp_path / "target"
        target_dir.mkdir()
        # Pretend prmtop+inpcrd already in place from Stage 2.
        for name in SYSTEM_NAMES:
            (target_dir / f"{name}.prmtop").write_text("# fake\n")
            (target_dir / f"{name}.inpcrd").write_text("# fake\n")

        inps = render_three_inps(EnergyProtocol(), MinimizationProtocol())

        captured: list = []

        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            captured.append(list(cmd))
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "[STEP4] fake output\n", "")

        with patch(
            "abmptools.genesis.mmgbsa.gbsa_runner.run_command",
            side_effect=fake_run,
        ):
            result = run_gbsa_one_target(
                target_dir=target_dir,
                inps=inps,
                atdyn_path="atdyn",
                mpirun_path="mpirun",
                mpi_processes=1,
            )

        # 3 inp files + 3 log files materialised.
        assert isinstance(result, GBSARunResult)
        for name in SYSTEM_NAMES:
            assert result.inp_files[name].exists()
            assert result.log_files[name].exists()
            assert "[STEP4]" in result.log_files[name].read_text()

        # 3 mpirun calls, each with -np 1 + atdyn + the right inp.
        assert len(captured) == 3
        for cmd in captured:
            assert cmd[0] == "mpirun"
            assert cmd[1:3] == ["-np", "1"]
            assert cmd[3] == "atdyn"

    def test_missing_key_raises(self, tmp_path):
        target_dir = tmp_path / "target"
        target_dir.mkdir()
        with pytest.raises(ValueError, match="missing keys"):
            run_gbsa_one_target(
                target_dir=target_dir,
                inps={"complex": "x", "ligand": "y"},  # receptor missing
            )

    def test_empty_log_raises(self, tmp_path):
        target_dir = tmp_path / "target"
        target_dir.mkdir()
        inps = render_three_inps(EnergyProtocol(), MinimizationProtocol())

        def fake_run_empty(cmd, cwd=None, capture=True, **kwargs):
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "", "")  # empty stdout

        from abmptools.genesis.mmgbsa._subprocess import CommandError
        with patch(
            "abmptools.genesis.mmgbsa.gbsa_runner.run_command",
            side_effect=fake_run_empty,
        ):
            with pytest.raises(CommandError, match="empty"):
                run_gbsa_one_target(target_dir=target_dir, inps=inps)

    def test_higher_mpi_processes_passed(self, tmp_path):
        target_dir = tmp_path / "target"
        target_dir.mkdir()
        inps = render_three_inps(EnergyProtocol(), MinimizationProtocol())
        captured: list = []

        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            captured.append(list(cmd))
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "[STEP4] ok\n", "")

        with patch(
            "abmptools.genesis.mmgbsa.gbsa_runner.run_command",
            side_effect=fake_run,
        ):
            run_gbsa_one_target(
                target_dir=target_dir,
                inps=inps,
                mpi_processes=4,
            )
        assert all(cmd[1:3] == ["-np", "4"] for cmd in captured)
