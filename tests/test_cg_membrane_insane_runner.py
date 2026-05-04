# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.insane_runner.run_insane."""
from __future__ import annotations

from pathlib import Path
from typing import List
from unittest.mock import MagicMock

import pytest

from abmptools.cg.membrane import insane_runner


@pytest.fixture
def captured_cmd(monkeypatch, tmp_path):
    """Patch run_command and capture argv. Materialize output files as stubs.

    Returns a list that the test can read after calling ``run_insane``.
    """
    captured: List[List[str]] = []

    def fake_run_command(cmd, cwd=None, **kw):
        captured.append(list(cmd))
        # Materialize -o gro and -p top so post-checks pass
        for i, tok in enumerate(cmd):
            if tok == "-o" and i + 1 < len(cmd):
                Path(cmd[i + 1]).parent.mkdir(parents=True, exist_ok=True)
                Path(cmd[i + 1]).write_text("; stub gro\n")
            elif tok == "-p" and i + 1 < len(cmd):
                Path(cmd[i + 1]).parent.mkdir(parents=True, exist_ok=True)
                Path(cmd[i + 1]).write_text("; stub top\n")
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(insane_runner, "run_command", fake_run_command)
    return captured


# ---------------------------------------------------------------------------
# Bilayer-only invocation (no peptide)
# ---------------------------------------------------------------------------

class TestBilayerOnly:
    def test_basic_argv(self, captured_cmd, tmp_path):
        gro, top = insane_runner.run_insane(
            output_dir=tmp_path,
            lipid_resname="POPC",
            n_per_leaflet=64,
            box_d_nm=8.0,
            box_z_nm=14.0,
        )
        assert gro.exists()
        assert top.exists()
        argv = captured_cmd[0]
        assert "insane" in argv[0]
        assert "-l" in argv and "POPC" in argv
        assert "-d" in argv
        assert "-dz" in argv
        assert "-sol" in argv and "W" in argv
        assert "-salt" in argv
        assert "-pbc" in argv and "hexagonal" in argv

    def test_no_peptide_args(self, captured_cmd, tmp_path):
        insane_runner.run_insane(
            output_dir=tmp_path,
            lipid_resname="POPC",
            n_per_leaflet=64,
            box_d_nm=8.0, box_z_nm=14.0,
        )
        argv = captured_cmd[0]
        # bilayer-only: no -f, no -dm, no -charge
        assert "-f" not in argv
        assert "-dm" not in argv
        assert "-charge" not in argv

    def test_custom_lipid_resname(self, captured_cmd, tmp_path):
        insane_runner.run_insane(
            output_dir=tmp_path,
            lipid_resname="DOPC",
            n_per_leaflet=32,
            box_d_nm=6.0, box_z_nm=10.0,
        )
        argv = captured_cmd[0]
        l_idx = argv.index("-l")
        assert argv[l_idx + 1] == "DOPC"

    def test_custom_solvent(self, captured_cmd, tmp_path):
        insane_runner.run_insane(
            output_dir=tmp_path,
            lipid_resname="POPC",
            n_per_leaflet=64,
            box_d_nm=8.0, box_z_nm=14.0,
            solvent="WN",       # negative-charge water bead variant
        )
        argv = captured_cmd[0]
        s_idx = argv.index("-sol")
        assert argv[s_idx + 1] == "WN"


# ---------------------------------------------------------------------------
# Peptide-in-bilayer invocation
# ---------------------------------------------------------------------------

class TestPeptideInBilayer:
    def test_peptide_argv_includes_f_dm_charge(
        self, captured_cmd, tmp_path,
    ):
        peptide = tmp_path / "peptide_cg.pdb"
        peptide.write_text("REMARK stub\n")
        insane_runner.run_insane(
            output_dir=tmp_path,
            lipid_resname="POPC",
            n_per_leaflet=64,
            box_d_nm=8.0, box_z_nm=14.0,
            peptide_pdb=peptide,
            peptide_z_offset_nm=4.0,
        )
        argv = captured_cmd[0]
        assert "-f" in argv
        assert str(peptide.resolve()) in argv
        assert "-dm" in argv
        dm_idx = argv.index("-dm")
        assert float(argv[dm_idx + 1]) == pytest.approx(4.0)
        assert "-charge" in argv
        ch_idx = argv.index("-charge")
        assert argv[ch_idx + 1] == "auto"

    def test_missing_peptide_pdb_raises(self, captured_cmd, tmp_path):
        with pytest.raises(FileNotFoundError):
            insane_runner.run_insane(
                output_dir=tmp_path,
                lipid_resname="POPC",
                n_per_leaflet=64,
                box_d_nm=8.0, box_z_nm=14.0,
                peptide_pdb=tmp_path / "does_not_exist.pdb",
            )


# ---------------------------------------------------------------------------
# Output file detection / errors
# ---------------------------------------------------------------------------

class TestOutputDetection:
    def test_returns_resolved_paths(self, captured_cmd, tmp_path):
        gro, top = insane_runner.run_insane(
            output_dir=tmp_path,
            lipid_resname="POPC",
            n_per_leaflet=64,
            box_d_nm=8.0, box_z_nm=14.0,
        )
        assert gro.is_absolute()
        assert top.is_absolute()
        assert gro.parent == tmp_path.resolve()

    def test_custom_output_filenames(self, captured_cmd, tmp_path):
        gro, top = insane_runner.run_insane(
            output_dir=tmp_path,
            lipid_resname="POPC",
            n_per_leaflet=64,
            box_d_nm=8.0, box_z_nm=14.0,
            output_gro_name="myinsane.gro",
            output_top_name="myinsane.top",
        )
        assert gro.name == "myinsane.gro"
        assert top.name == "myinsane.top"

    def test_failure_raises_when_output_absent(self, monkeypatch, tmp_path):
        # run_command returns success but doesn't produce files
        def fake(cmd, cwd=None, **kw):
            return MagicMock(returncode=0, stdout="", stderr="")
        monkeypatch.setattr(insane_runner, "run_command", fake)
        with pytest.raises(RuntimeError, match="insane did not produce"):
            insane_runner.run_insane(
                output_dir=tmp_path,
                lipid_resname="POPC",
                n_per_leaflet=64,
                box_d_nm=8.0, box_z_nm=14.0,
            )


# ---------------------------------------------------------------------------
# extra_args
# ---------------------------------------------------------------------------

def test_extra_args_appended_verbatim(captured_cmd, tmp_path):
    insane_runner.run_insane(
        output_dir=tmp_path,
        lipid_resname="POPC",
        n_per_leaflet=64,
        box_d_nm=8.0, box_z_nm=14.0,
        extra_args=["-asym", "0", "-rand", "0.1"],
    )
    argv = captured_cmd[0]
    assert "-asym" in argv
    assert "-rand" in argv
    assert "0.1" in argv


def test_custom_insane_path(captured_cmd, tmp_path):
    insane_runner.run_insane(
        output_dir=tmp_path,
        lipid_resname="POPC",
        n_per_leaflet=64,
        box_d_nm=8.0, box_z_nm=14.0,
        insane_path="/opt/bin/insane",
    )
    argv = captured_cmd[0]
    assert argv[0] == "/opt/bin/insane"
