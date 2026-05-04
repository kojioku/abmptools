# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.martinize_runner module."""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from abmptools.cg.peptide.martinize_runner import (
    _count_residues_pdb,
    _rename_itp,
    run_martinize2,
)


SAMPLE_PDB = """\
TITLE     test
MODEL     1
ATOM      1  N   ALA A   1       0.000   0.000   0.000
ATOM      2  CA  ALA A   1       1.500   0.000   0.000
ATOM      3  N   GLY A   2       3.000   0.000   0.000
ATOM      4  CA  GLY A   2       4.500   0.000   0.000
ENDMDL
END
"""


class TestCountResidues:
    def test_count_unique_resid(self, tmp_path):
        p = tmp_path / "x.pdb"
        p.write_text(SAMPLE_PDB)
        assert _count_residues_pdb(p) == 2


class TestRenameItp:
    def test_renames_molecule_zero(self, tmp_path):
        target = tmp_path / "kgg.itp"
        src = tmp_path / "molecule_0.itp"
        src.write_text("[ moleculetype ]")
        _rename_itp(tmp_path, "kgg", target)
        assert target.exists()
        assert not src.exists()

    def test_renames_other_itp_when_no_molecule_zero(self, tmp_path):
        target = tmp_path / "x.itp"
        src = tmp_path / "Protein_X.itp"
        src.write_text("; stub")
        _rename_itp(tmp_path, "x", target)
        assert target.exists()
        assert not src.exists()

    def test_raises_when_no_itp(self, tmp_path):
        with pytest.raises(RuntimeError, match="No .itp"):
            _rename_itp(tmp_path, "x", tmp_path / "x.itp")

    def test_idempotent_overwrite(self, tmp_path):
        # Pre-existing target overwritten with martinize2's fresh output
        target = tmp_path / "kgg.itp"
        target.write_text("STALE")
        src = tmp_path / "molecule_0.itp"
        src.write_text("FRESH")
        _rename_itp(tmp_path, "kgg", target)
        assert target.read_text() == "FRESH"


class TestRunMartinize2:

    def _make_fake_run(self, tmp_path, cg_name, itp_name="molecule_0.itp"):
        cg_pdb = tmp_path / cg_name
        molecule_itp = tmp_path / itp_name

        def fake_run(cmd, **kw):
            cg_pdb.write_text("ATOM 1 BB ALA 1\n")
            molecule_itp.write_text("[ moleculetype ]\n")
            return MagicMock(returncode=0)

        return fake_run, cg_pdb, molecule_itp

    @patch("abmptools.cg.peptide.martinize_runner.run_command")
    def test_invokes_martinize2_with_m3_flags(self, mock_run, tmp_path):
        pdb = tmp_path / "in.pdb"
        pdb.write_text(SAMPLE_PDB)

        fake_run, _, _ = self._make_fake_run(tmp_path, "kgg_cg.pdb")
        mock_run.side_effect = fake_run

        itp_path, cg_path = run_martinize2(pdb, tmp_path, "kgg")

        cmd = mock_run.call_args.args[0]
        assert cmd[0] == "martinize2"
        assert "martini3001" in cmd
        # All-coil SS string of length n_res
        ss_idx = cmd.index("-ss")
        assert cmd[ss_idx + 1] == "CC"  # 2 residues
        assert "-maxwarn" in cmd
        assert "-elastic" not in cmd

        assert itp_path.name == "kgg.itp"
        assert cg_path.name == "kgg_cg.pdb"
        assert itp_path.exists()
        assert cg_path.exists()

    @patch("abmptools.cg.peptide.martinize_runner.run_command")
    def test_elastic_network_flag(self, mock_run, tmp_path):
        pdb = tmp_path / "in.pdb"
        pdb.write_text(SAMPLE_PDB)
        fake_run, _, _ = self._make_fake_run(tmp_path, "x_cg.pdb")
        mock_run.side_effect = fake_run

        run_martinize2(pdb, tmp_path, "x", elastic_network=True)
        cmd = mock_run.call_args.args[0]
        assert "-elastic" in cmd

    @patch("abmptools.cg.peptide.martinize_runner.run_command")
    def test_custom_martinize2_path(self, mock_run, tmp_path):
        pdb = tmp_path / "in.pdb"
        pdb.write_text(SAMPLE_PDB)
        fake_run, _, _ = self._make_fake_run(tmp_path, "x_cg.pdb")
        mock_run.side_effect = fake_run

        run_martinize2(pdb, tmp_path, "x", martinize2_path="/opt/martinize2")
        cmd = mock_run.call_args.args[0]
        assert cmd[0] == "/opt/martinize2"

    @patch("abmptools.cg.peptide.martinize_runner.run_command")
    def test_no_cg_pdb_raises(self, mock_run, tmp_path):
        pdb = tmp_path / "in.pdb"
        pdb.write_text(SAMPLE_PDB)

        # martinize2 succeeds (returncode=0) but doesn't produce cg_pdb
        # but does produce molecule_0.itp; we expect the cg_pdb missing
        # check to fire.
        molecule_itp = tmp_path / "molecule_0.itp"

        def fake_run(cmd, **kw):
            molecule_itp.write_text("[ moleculetype ]\n")
            return MagicMock(returncode=0)

        mock_run.side_effect = fake_run
        with pytest.raises(RuntimeError, match="martinize2 did not produce"):
            run_martinize2(pdb, tmp_path, "x")

    @patch("abmptools.cg.peptide.martinize_runner.run_command")
    def test_custom_maxwarn(self, mock_run, tmp_path):
        pdb = tmp_path / "in.pdb"
        pdb.write_text(SAMPLE_PDB)
        fake_run, _, _ = self._make_fake_run(tmp_path, "x_cg.pdb")
        mock_run.side_effect = fake_run

        run_martinize2(pdb, tmp_path, "x", maxwarn=42)
        cmd = mock_run.call_args.args[0]
        i = cmd.index("-maxwarn")
        assert cmd[i + 1] == "42"
