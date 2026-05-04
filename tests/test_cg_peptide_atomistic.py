# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.peptide_atomistic module."""
from __future__ import annotations

import logging
from unittest.mock import MagicMock, patch

from abmptools.cg.peptide.peptide_atomistic import (
    ARTIFACT_PRONE_RESIDUES,
    _atoms_to_pdb,
    _build_extended_backbone,
    build_atomistic_pdb,
)


class TestExtendedBackbone:
    def test_kgg_atom_set(self):
        # K (5: N CA C O CB), G (4), G (4) => 13 atoms
        atoms = _build_extended_backbone("KGG")
        names = [a["name"] for a in atoms]
        assert names == [
            "N", "CA", "C", "O", "CB",
            "N", "CA", "C", "O",
            "N", "CA", "C", "O",
        ]

    def test_resids_increment(self):
        atoms = _build_extended_backbone("AAA")
        resids = sorted(set(a["resid"] for a in atoms))
        assert resids == [1, 2, 3]

    def test_resname_3letter(self):
        atoms = _build_extended_backbone("A")
        assert all(a["resname"] == "ALA" for a in atoms)

    def test_pdb_format(self):
        atoms = _build_extended_backbone("A")
        pdb = _atoms_to_pdb(atoms, title="X")
        assert "ATOM" in pdb
        assert "MODEL" in pdb
        assert "END" in pdb
        assert "ALA" in pdb


class TestBuildAtomisticPdb:
    """Tests for the public build_atomistic_pdb() entry point."""

    @patch("abmptools.cg.peptide.peptide_atomistic.run_command")
    @patch("abmptools.cg.peptide.peptide_atomistic.shutil.which")
    def test_uses_tleap_when_available(
        self, mock_which, mock_run, tmp_path,
    ):
        mock_which.return_value = "/u/bin/tleap"
        pdb_path = tmp_path / "kgg_atomistic.pdb"

        def fake_run(cmd, **kw):
            # Simulate tleap writing the PDB
            pdb_path.write_text("ATOM      1  N   ALA   1   0 0 0\n")
            return MagicMock(returncode=0)

        mock_run.side_effect = fake_run

        result = build_atomistic_pdb("KGG", "kgg", tmp_path)
        assert result == pdb_path
        assert pdb_path.exists()
        # tleap CLI was invoked
        cmd = mock_run.call_args.args[0]
        assert cmd[0] == "tleap"
        assert "-f" in cmd

    @patch("abmptools.cg.peptide.peptide_atomistic.shutil.which")
    def test_falls_back_when_tleap_missing(
        self, mock_which, tmp_path, caplog,
    ):
        mock_which.return_value = None
        with caplog.at_level(logging.WARNING,
                             logger="abmptools.cg.peptide.peptide_atomistic"):
            result = build_atomistic_pdb("KGG", "kgg", tmp_path)
        assert result.exists()
        assert "fallback" in caplog.text.lower() \
            or "extended" in caplog.text.lower()

    @patch("abmptools.cg.peptide.peptide_atomistic.shutil.which")
    def test_artifact_warning_for_aromatic(
        self, mock_which, tmp_path, caplog,
    ):
        mock_which.return_value = None  # force fallback
        with caplog.at_level(logging.WARNING,
                             logger="abmptools.cg.peptide.peptide_atomistic"):
            build_atomistic_pdb("KGW", "kgw", tmp_path)
        # NaN warning fires for W
        assert "NaN" in caplog.text

    @patch("abmptools.cg.peptide.peptide_atomistic.shutil.which")
    def test_no_artifact_warning_for_simple(
        self, mock_which, tmp_path, caplog,
    ):
        mock_which.return_value = None
        with caplog.at_level(logging.WARNING,
                             logger="abmptools.cg.peptide.peptide_atomistic"):
            build_atomistic_pdb("AAGG", "ag", tmp_path)
        # No W/F/Y -> no NaN warning
        assert "NaN" not in caplog.text

    @patch("abmptools.cg.peptide.peptide_atomistic.shutil.which")
    def test_prefer_tleap_false_uses_fallback(
        self, mock_which, tmp_path,
    ):
        mock_which.return_value = "/u/bin/tleap"  # tleap present but ignored
        result = build_atomistic_pdb(
            "AAA", "aaa", tmp_path, prefer_tleap=False,
        )
        assert result.exists()
        with open(result) as f:
            content = f.read()
        # Extended backbone path produces TITLE line
        assert "TITLE" in content

    def test_artifact_residue_set(self):
        # Aromatic only — Plan で確定した範囲
        assert ARTIFACT_PRONE_RESIDUES == set("WFY")
