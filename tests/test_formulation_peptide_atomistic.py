# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.peptide_atomistic."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.formulation import peptide_atomistic as pa


def test_sequence_to_tleap_residues_natural():
    out = pa.sequence_to_tleap_residues("AAAAA")
    assert out == "{ ALA ALA ALA ALA ALA }"


def test_sequence_to_tleap_residues_with_caps():
    out = pa.sequence_to_tleap_residues("AAA", cap_n="ACE", cap_c="NME")
    assert out == "{ ACE ALA ALA ALA NME }"


def test_sequence_to_tleap_residues_rejects_non_natural():
    with pytest.raises(ValueError, match="non-natural"):
        pa.sequence_to_tleap_residues("AXAA")


def test_render_tleap_input_basic():
    text = pa.render_tleap_input_from_sequence(
        sequence="AAAAA", out_pdb="/tmp/x.pdb",
    )
    assert "source leaprc.protein.ff14SB" in text
    assert "source leaprc.water.tip3p" in text
    assert "pep = sequence { ALA ALA ALA ALA ALA }" in text
    assert "savepdb pep /tmp/x.pdb" in text
    assert "quit" in text


def test_render_tleap_input_with_bond_directives():
    text = pa.render_tleap_input_from_sequence(
        sequence="CCC", out_pdb="/tmp/x.pdb",
        extra_bonds=[("1.SG", "3.SG")],
    )
    assert "bond pep.1.SG pep.3.SG" in text


def test_build_peptide_from_sequence_mocked_tleap(tmp_path):
    expected_pdb = tmp_path / "out" / "pep.pdb"

    def fake_run(argv, **kwargs):
        # Simulate tleap creating the saved PDB
        expected_pdb.parent.mkdir(parents=True, exist_ok=True)
        expected_pdb.write_text("HETATM    1  N   ALA A   1       0.000   0.000   0.000\nEND\n")
        from subprocess import CompletedProcess
        return CompletedProcess(argv, returncode=0, stdout="ok", stderr="")

    with patch("abmptools.formulation.peptide_atomistic.run_command", side_effect=fake_run):
        out = pa.build_peptide_from_sequence(
            sequence="AAAAA",
            out_pdb=str(expected_pdb),
            workdir=str(tmp_path / "wd"),
        )
        assert Path(out).is_file()


def test_build_peptide_from_sequence_tleap_failure(tmp_path):
    def fake_run(argv, **kwargs):
        from subprocess import CompletedProcess
        return CompletedProcess(argv, returncode=1, stdout="", stderr="tleap error")

    with patch("abmptools.formulation.peptide_atomistic.run_command", side_effect=fake_run):
        with pytest.raises(pa.CommandError):
            pa.build_peptide_from_sequence(
                sequence="AAAAA", out_pdb=str(tmp_path / "x.pdb"),
                workdir=str(tmp_path / "wd"),
            )


def test_normalize_existing_pdb_strips_cryst_and_waters(tmp_path):
    src = tmp_path / "raw.pdb"
    src.write_text(
        "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n"
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000\n"
        "HETATM    3  O   HOH B   1       5.000   5.000   5.000\n"
        "CONECT    1    2\n"
        "END\n"
    )
    out = tmp_path / "clean.pdb"
    pa.normalize_existing_pdb(src_pdb=str(src), out_pdb=str(out))
    text = out.read_text()
    assert "CRYST1" not in text
    assert "CONECT" not in text
    assert "HOH" not in text
    assert "ALA" in text
