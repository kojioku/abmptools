# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.system_packer module."""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from abmptools.cg.peptide.system_packer import (
    add_ions,
    insert_peptides,
    make_index,
    solvate,
)


# ---------------------------------------------------------------------------
# insert_peptides
# ---------------------------------------------------------------------------

class TestInsertPeptides:
    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_single_species_cubic_box(self, mock_run, tmp_path):
        cg = tmp_path / "kgg_cg.pdb"
        cg.write_text("ATOM 1 BB ALA 1\n")

        # Create a chain of fake outputs that the real gmx would produce
        def fake_run(cmd, **kw):
            # find -o argument
            i = cmd.index("-o")
            (tmp_path / cmd[i + 1]).touch()
            return MagicMock(returncode=0)

        mock_run.side_effect = fake_run
        out = insert_peptides(
            [("kgg", cg, 5)], tmp_path, box_size_nm=8.0,
        )
        assert out.name == "packed.gro"
        cmd = mock_run.call_args.args[0]
        assert cmd[0] == "gmx"
        assert "insert-molecules" in cmd
        assert "-nmol" in cmd
        ix = cmd.index("-nmol")
        assert cmd[ix + 1] == "5"
        # Box args
        bi = cmd.index("-box")
        assert cmd[bi + 1] == "8.000"
        assert cmd[bi + 2] == "8.000"
        assert cmd[bi + 3] == "8.000"

    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_rectangular_box(self, mock_run, tmp_path):
        cg = tmp_path / "x_cg.pdb"
        cg.write_text("ATOM 1 BB ALA 1\n")

        def fake_run(cmd, **kw):
            (tmp_path / cmd[cmd.index("-o") + 1]).touch()
            return MagicMock(returncode=0)
        mock_run.side_effect = fake_run

        insert_peptides(
            [("x", cg, 1)], tmp_path,
            box_lengths_nm=[10.0, 12.0, 15.0],
        )
        cmd = mock_run.call_args.args[0]
        bi = cmd.index("-box")
        assert cmd[bi + 1] == "10.000"
        assert cmd[bi + 2] == "12.000"
        assert cmd[bi + 3] == "15.000"

    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_chained_two_species(self, mock_run, tmp_path):
        a = tmp_path / "a_cg.pdb"
        a.write_text("ATOM 1\n")
        b = tmp_path / "b_cg.pdb"
        b.write_text("ATOM 1\n")

        def fake_run(cmd, **kw):
            (tmp_path / cmd[cmd.index("-o") + 1]).touch()
            return MagicMock(returncode=0)
        mock_run.side_effect = fake_run

        insert_peptides(
            [("a", a, 3), ("b", b, 2)], tmp_path, box_size_nm=8.0,
        )
        # Two calls: first uses -box, second uses -f packed_00_a.gro
        assert mock_run.call_count == 2
        first_cmd = mock_run.call_args_list[0].args[0]
        second_cmd = mock_run.call_args_list[1].args[0]
        assert "-box" in first_cmd
        assert "-f" in second_cmd

    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_skips_zero_count(self, mock_run, tmp_path):
        a = tmp_path / "a.pdb"
        a.write_text("ATOM 1\n")
        b = tmp_path / "b.pdb"
        b.write_text("ATOM 1\n")

        def fake_run(cmd, **kw):
            (tmp_path / cmd[cmd.index("-o") + 1]).touch()
            return MagicMock(returncode=0)
        mock_run.side_effect = fake_run

        insert_peptides(
            [("a", a, 0), ("b", b, 3)], tmp_path, box_size_nm=8.0,
        )
        assert mock_run.call_count == 1  # only b inserted

    def test_no_box_specified_raises(self, tmp_path):
        a = tmp_path / "a.pdb"
        a.write_text("ATOM\n")
        with pytest.raises(ValueError, match="box"):
            insert_peptides([("a", a, 1)], tmp_path)

    def test_all_zero_counts_raises(self, tmp_path):
        a = tmp_path / "a.pdb"
        a.write_text("ATOM\n")
        with pytest.raises(ValueError, match="No peptides"):
            insert_peptides([("a", a, 0)], tmp_path, box_size_nm=8.0)

    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_seed_passed(self, mock_run, tmp_path):
        a = tmp_path / "a.pdb"
        a.write_text("ATOM\n")

        def fake_run(cmd, **kw):
            (tmp_path / cmd[cmd.index("-o") + 1]).touch()
            return MagicMock(returncode=0)
        mock_run.side_effect = fake_run

        insert_peptides(
            [("a", a, 1)], tmp_path, box_size_nm=8.0, seed=42,
        )
        cmd = mock_run.call_args.args[0]
        assert "-seed" in cmd
        assert cmd[cmd.index("-seed") + 1] == "42"


# ---------------------------------------------------------------------------
# solvate
# ---------------------------------------------------------------------------

class TestSolvate:
    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_invokes_gmx_solvate(self, mock_run, tmp_path):
        packed = tmp_path / "packed.gro"
        packed.write_text("title\n0\n   8 8 8\n")
        topol = tmp_path / "topol.top"
        topol.write_text("[ molecules ]\nKGG 1\n")
        water = tmp_path / "martini_v3.0.0_water.gro"
        water.write_text("water box\n")

        mock_run.return_value = MagicMock(returncode=0)
        out = solvate(packed, topol, tmp_path, water)
        assert out.name == "system_solv.gro"
        cmd = mock_run.call_args.args[0]
        assert cmd[0] == "gmx"
        assert cmd[1] == "solvate"
        assert "-cs" in cmd


# ---------------------------------------------------------------------------
# add_ions
# ---------------------------------------------------------------------------

class TestAddIons:
    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_grompp_then_genion(self, mock_run, tmp_path):
        solv = tmp_path / "solv.gro"
        solv.write_text("title\n")
        topol = tmp_path / "topol.top"
        topol.write_text("[ molecules ]\nKGG 1\n")
        mock_run.return_value = MagicMock(returncode=0)

        result = add_ions(solv, topol, tmp_path)

        # Two calls: grompp + genion
        assert mock_run.call_count == 2
        cmd1 = mock_run.call_args_list[0].args[0]
        cmd2 = mock_run.call_args_list[1].args[0]
        assert cmd1[1] == "grompp"
        assert cmd2[1] == "genion"
        # genion gets stdin "W\n"
        kwargs2 = mock_run.call_args_list[1].kwargs
        assert kwargs2.get("stdin_text") == "W\n"
        # genion has -neutral and -conc
        assert "-neutral" in cmd2
        ci = cmd2.index("-conc")
        assert cmd2[ci + 1] == "0.15"
        # ions.mdp written
        assert result.ions_mdp.exists()
        # ions.mdp contains expected MDP integrator
        assert "integrator" in result.ions_mdp.read_text()

    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_no_neutral_no_conc(self, mock_run, tmp_path):
        solv = tmp_path / "solv.gro"
        solv.write_text("title\n")
        topol = tmp_path / "topol.top"
        topol.write_text("[ molecules ]\n")
        mock_run.return_value = MagicMock(returncode=0)

        add_ions(solv, topol, tmp_path,
                 neutralize=False, conc_molar=0.0)
        cmd2 = mock_run.call_args_list[1].args[0]
        assert "-neutral" not in cmd2
        assert "-conc" not in cmd2


# ---------------------------------------------------------------------------
# make_index
# ---------------------------------------------------------------------------

class TestMakeIndex:
    @patch("abmptools.cg.peptide.system_packer.run_command")
    def test_invokes_make_ndx_with_quit_stdin(self, mock_run, tmp_path):
        gro = tmp_path / "x.gro"
        gro.write_text("title\n")
        mock_run.return_value = MagicMock(returncode=0)
        out = make_index(gro, tmp_path)
        assert out.name == "index.ndx"
        cmd = mock_run.call_args.args[0]
        assert cmd[1] == "make_ndx"
        kwargs = mock_run.call_args.kwargs
        assert kwargs["stdin_text"] == "q\n"
