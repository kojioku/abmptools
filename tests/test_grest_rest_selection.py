# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.rest_selection."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.genesis.grest.models import RESTSelectionSpec
from abmptools.genesis.grest.rest_selection import (
    RESTSelectionResult,
    format_genesis_selection,
    parse_cpptraj_resinfo,
    parse_explicit_residues,
    render_cpptraj_around_script,
    resolve_around,
    resolve_rest_selection,
)


# ---------------------------------------------------------------------------
# parse_explicit_residues
# ---------------------------------------------------------------------------

class TestParseExplicitResidues:
    def test_single_residue(self):
        assert parse_explicit_residues(["21"]) == [21]

    def test_hyphen_range(self):
        assert parse_explicit_residues(["1-5"]) == [1, 2, 3, 4, 5]

    def test_comma_list(self):
        assert parse_explicit_residues(["21,96,274"]) == [21, 96, 274]

    def test_mixed_specs(self):
        assert parse_explicit_residues(["1-3,5,7-9"]) == [1, 2, 3, 5, 7, 8, 9]

    def test_multiple_string_entries(self):
        assert parse_explicit_residues(["1-3", "5", "7-9"]) == [1, 2, 3, 5, 7, 8, 9]

    def test_dedup_and_sort(self):
        # Overlap + repeat across entries.
        assert parse_explicit_residues(["3,1", "1-2"]) == [1, 2, 3]

    def test_whitespace_tolerant(self):
        assert parse_explicit_residues(["  1 - 3 , 5 "]) == [1, 2, 3, 5]

    def test_empty_token_skipped(self):
        # Trailing comma should not raise.
        assert parse_explicit_residues(["1-3,"]) == [1, 2, 3]

    def test_malformed_raises(self):
        with pytest.raises(ValueError, match="Cannot parse"):
            parse_explicit_residues(["1-x"])

    def test_descending_range_raises(self):
        with pytest.raises(ValueError, match="start > end"):
            parse_explicit_residues(["5-1"])

    def test_pure_python_138_residues(self):
        # POC scale: 138 residues from a single-chain protein.
        result = parse_explicit_residues(["1-138"])
        assert len(result) == 138
        assert result[0] == 1
        assert result[-1] == 138


# ---------------------------------------------------------------------------
# format_genesis_selection
# ---------------------------------------------------------------------------

class TestFormatGenesisSelection:
    def test_empty(self):
        assert format_genesis_selection([]) == ""

    def test_single(self):
        assert format_genesis_selection([21]) == "rno:21"

    def test_contiguous_run_compacted(self):
        assert format_genesis_selection([1, 2, 3]) == "rno:1-3"

    def test_mixed_compacted(self):
        # 21 alone, 96 alone, 274-275 contiguous.
        assert format_genesis_selection([21, 96, 274, 275]) == "rno:21,96,274-275"

    def test_unsorted_input_normalised(self):
        assert format_genesis_selection([3, 1, 2]) == "rno:1-3"

    def test_dedup(self):
        assert format_genesis_selection([1, 1, 2, 2, 3]) == "rno:1-3"


# ---------------------------------------------------------------------------
# render_cpptraj_around_script
# ---------------------------------------------------------------------------

class TestRenderCpptrajAroundScript:
    def test_rno_prefix_normalised(self, tmp_path):
        prmtop = tmp_path / "x.prmtop"
        out_txt = tmp_path / "rest.txt"
        script = render_cpptraj_around_script(
            prmtop=prmtop,
            center="rno:96",
            radius_A=5.0,
            out_txt=out_txt,
        )
        assert "parm" in script
        # GENESIS-style "rno:96" should be normalised to cpptraj ":96".
        assert "(:96)<@5.0" in script
        # Solvent / counterion exclusion.
        assert ":WAT" in script
        assert ":Na+" in script
        assert "run" in script
        assert "quit" in script

    def test_colon_prefix_passthrough(self, tmp_path):
        script = render_cpptraj_around_script(
            prmtop=tmp_path / "x.prmtop",
            center=":42",
            radius_A=4.5,
            out_txt=tmp_path / "out.txt",
        )
        assert "(:42)<@4.5" in script

    def test_bare_resno(self, tmp_path):
        script = render_cpptraj_around_script(
            prmtop=tmp_path / "x.prmtop",
            center="42",
            radius_A=3.0,
            out_txt=tmp_path / "out.txt",
        )
        assert "(:42)<@3.0" in script


# ---------------------------------------------------------------------------
# parse_cpptraj_resinfo
# ---------------------------------------------------------------------------

class TestParseCpptrajResinfo:
    def test_basic_dump(self):
        # Excerpt-shaped lines as written by cpptraj's resinfo.
        text = """\
#Res  Name First  Last Natom #Orig #Mol C I
   1 LYS      1    24    24     1     1
   5 ARG     68    91    24     5     1
  96 LYS   1418  1439    22    96     1
"""
        result = parse_cpptraj_resinfo(text)
        assert result == [1, 5, 96]

    def test_skips_non_data_lines(self):
        text = """\
some preamble
#Res  Name First  Last Natom #Orig
  21 ARG    100   123    23    21     1

other noise
  96 LYS  1000  1020    20    96     1
"""
        assert parse_cpptraj_resinfo(text) == [21, 96]

    def test_empty_returns_empty(self):
        assert parse_cpptraj_resinfo("") == []
        assert parse_cpptraj_resinfo("# only a comment\n") == []

    def test_dedup_and_sort(self):
        text = """\
   5 ARG  100  120  20   5  1
   1 LYS    1   24  24   1  1
   5 ARG  100  120  20   5  1
"""
        assert parse_cpptraj_resinfo(text) == [1, 5]


# ---------------------------------------------------------------------------
# resolve_around (mocked subprocess + file write)
# ---------------------------------------------------------------------------

class TestResolveAround:
    def test_explicit_mode_rejected(self, tmp_path):
        spec = RESTSelectionSpec(mode="explicit", residues=["1"])
        with pytest.raises(ValueError, match="mode='around'"):
            resolve_around(spec, tmp_path / "x.prmtop", tmp_path)

    def test_invokes_cpptraj_and_parses(self, tmp_path):
        spec = RESTSelectionSpec(
            mode="around", center="rno:96", radius_A=5.0
        )
        prmtop = tmp_path / "x.prmtop"
        prmtop.write_text("# fake prmtop")
        workdir = tmp_path / "rest_workdir"

        # Fake cpptraj behaviour: write resinfo to the expected out_txt
        # via side_effect, then return success.
        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            # cpptraj writes its output to out_txt declared inside the .cpptraj
            # script. Mirror that here.
            out_txt = Path(workdir) / "rest_residues.txt"
            out_txt.write_text(
                "#Res  Name First  Last Natom #Orig\n"
                "  21 ARG  100  123  23  21  1\n"
                "  96 LYS 1418 1439  22  96  1\n"
                " 274 ASP 4029 4040  12 274  1\n"
            )
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "", "")

        with patch(
            "abmptools.genesis.grest.rest_selection.run_command",
            side_effect=fake_run,
        ):
            residues = resolve_around(
                spec, prmtop, workdir, cpptraj_path="cpptraj"
            )

        assert residues == [21, 96, 274]
        # cpptraj input script was written.
        assert (workdir / "rest_around.cpptraj").exists()
        script = (workdir / "rest_around.cpptraj").read_text()
        assert "(:96)<@5.0" in script
        assert str(prmtop) in script


# ---------------------------------------------------------------------------
# resolve_rest_selection (top-level)
# ---------------------------------------------------------------------------

class TestResolveRESTSelection:
    def test_explicit_mode_no_external_tools(self):
        spec = RESTSelectionSpec(
            mode="explicit", residues=["1-3", "5"]
        )
        result = resolve_rest_selection(spec)
        assert isinstance(result, RESTSelectionResult)
        assert result.residues == [1, 2, 3, 5]
        assert result.n_residues == 4
        assert result.selection_string == "rno:1-3,5"

    def test_around_mode_requires_prmtop(self):
        spec = RESTSelectionSpec(
            mode="around", center="rno:96", radius_A=5.0
        )
        with pytest.raises(ValueError, match="prmtop and workdir"):
            resolve_rest_selection(spec)

    def test_around_mode_dispatches(self, tmp_path):
        spec = RESTSelectionSpec(
            mode="around", center="rno:96", radius_A=5.0
        )
        prmtop = tmp_path / "x.prmtop"
        prmtop.write_text("# fake prmtop")
        workdir = tmp_path / "wd"

        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            out_txt = Path(workdir) / "rest_residues.txt"
            out_txt.write_text(
                "  96 LYS 1418 1439 22 96 1\n"
                " 274 ASP 4029 4040 12 274 1\n"
            )
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "", "")

        with patch(
            "abmptools.genesis.grest.rest_selection.run_command",
            side_effect=fake_run,
        ):
            result = resolve_rest_selection(
                spec, prmtop=prmtop, workdir=workdir
            )

        assert result.residues == [96, 274]
        assert result.selection_string == "rno:96,274"
