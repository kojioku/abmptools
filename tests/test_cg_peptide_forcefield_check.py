# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.forcefield_check module."""
from __future__ import annotations

from unittest.mock import patch

from abmptools.cg.peptide.forcefield_check import (
    OPTIONAL_MARTINI_FILES,
    REQUIRED_MARTINI_FILES,
    check_external_tools,
    check_martini_files,
    report,
)


class TestCheckExternalTools:
    @patch("abmptools.cg.peptide.forcefield_check.shutil.which")
    def test_all_present(self, mock_which):
        mock_which.side_effect = lambda x: f"/usr/bin/{x}"
        tools = check_external_tools()
        assert {t.name for t in tools} == {"martinize2", "gmx", "tleap"}
        assert all(t.found for t in tools)

    @patch("abmptools.cg.peptide.forcefield_check.shutil.which")
    def test_tleap_is_optional(self, mock_which):
        mock_which.side_effect = lambda x: None if x == "tleap" else f"/u/{x}"
        tools = check_external_tools()
        tleap = next(t for t in tools if t.name == "tleap")
        assert tleap.found is False
        assert tleap.required is False

    @patch("abmptools.cg.peptide.forcefield_check.shutil.which")
    def test_martinize2_required_missing(self, mock_which):
        mock_which.side_effect = (
            lambda x: None if x == "martinize2" else f"/u/{x}"
        )
        tools = check_external_tools()
        m2 = next(t for t in tools if t.name == "martinize2")
        assert m2.found is False
        assert m2.required is True

    @patch("abmptools.cg.peptide.forcefield_check.shutil.which")
    def test_custom_paths(self, mock_which):
        mock_which.side_effect = lambda x: f"/custom/{x}"
        tools = check_external_tools(
            martinize2="/opt/martinize2",
            gmx="/opt/gmx",
            tleap="/opt/tleap",
        )
        # which() received the custom string
        called_args = [c.args[0] for c in mock_which.call_args_list]
        assert "/opt/martinize2" in called_args
        assert "/opt/gmx" in called_args
        assert "/opt/tleap" in called_args
        # All resolved (mock returns truthy)
        assert all(t.found for t in tools)


class TestCheckMartiniFiles:
    def test_all_required_present(self, tmp_path):
        for fname in REQUIRED_MARTINI_FILES:
            (tmp_path / fname).write_text("; stub")
        files = check_martini_files(str(tmp_path))
        required = [f for f in files if not f.optional]
        assert all(f.found for f in required)
        assert {f.name for f in required} == set(REQUIRED_MARTINI_FILES)

    def test_optional_water_gro_listed_separately(self, tmp_path):
        files = check_martini_files(str(tmp_path))
        optional = [f for f in files if f.optional]
        assert {f.name for f in optional} == set(OPTIONAL_MARTINI_FILES)
        assert all(not f.found for f in optional)  # not placed

    def test_optional_can_be_found(self, tmp_path):
        for fname in OPTIONAL_MARTINI_FILES:
            (tmp_path / fname).write_text("; stub")
        files = check_martini_files(str(tmp_path))
        optional = [f for f in files if f.optional]
        assert all(f.found for f in optional)

    def test_some_files_missing(self, tmp_path):
        (tmp_path / REQUIRED_MARTINI_FILES[0]).write_text("; stub")
        files = check_martini_files(str(tmp_path))
        # First required found
        assert files[0].found is True
        # Others (required + optional) missing
        for f in files[1:]:
            assert f.found is False

    def test_empty_dir_returns_all_missing(self):
        files = check_martini_files("")
        assert all(not f.found for f in files)
        assert len(files) == (
            len(REQUIRED_MARTINI_FILES) + len(OPTIONAL_MARTINI_FILES)
        )


class TestReport:
    @patch("abmptools.cg.peptide.forcefield_check.shutil.which")
    def test_buildable_when_required_ok(self, mock_which, tmp_path, capsys):
        mock_which.side_effect = lambda x: f"/u/{x}"
        # Place only REQUIRED files (water.gro is optional)
        for fname in REQUIRED_MARTINI_FILES:
            (tmp_path / fname).write_text("; stub")
        tools = check_external_tools()
        files = check_martini_files(str(tmp_path))
        ok = report(tools, files, str(tmp_path))
        assert ok is True
        captured = capsys.readouterr().out
        assert "[OK]" in captured
        # Optional water.gro shown but doesn't block build
        assert "[opt]" in captured

    @patch("abmptools.cg.peptide.forcefield_check.shutil.which")
    def test_not_buildable_when_files_missing(
        self, mock_which, tmp_path, capsys,
    ):
        mock_which.side_effect = lambda x: f"/u/{x}"
        tools = check_external_tools()
        files = check_martini_files(str(tmp_path))  # empty dir
        ok = report(tools, files, str(tmp_path))
        assert ok is False
        captured = capsys.readouterr().out
        assert "MISSING" in captured
        assert "Download" in captured

    @patch("abmptools.cg.peptide.forcefield_check.shutil.which")
    def test_not_buildable_when_required_tool_missing(
        self, mock_which, tmp_path, capsys,
    ):
        mock_which.side_effect = (
            lambda x: None if x == "martinize2" else f"/u/{x}"
        )
        for fname in REQUIRED_MARTINI_FILES:
            (tmp_path / fname).write_text("; stub")
        tools = check_external_tools()
        files = check_martini_files(str(tmp_path))
        ok = report(tools, files, str(tmp_path))
        assert ok is False
        captured = capsys.readouterr().out
        assert "[MISS]" in captured
