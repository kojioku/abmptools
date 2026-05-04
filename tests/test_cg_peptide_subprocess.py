# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide._subprocess module."""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from abmptools.cg.peptide._subprocess import (
    CommandError,
    ensure_dir,
    run_command,
    write_text,
)


class TestRunCommand:
    """Tests for run_command()."""

    @patch("abmptools.cg.peptide._subprocess.subprocess.run")
    def test_success_returns_completed_process(self, mock_run):
        mock_run.return_value = MagicMock(
            returncode=0, stdout="ok", stderr=""
        )
        result = run_command(["echo", "hi"])
        assert result.returncode == 0
        mock_run.assert_called_once()

    @patch("abmptools.cg.peptide._subprocess.subprocess.run")
    def test_passes_cwd(self, mock_run, tmp_path):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        run_command(["ls"], cwd=tmp_path)
        kwargs = mock_run.call_args.kwargs
        assert kwargs["cwd"] == tmp_path

    @patch("abmptools.cg.peptide._subprocess.subprocess.run")
    def test_passes_stdin(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        run_command(["cat"], stdin_text="hello\n")
        kwargs = mock_run.call_args.kwargs
        assert kwargs["input"] == "hello\n"

    @patch("abmptools.cg.peptide._subprocess.subprocess.run")
    def test_failure_raises_command_error(self, mock_run):
        mock_run.return_value = MagicMock(
            returncode=2, stdout="", stderr="boom"
        )
        with pytest.raises(CommandError) as excinfo:
            run_command(["false"])
        assert excinfo.value.returncode == 2
        assert "boom" in excinfo.value.stderr
        assert "false" in str(excinfo.value)

    @patch("abmptools.cg.peptide._subprocess.subprocess.run")
    def test_check_false_does_not_raise(self, mock_run):
        mock_run.return_value = MagicMock(returncode=99, stderr="oops")
        result = run_command(["false"], check=False)
        assert result.returncode == 99


class TestEnsureDir:
    def test_creates_dir(self, tmp_path):
        target = tmp_path / "a" / "b" / "c"
        result = ensure_dir(target)
        assert result.exists() and result.is_dir()

    def test_existing_dir_idempotent(self, tmp_path):
        ensure_dir(tmp_path)
        # second call must not raise
        ensure_dir(tmp_path)


class TestWriteText:
    def test_writes_text(self, tmp_path):
        target = tmp_path / "x" / "y" / "out.txt"
        write_text(target, "hello\n")
        assert target.read_text() == "hello\n"

    def test_creates_parent_dirs(self, tmp_path):
        target = tmp_path / "deep" / "deeper" / "f.txt"
        write_text(target, "x")
        assert target.parent.exists()
