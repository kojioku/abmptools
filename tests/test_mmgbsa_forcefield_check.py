# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.mmgbsa.forcefield_check."""
from __future__ import annotations

import io
from contextlib import redirect_stdout
from unittest.mock import patch

import pytest

from abmptools.genesis.mmgbsa.forcefield_check import (
    check_external_tools,
    check_python_modules,
    report,
)
from abmptools.genesis.mmgbsa.models import (
    MMGBSABuildConfig,
    TargetSpec,
)


# ---------------------------------------------------------------------------
# check_external_tools
# ---------------------------------------------------------------------------

class TestCheckExternalTools:
    def test_returns_four_tools(self):
        tools = check_external_tools()
        names = [t.name for t in tools]
        assert names == ["atdyn", "tleap", "acpype", "mpirun"]

    def test_all_required(self):
        tools = check_external_tools()
        for t in tools:
            assert t.required is True

    def test_resolves_via_shutil_which(self):
        with patch(
            "abmptools.genesis.mmgbsa.forcefield_check.shutil.which",
            return_value="/fake/bin/atdyn",
        ):
            tools = check_external_tools()
        for t in tools:
            assert t.found is True
            assert t.path == "/fake/bin/atdyn"

    def test_missing_tool_recorded(self):
        with patch(
            "abmptools.genesis.mmgbsa.forcefield_check.shutil.which",
            return_value=None,
        ):
            tools = check_external_tools()
        for t in tools:
            assert t.found is False
            assert t.path is None


# ---------------------------------------------------------------------------
# check_python_modules
# ---------------------------------------------------------------------------

class TestCheckPythonModules:
    def test_returns_two_modules(self):
        py = check_python_modules()
        names = [m.name for m in py]
        assert names == ["Biopython", "matplotlib"]

    def test_both_required(self):
        for m in check_python_modules():
            assert m.required is True

    def test_real_environment(self):
        # In the abmptoolsenv test environment, both should be importable.
        py = check_python_modules()
        for m in py:
            # If Biopython / matplotlib happen to be installed, found=True;
            # if not, they show up as MISS. Just verify shape, not value.
            assert isinstance(m.found, bool)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def _valid_cfg() -> MMGBSABuildConfig:
    return MMGBSABuildConfig(
        targets=[
            TargetSpec(pdb="3772L.pdb", ligand_resno=201),
            TargetSpec(pdb="9MM52.pdb", ligand_resno=201, chain="A"),
        ],
        project_name="poc_run",
    )


class TestReport:
    def test_ok_when_all_present(self):
        cfg = _valid_cfg()
        with patch(
            "abmptools.genesis.mmgbsa.forcefield_check.shutil.which",
            return_value="/fake/bin/x",
        ), patch(
            "abmptools.genesis.mmgbsa.forcefield_check.importlib.import_module",
        ) as m_import:
            m_import.return_value = object()
            buf = io.StringIO()
            with redirect_stdout(buf):
                ok = report(cfg)
        out = buf.getvalue()
        assert ok is True
        assert "Configuration: OK" in out
        assert "[OK]    atdyn" in out
        assert "[OK]    Biopython" in out
        assert "Targets preview" in out
        assert "T01: 3772L.pdb  ligand_resno=201" in out

    def test_fail_when_tool_missing(self):
        cfg = _valid_cfg()
        with patch(
            "abmptools.genesis.mmgbsa.forcefield_check.shutil.which",
            return_value=None,  # all tools missing
        ), patch(
            "abmptools.genesis.mmgbsa.forcefield_check.importlib.import_module",
        ) as m_import:
            m_import.return_value = object()
            buf = io.StringIO()
            with redirect_stdout(buf):
                ok = report(cfg)
        assert ok is False
        out = buf.getvalue()
        assert "[MISS]  atdyn" in out
        assert "GENESIS build instructions" in out

    def test_fail_when_python_module_missing(self):
        cfg = _valid_cfg()
        with patch(
            "abmptools.genesis.mmgbsa.forcefield_check.shutil.which",
            return_value="/fake/bin/x",
        ), patch(
            "abmptools.genesis.mmgbsa.forcefield_check.importlib.import_module",
            side_effect=ImportError("missing"),
        ):
            buf = io.StringIO()
            with redirect_stdout(buf):
                ok = report(cfg)
        assert ok is False
        out = buf.getvalue()
        assert "[MISS]  Biopython" in out
        assert "pip install abmptools[mmgbsa]" in out

    def test_folder_mode_preview(self):
        cfg = MMGBSABuildConfig(input_dir="./input", project_name="folder")
        with patch(
            "abmptools.genesis.mmgbsa.forcefield_check.shutil.which",
            return_value="/fake/bin/x",
        ), patch(
            "abmptools.genesis.mmgbsa.forcefield_check.importlib.import_module",
        ) as m_import:
            m_import.return_value = object()
            buf = io.StringIO()
            with redirect_stdout(buf):
                report(cfg)
        out = buf.getvalue()
        assert "folder mode: ./input/*.pdb" in out

    def test_chain_preview(self):
        cfg = _valid_cfg()
        buf = io.StringIO()
        with patch(
            "abmptools.genesis.mmgbsa.forcefield_check.shutil.which",
            return_value="/fake/bin/x",
        ), patch(
            "abmptools.genesis.mmgbsa.forcefield_check.importlib.import_module",
        ) as m_import:
            m_import.return_value = object()
            with redirect_stdout(buf):
                report(cfg)
        out = buf.getvalue()
        assert "chain A" in out  # second target has chain="A"
        assert "any chain" in out  # first target has chain=None
