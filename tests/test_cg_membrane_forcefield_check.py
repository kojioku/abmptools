# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.forcefield_check."""
from __future__ import annotations

from abmptools.cg.membrane import forcefield_check


# ---------------------------------------------------------------------------
# REQUIRED_MARTINI_FILES
# ---------------------------------------------------------------------------

def test_required_files_count_is_4():
    assert len(forcefield_check.REQUIRED_MARTINI_FILES) == 4


def test_required_files_includes_phospholipids():
    assert (
        "martini_v3.0.0_phospholipids_v1.itp"
        in forcefield_check.REQUIRED_MARTINI_FILES
    )


def test_required_files_includes_main_three():
    assert "martini_v3.0.0.itp" in forcefield_check.REQUIRED_MARTINI_FILES
    assert (
        "martini_v3.0.0_solvents_v1.itp"
        in forcefield_check.REQUIRED_MARTINI_FILES
    )
    assert (
        "martini_v3.0.0_ions_v1.itp"
        in forcefield_check.REQUIRED_MARTINI_FILES
    )


# ---------------------------------------------------------------------------
# check_martini_files
# ---------------------------------------------------------------------------

def test_all_present(tmp_path):
    for fname in forcefield_check.REQUIRED_MARTINI_FILES:
        (tmp_path / fname).write_text("; stub itp\n")
    statuses = forcefield_check.check_martini_files(str(tmp_path))
    assert len(statuses) == 4
    assert all(s.found for s in statuses)
    assert all(not s.optional for s in statuses)


def test_missing_phospholipids_only(tmp_path):
    for fname in forcefield_check.REQUIRED_MARTINI_FILES:
        if fname == forcefield_check.PHOSPHOLIPIDS_ITP:
            continue
        (tmp_path / fname).write_text("; stub itp\n")
    statuses = forcefield_check.check_martini_files(str(tmp_path))
    found = {s.name: s.found for s in statuses}
    assert found["martini_v3.0.0.itp"] is True
    assert found[forcefield_check.PHOSPHOLIPIDS_ITP] is False


def test_empty_dir_string_yields_all_missing():
    statuses = forcefield_check.check_martini_files("")
    assert all(not s.found for s in statuses)


# ---------------------------------------------------------------------------
# check_external_tools
# ---------------------------------------------------------------------------

def test_check_external_tools_reports_insane_required(monkeypatch):
    # Pretend nothing is on PATH
    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        lambda x: None,
    )
    tools = forcefield_check.check_external_tools()
    by_name = {t.name: t for t in tools}
    assert "insane" in by_name
    assert by_name["insane"].required is True
    assert by_name["insane"].found is False


def test_check_external_tools_tleap_optional(monkeypatch):
    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        lambda x: None,
    )
    tools = forcefield_check.check_external_tools()
    tleap = next(t for t in tools if t.name == "tleap")
    assert tleap.required is False


def test_check_external_tools_finds_present_tool(monkeypatch):
    fake_paths = {
        "insane": "/opt/bin/insane",
        "martinize2": "/opt/bin/martinize2",
        "gmx": "/opt/bin/gmx",
        "tleap": "/opt/bin/tleap",
    }
    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        lambda x: fake_paths.get(x),
    )
    tools = forcefield_check.check_external_tools()
    by_name = {t.name: t for t in tools}
    for name, path in fake_paths.items():
        assert by_name[name].found is True
        assert by_name[name].path == path


def test_check_external_tools_custom_executable_paths(monkeypatch):
    seen = []

    def fake_which(arg):
        seen.append(arg)
        return None

    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        fake_which,
    )
    forcefield_check.check_external_tools(
        insane="/tmp/insane",
        gmx="/tmp/gmx",
    )
    assert "/tmp/insane" in seen
    assert "/tmp/gmx" in seen


# ---------------------------------------------------------------------------
# report (smoke; just ensure it doesn't blow up and returns bool)
# ---------------------------------------------------------------------------

def test_report_all_ok_returns_true(tmp_path, capsys, monkeypatch):
    for fname in forcefield_check.REQUIRED_MARTINI_FILES:
        (tmp_path / fname).write_text("; stub itp\n")
    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        lambda x: f"/opt/bin/{x}",
    )
    tools = forcefield_check.check_external_tools()
    files = forcefield_check.check_martini_files(str(tmp_path))
    result = forcefield_check.report(tools, files, str(tmp_path))
    assert result is True


def test_report_missing_files_returns_false(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        lambda x: f"/opt/bin/{x}",
    )
    tools = forcefield_check.check_external_tools()
    files = forcefield_check.check_martini_files(str(tmp_path))   # empty dir
    result = forcefield_check.report(tools, files, str(tmp_path))
    captured = capsys.readouterr()
    assert result is False
    assert "MISSING" in captured.out
    assert "Download" in captured.out
