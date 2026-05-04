# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.cli (argparse subcommands)."""
from __future__ import annotations

import json
from io import StringIO
from pathlib import Path

import pytest

from abmptools.cg.membrane import cli


# ---------------------------------------------------------------------------
# example
# ---------------------------------------------------------------------------

def test_example_emits_valid_json(capsys):
    rc = cli.main(["example"])
    out = capsys.readouterr().out
    data = json.loads(out)
    assert rc == 0
    assert "lipids" in data
    assert "peptide" in data
    assert "umbrella" in data
    assert data["lipids"][0]["resname"] == "POPC"
    assert data["peptide"]["sequence"] == "KGG"


def test_example_includes_all_protocols(capsys):
    cli.main(["example"])
    data = json.loads(capsys.readouterr().out)
    assert "equilibration" in data
    assert "pulling" in data
    assert "umbrella" in data


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

def test_validate_returns_0_when_all_ok(capsys, tmp_path, monkeypatch):
    # Build a minimal valid config.
    cfg_path = tmp_path / "cfg.json"
    cli.main(["example"])
    capsys.readouterr()  # consume

    # Re-emit example using the API directly so we can write to file.
    from abmptools.cg.membrane.cli import _default_example
    cfg = _default_example()
    cfg.to_json(str(cfg_path))

    # Materialize 4 ITPs in ff_dir + monkeypatch shutil.which to "find" all tools.
    ff = tmp_path / "ff"
    ff.mkdir()
    for fname in [
        "martini_v3.0.0.itp",
        "martini_v3.0.0_solvents_v1.itp",
        "martini_v3.0.0_ions_v1.itp",
        "martini_v3.0.0_phospholipids_v1.itp",
    ]:
        (ff / fname).write_text("; stub\n")

    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        lambda x: f"/opt/bin/{x}",
    )

    rc = cli.main(["validate", "--config", str(cfg_path), "--ff-dir", str(ff)])
    assert rc == 0


def test_validate_returns_nonzero_when_itps_missing(capsys, tmp_path, monkeypatch):
    cfg_path = tmp_path / "cfg.json"
    from abmptools.cg.membrane.cli import _default_example
    cfg = _default_example()
    cfg.to_json(str(cfg_path))

    monkeypatch.setattr(
        "abmptools.cg.membrane.forcefield_check.shutil.which",
        lambda x: f"/opt/bin/{x}",
    )
    rc = cli.main([
        "validate", "--config", str(cfg_path),
        "--ff-dir", str(tmp_path / "no_ff"),
    ])
    assert rc != 0


# ---------------------------------------------------------------------------
# build (parser-level only; actual build is exercised in builder_mocked tests)
# ---------------------------------------------------------------------------

def test_build_parser_accepts_required_args(tmp_path):
    parser = cli._build_parser()
    args = parser.parse_args([
        "build", "--config", "cfg.json",
        "-o", "/tmp/out", "--ff-dir", "/tmp/ff",
    ])
    assert args.cmd == "build"
    assert args.config == "cfg.json"
    assert args.output_dir == "/tmp/out"
    assert args.ff_dir == "/tmp/ff"


def test_build_requires_config():
    parser = cli._build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["build"])


# ---------------------------------------------------------------------------
# make-windows (parser-level)
# ---------------------------------------------------------------------------

def test_make_windows_parser():
    parser = cli._build_parser()
    args = parser.parse_args([
        "make-windows",
        "--config", "cfg.json",
        "--pull-tpr", "pull/pull.tpr",
        "--pull-xtc", "pull/pull.xtc",
        "--pull-xvg", "pull/pullx.xvg",
        "--windows-dir", "windows",
    ])
    assert args.cmd == "make-windows"
    assert args.pull_tpr == "pull/pull.tpr"


# ---------------------------------------------------------------------------
# wham (parser-level)
# ---------------------------------------------------------------------------

def test_wham_parser():
    parser = cli._build_parser()
    args = parser.parse_args([
        "wham",
        "--config", "cfg.json",
        "--windows-dir", "windows",
        "--analysis-dir", "analysis",
        "--bootstrap-n", "100",
        "--temperature", "298.15",
    ])
    assert args.cmd == "wham"
    assert args.bootstrap_n == 100
    assert args.temperature == 298.15
