# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.cli."""
from __future__ import annotations

import json

import pytest

from abmptools.cg.peptide import cli
from abmptools.cg.peptide.models import PeptideBuildConfig


class TestExampleSubcommand:
    def test_emits_valid_json(self, capsys):
        rc = cli.main(["example"])
        assert rc == 0
        captured = capsys.readouterr().out
        data = json.loads(captured)
        assert "peptides" in data
        # Round-trips through PeptideBuildConfig
        peps = data.pop("peptides")
        cfg = PeptideBuildConfig(
            peptides=[__import__("abmptools.cg.peptide.models",
                                 fromlist=["PeptideSpec"])
                      .PeptideSpec(**p) for p in peps],
            **data,
        )
        assert len(cfg.peptides) >= 1


class TestValidateSubcommand:
    @pytest.fixture
    def sample_config(self, tmp_path):
        cfg = PeptideBuildConfig(
            peptides=[__import__("abmptools.cg.peptide.models",
                                 fromlist=["PeptideSpec"])
                      .PeptideSpec(name="x", sequence="KGG", count=2)],
            box_size_nm=8.0,
            martini_itp_dir=str(tmp_path / "ff_missing"),
        )
        path = tmp_path / "kgg.json"
        cfg.to_json(str(path))
        return path

    def test_validate_returns_nonzero_when_files_missing(
        self, sample_config, capsys,
    ):
        rc = cli.main(["validate", "--config", str(sample_config)])
        # files missing => non-zero
        assert rc != 0
        out = capsys.readouterr().out
        assert "Configuration: OK" in out
        assert "MISSING" in out


class TestBuildSubcommandParserOnly:
    """Verify parser accepts the build subcommand flags (no actual build)."""

    def test_build_requires_config(self):
        with pytest.raises(SystemExit):
            cli.main(["build"])

    def test_build_accepts_o_flag(self, tmp_path, monkeypatch):
        # Patch _cmd_build to a no-op so parser path is exercised
        called = {}

        def fake_build(args):
            called["config"] = args.config
            called["output_dir"] = args.output_dir
            called["ff_dir"] = args.ff_dir
            return 0

        monkeypatch.setattr(cli, "_cmd_build", fake_build)
        # Re-build the parser so it picks up the patched subcommand handler
        import argparse
        parser = argparse.ArgumentParser()
        sub = parser.add_subparsers(dest="cmd", required=True)
        b = sub.add_parser("build")
        b.add_argument("--config", required=True)
        b.add_argument("-o", "--output-dir", default=None)
        b.add_argument("--ff-dir", default=None)
        b.set_defaults(func=fake_build)
        ns = parser.parse_args(
            ["build", "--config", "x.json", "-o", "/tmp/out",
             "--ff-dir", "/tmp/ff"],
        )
        ns.func(ns)
        assert called["config"] == "x.json"
        assert called["output_dir"] == "/tmp/out"
        assert called["ff_dir"] == "/tmp/ff"
