# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.mmgbsa.cli."""
from __future__ import annotations

import io
import json
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from abmptools.genesis.mmgbsa.cli import main as cli_main
from abmptools.genesis.mmgbsa.models import MMGBSABuildConfig


FIXTURE_LOG_DIR = (
    Path(__file__).parent / "fixtures" / "gbsa_logs" / "3772L-rename"
)


# ---------------------------------------------------------------------------
# example
# ---------------------------------------------------------------------------

class TestExample:
    def test_outputs_valid_json(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = cli_main(["example"])
        assert rc == 0
        data = json.loads(buf.getvalue())
        # Default example has 4 POC targets at ligand_resno=201.
        assert len(data["targets"]) == 4
        assert data["targets"][0]["ligand_resno"] == 201

    def test_roundtrip_through_from_json(self, tmp_path):
        buf = io.StringIO()
        with redirect_stdout(buf):
            cli_main(["example"])
        cfg_path = tmp_path / "ex.json"
        cfg_path.write_text(buf.getvalue())
        cfg = MMGBSABuildConfig.from_json(str(cfg_path))
        assert cfg.project_name == "poc_reproduction"
        assert len(cfg.targets) == 4


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

class TestValidate:
    def test_validate_default_example(self, tmp_path):
        gen_buf = io.StringIO()
        with redirect_stdout(gen_buf):
            cli_main(["example"])
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(gen_buf.getvalue())

        out_buf = io.StringIO()
        with redirect_stdout(out_buf):
            rc = cli_main(["validate", "--config", str(cfg_path)])
        out = out_buf.getvalue()
        assert "Configuration: OK" in out
        assert "External tools:" in out
        assert "Targets preview:" in out
        # CI environments might or might not have all tools; accept 0 or 1.
        assert rc in (0, 1)


# ---------------------------------------------------------------------------
# analyze (replay against fixture logs)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not FIXTURE_LOG_DIR.is_dir(),
    reason="POC log fixtures missing",
)
class TestAnalyzeReplay:
    def test_analyze_via_target_dirs(self, tmp_path):
        out_dir = tmp_path / "ana"
        rc = cli_main([
            "analyze",
            "--target-dirs", str(FIXTURE_LOG_DIR),
            "--out-dir", str(out_dir),
        ])
        assert rc == 0
        assert (out_dir / "analysis_results.csv").exists()
        assert (out_dir / "dg_bind_plot.png").exists()


# ---------------------------------------------------------------------------
# Parser shape
# ---------------------------------------------------------------------------

class TestParser:
    def test_no_subcommand_errors(self):
        with pytest.raises(SystemExit):
            cli_main([])

    def test_pipeline_requires_config_or_input(self, tmp_path):
        err_buf = io.StringIO()
        with redirect_stderr(err_buf):
            rc = cli_main(["pipeline"])
        assert rc == 1
        assert "--config" in err_buf.getvalue()

    def test_divide_subcommand_accepts_config(self, tmp_path):
        # Just make sure argparse accepts it (handler runs ahead of
        # subprocess — without atdyn / Biopython it will fail, but
        # parser shouldn't reject).
        gen_buf = io.StringIO()
        with redirect_stdout(gen_buf):
            cli_main(["example"])
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(gen_buf.getvalue())

        # divide will fail because the input/ folder doesn't exist
        # in the example; we only verify parsing works.
        with pytest.raises((FileNotFoundError, SystemExit)):
            cli_main(["divide", "--config", str(cfg_path)])
