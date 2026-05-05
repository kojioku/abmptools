# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.cli."""
from __future__ import annotations

import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout

import pytest

from abmptools.genesis.grest.cli import main as cli_main
from abmptools.genesis.grest.models import GrestBuildConfig


class TestExampleSubcommand:
    def test_example_outputs_valid_json(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = cli_main(["example"])
        assert rc == 0
        text = buf.getvalue()
        # Output is valid JSON.
        data = json.loads(text)
        assert data["project_name"] == "grest_run"
        assert data["rest_selection"]["mode"] == "explicit"
        # Default ladder is the POC's 4 replicas.
        assert len(data["replica_temperatures"]["temperatures"]) == 4

    def test_example_roundtrip_through_GrestBuildConfig(self, tmp_path):
        # Capture example output, write to file, load via from_json -> no errors.
        buf = io.StringIO()
        with redirect_stdout(buf):
            cli_main(["example"])
        out = tmp_path / "example.json"
        out.write_text(buf.getvalue())
        cfg = GrestBuildConfig.from_json(str(out))
        assert cfg.input_pdb == "protein.pdb"
        assert cfg.rest_selection.residues == ["1-138"]


class TestValidateSubcommand:
    def test_validate_explicit_mode(self, tmp_path):
        # Generate the default example and feed it back to validate.
        gen_buf = io.StringIO()
        with redirect_stdout(gen_buf):
            cli_main(["example"])
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(gen_buf.getvalue())

        val_buf = io.StringIO()
        with redirect_stdout(val_buf):
            rc = cli_main(["validate", "--config", str(cfg_path)])
        out = val_buf.getvalue()

        # Configuration line exists.
        assert "Configuration: OK" in out
        # Tools section + REST preview present.
        assert "External tools:" in out
        assert "REST selection preview:" in out
        assert "Replica temperature ladder" in out
        assert "rep01" in out
        # Return code reflects whether all required tools were found.
        # We don't assume tleap/spdyn/atdyn/remd_convert are installed in CI,
        # so accept either 0 or 1.
        assert rc in (0, 1)


class TestParserRejectsUnknown:
    def test_no_subcommand_errors(self):
        # argparse exits with SystemExit when subcommand is required.
        with pytest.raises(SystemExit):
            cli_main([])

    def test_unknown_what_token_rejected(self):
        with pytest.raises(SystemExit):
            cli_main(
                ["analyze",
                 "--config", "/dev/null",
                 "--run-dir", "/tmp",
                 "--what", "unknown_task"]
            )
