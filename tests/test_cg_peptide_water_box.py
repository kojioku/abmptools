# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.water_box module."""
from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from abmptools.cg.peptide.water_box import (
    MARTINI_W_DENSITY,
    make_martini_water_box,
)


class TestMakeMartiniWaterBox:

    @patch("abmptools.cg.peptide.water_box.run_command")
    def test_calls_gmx_insert_molecules(self, mock_run, tmp_path):
        out = tmp_path / "water.gro"

        def fake_run(cmd, **kw):
            # Simulate gmx insert-molecules writing the box
            out.write_text("water box\n   1000\n...\n")
            return MagicMock(returncode=0)

        mock_run.side_effect = fake_run

        result = make_martini_water_box(out, box_nm=5.0)
        assert result == out
        assert out.exists()

        cmd = mock_run.call_args.args[0]
        assert cmd[0] == "gmx"
        assert "insert-molecules" in cmd
        assert "-nmol" in cmd
        # 5 nm cubic box passed
        bi = cmd.index("-box")
        assert cmd[bi + 1] == "5.000"
        assert cmd[bi + 2] == "5.000"
        assert cmd[bi + 3] == "5.000"
        # nmol target near MARTINI_W_DENSITY * 5^3 = 1045
        ni = cmd.index("-nmol")
        nmol = int(cmd[ni + 1])
        assert 900 < nmol < 1100

    @patch("abmptools.cg.peptide.water_box.run_command")
    def test_seed_gro_cleaned_up_on_success(self, mock_run, tmp_path):
        out = tmp_path / "water.gro"

        def fake_run(cmd, **kw):
            out.write_text("X")
            return MagicMock(returncode=0)

        mock_run.side_effect = fake_run
        make_martini_water_box(out, box_nm=4.0)

        leftovers = list(tmp_path.glob("_seed_*"))
        assert leftovers == []

    @patch("abmptools.cg.peptide.water_box.run_command")
    def test_seed_gro_cleaned_up_on_failure(self, mock_run, tmp_path):
        from abmptools.cg.peptide._subprocess import CommandError

        out = tmp_path / "water.gro"
        mock_run.side_effect = CommandError(
            "gmx insert-molecules ...", 1, "synthetic error"
        )

        with pytest.raises(CommandError):
            make_martini_water_box(out)

        leftovers = list(tmp_path.glob("_seed_*"))
        assert leftovers == []

    @patch("abmptools.cg.peptide.water_box.run_command")
    def test_failure_no_output_raises_runtime_error(self, mock_run, tmp_path):
        out = tmp_path / "water.gro"
        # gmx returns 0 but doesn't actually produce the output file
        mock_run.return_value = MagicMock(returncode=0)
        with pytest.raises(RuntimeError, match="failed to produce"):
            make_martini_water_box(out)

    @patch("abmptools.cg.peptide.water_box.run_command")
    def test_custom_density(self, mock_run, tmp_path):
        out = tmp_path / "water.gro"

        def fake_run(cmd, **kw):
            out.write_text("X")
            return MagicMock(returncode=0)

        mock_run.side_effect = fake_run
        make_martini_water_box(out, box_nm=5.0, n_w_per_nm3=10.0)
        cmd = mock_run.call_args.args[0]
        ni = cmd.index("-nmol")
        # 5^3 * 10 = 1250
        assert int(cmd[ni + 1]) == 1250

    @patch("abmptools.cg.peptide.water_box.run_command")
    def test_custom_gmx_path(self, mock_run, tmp_path):
        out = tmp_path / "water.gro"

        def fake_run(cmd, **kw):
            out.write_text("X")
            return MagicMock(returncode=0)

        mock_run.side_effect = fake_run
        make_martini_water_box(out, gmx_path="/opt/gmx-2024/bin/gmx")
        cmd = mock_run.call_args.args[0]
        assert cmd[0] == "/opt/gmx-2024/bin/gmx"

    def test_invalid_box_nm_raises(self, tmp_path):
        with pytest.raises(ValueError, match="box_nm"):
            make_martini_water_box(tmp_path / "x.gro", box_nm=0.0)

    def test_invalid_density_raises(self, tmp_path):
        with pytest.raises(ValueError, match="n_w_per_nm3"):
            make_martini_water_box(tmp_path / "x.gro", n_w_per_nm3=-1.0)

    def test_density_constant_is_realistic(self):
        # 1 W = 4 H2O = 72 g/mol; rho=1 g/cm^3 => ~8.36 W/nm^3
        assert 8.0 < MARTINI_W_DENSITY < 9.0
