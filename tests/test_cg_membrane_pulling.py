# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.pulling."""
from __future__ import annotations

import pytest

from abmptools.cg.membrane import pulling
from abmptools.cg.membrane.models import (
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
)


def _cfg():
    return MembraneCGBuildConfig(
        lipids=[LipidMix()],
        peptide=PeptideMembraneSpec(),
        martini_itp_dir="./ff",
    )


# ---------------------------------------------------------------------------
# Re-exports work
# ---------------------------------------------------------------------------

def test_parse_pullx_xvg_reexported():
    assert callable(pulling.parse_pullx_xvg)


def test_find_pbc_center_atom_reexported():
    assert callable(pulling.find_pbc_center_atom)


def test_estimate_initial_pull_coord_reexported():
    assert callable(pulling.estimate_initial_pull_coord)


def test_extract_window_frames_reexported():
    assert callable(pulling.extract_window_frames)


# ---------------------------------------------------------------------------
# write_pulling_mdp_cg
# ---------------------------------------------------------------------------

class TestWritePullingMdpCG:
    def test_writes_pull_mdp(self, tmp_path):
        out = pulling.write_pulling_mdp_cg(
            config=_cfg(), pull_dir=tmp_path,
            pull_init_nm=3.0,
        )
        assert out.exists()
        assert out.name == "pull.mdp"

    def test_content_uses_direction_periodic(self, tmp_path):
        out = pulling.write_pulling_mdp_cg(
            config=_cfg(), pull_dir=tmp_path,
            pull_init_nm=3.0,
        )
        text = out.read_text()
        assert "direction-periodic" in text

    def test_pbc_atom_g1_inserted(self, tmp_path):
        out = pulling.write_pulling_mdp_cg(
            config=_cfg(), pull_dir=tmp_path,
            pull_init_nm=3.0,
            pbc_atom_g1=1234,
        )
        text = out.read_text()
        assert "pull-group1-pbcatom        = 1234" in text

    def test_custom_filename(self, tmp_path):
        out = pulling.write_pulling_mdp_cg(
            config=_cfg(), pull_dir=tmp_path,
            pull_init_nm=3.0,
            filename="custom_pull.mdp",
        )
        assert out.name == "custom_pull.mdp"


# ---------------------------------------------------------------------------
# Smoke: parse_pullx_xvg with a fake file
# ---------------------------------------------------------------------------

def test_parse_pullx_xvg_basic(tmp_path):
    p = tmp_path / "pullx.xvg"
    p.write_text(
        "@ title \"pull\"\n"
        "# comment\n"
        "0.000  -1.500\n"
        "10.0   -1.000\n"
        "20.0   -0.500\n"
    )
    times, zs = pulling.parse_pullx_xvg(str(p))
    assert times == [0.0, 10.0, 20.0]
    assert zs == [-1.5, -1.0, -0.5]


def test_parse_pullx_xvg_skips_short_lines(tmp_path):
    p = tmp_path / "pullx.xvg"
    p.write_text(
        "@ title\n"
        "0.0\n"           # too short, skipped
        "1.0  2.0\n"
    )
    times, _zs = pulling.parse_pullx_xvg(str(p))
    assert times == [1.0]
