# -*- coding: utf-8 -*-
"""Tests for abmptools.amorphous.mdp_protocol."""
import os
from pathlib import Path

import pytest

from abmptools.core.system_model import AnnealProtocol
from abmptools.amorphous.mdp_protocol import (
    _format_mdp,
    generate_em_mdp,
    generate_nvt_high_mdp,
    generate_npt_high_mdp,
    generate_anneal_mdp,
    generate_npt_final_mdp,
    write_all_mdp,
)


@pytest.fixture
def protocol():
    return AnnealProtocol(T_high=600.0, T_low=300.0)


# ---------------------------------------------------------------------------
# Stage generators
# ---------------------------------------------------------------------------

def test_generate_em_mdp_contains_integrator(protocol):
    text = generate_em_mdp(protocol)
    assert "integrator" in text
    assert "steep" in text


def test_generate_nvt_high_mdp_contains_gen_vel(protocol):
    text = generate_nvt_high_mdp(protocol)
    assert "gen-vel" in text
    assert "yes" in text


def test_generate_npt_high_mdp_contains_barostat(protocol):
    text = generate_npt_high_mdp(protocol)
    assert "pcoupl" in text
    assert "Parrinello-Rahman" in text


def test_generate_anneal_mdp_contains_annealing(protocol):
    text = generate_anneal_mdp(protocol)
    assert "annealing" in text
    assert "single" in text


def test_generate_npt_final_mdp_ref_t_uses_T_low(protocol):
    text = generate_npt_final_mdp(protocol)
    # ref-t should be T_low (300.0)
    assert "ref-t" in text
    assert str(protocol.T_low) in text


# ---------------------------------------------------------------------------
# write_all_mdp
# ---------------------------------------------------------------------------

def test_write_all_mdp_creates_5_files(protocol, tmp_path):
    paths = write_all_mdp(protocol, str(tmp_path))
    assert len(paths) == 5
    for p in paths:
        assert Path(p).exists()
        assert Path(p).stat().st_size > 0


# ---------------------------------------------------------------------------
# _format_mdp
# ---------------------------------------------------------------------------

def test_format_mdp_with_title():
    text = _format_mdp({"nsteps": 1000}, title="My Title")
    assert "title = My Title" in text
    assert "nsteps" in text


def test_format_mdp_bool_conversion():
    text = _format_mdp({"gen-vel": True, "continuation": False})
    assert "yes" in text
    assert "no" in text
