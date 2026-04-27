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
    write_wrap_script,
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


# ---------------------------------------------------------------------------
# write_wrap_script: stages override (hybrid 4-stage protocol)
# ---------------------------------------------------------------------------

def test_write_wrap_script_default_uses_openff_stages(tmp_path):
    path = write_wrap_script(str(tmp_path))
    text = Path(path).read_text()
    assert "STAGES=(02_nvt_highT 03_npt_highT 04_anneal 05_npt_final)" in text
    assert "05_npt_final.gro" in text
    assert "05_npt_final.tpr" in text


def test_write_wrap_script_supports_4stage_hybrid(tmp_path):
    """fcews-manybody の 4-stage protocol で stages を上書きできる。"""
    path = write_wrap_script(
        str(tmp_path),
        ndx=None,
        stages=["01_nvt_eq", "02_npt_high",
                "03_npt_low", "04_nvt_sampling"],
        final_stage="04_nvt_sampling",
    )
    text = Path(path).read_text()
    assert "STAGES=(01_nvt_eq 02_npt_high 03_npt_low 04_nvt_sampling)" in text
    assert "04_nvt_sampling.gro" in text
    assert "04_nvt_sampling_pbc.gro" in text
    # default OpenFF final stage shouldn't leak in
    assert "05_npt_final" not in text


def test_write_wrap_script_final_stage_defaults_to_last(tmp_path):
    """final_stage を省略すると stages の末尾が使われる。"""
    path = write_wrap_script(
        str(tmp_path),
        ndx=None,
        stages=["a", "b", "c"],
    )
    text = Path(path).read_text()
    assert 'if [ -f "c.gro" ]' in text
    assert "c_pbc.gro" in text


def test_write_wrap_script_omits_ndx_when_none(tmp_path):
    """ndx=None なら -n flag が出ない。"""
    path = write_wrap_script(str(tmp_path), ndx=None)
    text = Path(path).read_text()
    assert " -n " not in text


def test_write_wrap_script_includes_ndx_when_given(tmp_path):
    """ndx 指定があれば -n flag が出る。"""
    path = write_wrap_script(str(tmp_path), ndx="system.ndx")
    text = Path(path).read_text()
    assert "-n " in text
    assert "system.ndx" in text
