# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_orchestrator_cli.py — CGDpdBuilder + CLI subcommand."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


def test_builder_from_files_with_calc_sett(cholesterol_cg, sample_aij_a_mode):
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"],
        aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    assert len(builder.spec.monomers) == 1
    assert builder.spec.monomers[0].name == "chol"
    assert builder.spec.calc_sett is not None
    assert builder.spec.aij.mode == "a"


def test_builder_from_files_without_calc_sett(cholesterol_cg, sample_aij_a_mode):
    """R2 (DPM) のみなら calc_sett 省略可。"""
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
    )
    assert builder.spec.calc_sett is None


def test_builder_assign_particle_names(cholesterol_cg, sample_aij_a_mode):
    from abmptools.cg.dpd import CGDpdBuilder
    custom_names = ["A","B","C","D","T"]
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
        particle_names=custom_names,
    )
    assert builder.spec.monomers[0].particle_names == custom_names


def test_segment_names_unique(cholesterol_cg, sample_aij_a_mode):
    """DpdSystemSpec.segment_names() は登場順で unique."""
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    names = builder.spec.segment_names()
    assert names == ["P0","P1","P2","P3","P4"]


# --- CLI ---------------------------------------------------------------------

def test_cli_build_udf(tmp_path, cholesterol_cg, sample_aij_a_mode):
    """`python -m abmptools.cg.dpd build-udf ...` で UDF 生成。"""
    out = tmp_path / "chol_uin.udf"
    result = subprocess.run(
        [
            sys.executable, "-m", "abmptools.cg.dpd", "build-udf",
            "--monomer", str(cholesterol_cg["monomer"]),
            "--aij", str(sample_aij_a_mode),
            "--calc-sett", str(cholesterol_cg["calc_sett"]),
            "--output", str(out),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert out.exists()
    assert 'cognac112.udf' in out.read_text()


def test_cli_build_dpm(
    tmp_path, cholesterol_cg, sample_aij_a_mode,
    dpm_template_path, virtual_mom_template_path,
):
    """`python -m abmptools.cg.dpd build-dpm ...` で dpm 生成。"""
    out_dir = tmp_path / "chol_proj"
    result = subprocess.run(
        [
            sys.executable, "-m", "abmptools.cg.dpd", "build-dpm",
            "--monomer", str(cholesterol_cg["monomer"]),
            "--aij", str(sample_aij_a_mode),
            "--calc-sett", str(cholesterol_cg["calc_sett"]),
            "--template", str(dpm_template_path),
            "--virtual-mom", str(virtual_mom_template_path),
            "--output-dir", str(out_dir),
            "--dpm-filename", "chol.dpm",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert (out_dir / "chol.dpm").exists()
    # monomer-lib + Virtual.mom 配置確認
    assert (out_dir / "chol" / "monomer-lib" / "P0" / "Virtual.mom").exists()
    assert (out_dir / "chol" / "#Message.txt").exists()
