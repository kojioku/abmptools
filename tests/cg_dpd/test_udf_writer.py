# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_udf_writer.py — Route R1 (Cognac DPD 入力 UDF plain text writer)."""
from __future__ import annotations

import re

import pytest


def test_udf_writer_cholesterol_e2e(tmp_path, cholesterol_cg, sample_aij_a_mode):
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"],
        aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    udf = builder.build_udf(tmp_path / "chol_uin.udf")
    assert udf.exists()
    text = udf.read_text()

    # \include{"cognac112.udf"} で class 定義を J-OCTA 環境に委譲 (権利配慮)
    assert 'cognac112.udf' in text
    # EngineType: COGNAC
    assert 'EngineType:"COGNAC"' in text

    # 4 top-level sections 完備
    for section in ("Simulation_Conditions:{", "Initial_Structure:{",
                    "Molecular_Attributes:{", "Set_of_Molecules:{"):
        assert section in text, f"missing section: {section}"

    # Atom_Type: 全 5 segment
    for i in range(5):
        assert f'"P{i}",1.00000000000000' in text or f'"P{i}",1.0' in text, \
            f"Atom_Type P{i} missing"

    # Bond_Potential: cholesterol bond12 = 4 entry (multi-line regex)
    bonds = re.findall(r'"P\d-P\d",\s+"Harmonic"', text)
    assert len(bonds) == 4, f"expected 4 Bond_Potential, got {len(bonds)}"

    # Angle_Potential: cholesterol angle13 = 3 entry
    angles = re.findall(r'"P\d-P\d-P\d",\s+"Theta"', text)
    assert len(angles) == 3, f"expected 3 Angle_Potential, got {len(angles)}"

    # Pair_Interaction: DPD type で aij.dat の全 15 pair
    pair_dpd = re.findall(r'"P\d-P\d",\s+"DPD",\s+"P\d",\s+"P\d"', text)
    assert len(pair_dpd) == 15, f"expected 15 Pair_Interaction DPD, got {len(pair_dpd)}"

    # Set_of_Molecules: 'chol' molecule entry
    assert '"chol"' in text

    # syntactic validity: brace と bracket の balance
    assert text.count("{") == text.count("}"), "brace imbalance"
    assert text.count("[") == text.count("]"), "bracket imbalance"


def test_udf_writer_chi_mode_auto_converts(tmp_path, cholesterol_cg, sample_aij_chi_mode):
    """chi モードの aij.dat でも UDF 内では a 値で出力される。"""
    from abmptools.cg.dpd import CGDpdBuilder, read_monomer, assign_particle_names
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"],
        aij=sample_aij_chi_mode,  # chi mode
        calc_sett=cholesterol_cg["calc_sett"],
    )
    udf = builder.build_udf(tmp_path / "chi_uin.udf")
    text = udf.read_text()
    # chi=-4.108, aii=25.0 → a≈11.58 が UDF に書き込まれているはず
    assert "11.5" in text, "chi → a 変換結果が UDF に反映されていない"


def test_udf_writer_custom_include(tmp_path, cholesterol_cg, sample_aij_a_mode):
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    udf = builder.build_udf(
        tmp_path / "custom.udf", include_file="cognac90.udf",
    )
    text = udf.read_text()
    assert 'cognac90.udf' in text
    assert 'cognac112.udf' not in text
