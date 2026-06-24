# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_udf_writer.py — Route R1 (Cognac DPD 入力 UDF writer)。

writer は UDFManager ベース (cognac112 スキーマ準拠) に変更されたため、
テストも **実際に UDFManager でロードできること** + **正しいスキーマ位置に
値が入っていること** を検証する (旧 plain-text regex でなく)。

UDFManager + cognac112.udf が無い環境 (OCTA 非導入 / CI) では skip。
"""
from __future__ import annotations

import pytest

from .conftest import requires_cognac


@requires_cognac
def test_udf_writer_loads_and_schema_correct(tmp_path, cholesterol_cg, sample_aij_a_mode):
    """build_udf 出力が UDFManager でロードでき、Pair_Interaction が
    Interactions 配下 (= cognac112 スキーマ準拠) にあることを検証。"""
    from abmptools.cg.dpd import CGDpdBuilder
    from UDFManager import UDFManager

    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"],
        aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    udf = builder.build_udf(tmp_path / "chol_uin.udf")
    assert udf.exists()

    # 実ロード (旧 writer はここで RuntimeError でコケていた)
    u = UDFManager(str(udf))

    # Pair_Interaction は Interactions 配下 (Molecular_Attributes ではない)
    assert u.size("Interactions.Pair_Interaction[]") == 15      # aij 15 pair
    # Interaction_Site_Type は Molecular_Attributes 配下に 5 segment
    assert u.size("Molecular_Attributes.Interaction_Site_Type[]") == 5
    assert u.size("Molecular_Attributes.Atom_Type[]") == 5
    assert u.size("Molecular_Attributes.Bond_Potential[]") == 4
    assert u.size("Molecular_Attributes.Angle_Potential[]") == 3
    assert u.size("Set_of_Molecules.molecule[]") == 1

    # DPD type + Site 名 + a 値が読める
    assert u.get("Interactions.Pair_Interaction[0].Potential_Type") == "DPD"
    s1 = u.get("Interactions.Pair_Interaction[0].Site1_Name")
    s2 = u.get("Interactions.Pair_Interaction[0].Site2_Name")
    assert s1.startswith("P") and s2.startswith("P")
    a0 = u.get("Interactions.Pair_Interaction[0].DPD.a")
    assert isinstance(a0, float) and a0 > 0.0


@requires_cognac
def test_udf_writer_chi_mode_auto_converts(tmp_path, cholesterol_cg, sample_aij_chi_mode):
    """chi モードの aij.dat でも UDF 内では a = aii + chi/0.306 で出力される。"""
    from abmptools.cg.dpd import CGDpdBuilder
    from UDFManager import UDFManager

    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"],
        aij=sample_aij_chi_mode,   # chi mode (A,B,W / aii=25)
        calc_sett=cholesterol_cg["calc_sett"],
    )
    udf = builder.build_udf(tmp_path / "chi_uin.udf")
    u = UDFManager(str(udf))
    n = u.size("Interactions.Pair_Interaction[]")
    assert n >= 1
    # A-B chi=-4.108 → a = 25 + (-4.108)/0.306 ≈ 11.58 が、いずれかの pair に入る
    a_vals = [u.get(f"Interactions.Pair_Interaction[{k}].DPD.a") for k in range(n)]
    assert any(abs(a - (25.0 + (-4.108) / 0.306)) < 1e-2 for a in a_vals), \
        f"chi→a 変換結果が見当たらない: {a_vals}"


@requires_cognac
def test_udf_writer_custom_include_line(tmp_path, cholesterol_cg, sample_aij_a_mode):
    """include_file 引数が出力 UDF の \\include 行に反映される。

    NOTE: UDFManager がスキーマ解決できる必要があるので、解決可能な
    cognac112.udf を別名で tmp に配置して include 名に使う。
    """
    import shutil
    from abmptools.cg.dpd import CGDpdBuilder

    # 解決可能な cognac112.udf を別名でコピー (環境の OCTA udf から探す)
    import os
    src = None
    for cand in ("/home/okuwaki/OCTA85/ENGINES/udf/cognac112.udf",):
        if os.path.exists(cand):
            src = cand
            break
    if src is None:
        pytest.skip("cognac112.udf source not found for custom-include test")
    alt = tmp_path / "mycognac.udf"
    shutil.copyfile(src, alt)

    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    udf = builder.build_udf(tmp_path / "custom.udf", include_file="mycognac.udf")
    text = udf.read_text()
    assert "mycognac.udf" in text
