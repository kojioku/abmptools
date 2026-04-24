# -*- coding: utf-8 -*-
"""
Tests for the HF-only fallback in ``abmptools.abinit_io.read_ifie``.

Before the fallback, MP2-IFIE section detection was hard-wired, so
method='HF' ABINIT-MP outputs produced an empty IFIE table (→ χ = 0
downstream). The fallback detects HF-only outputs, starts reading from
``## HF-IFIE`` instead, and zeros out cols 5-6 so the MP2-shaped
``[..., HF/Hartree, MP2/Hartree, ...]`` contract still holds
(MP2 contribution = 0 in HF-only mode).
"""
from __future__ import annotations

from abmptools.abinit_io import abinit_io


MP2_SNIPPET = """\
     ## MP2-IFIE

                 IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE   PR-TYPE1
                            / A      APPROX.   / Hartree  / Hartree  / Hartree
        ----------------------------------------------------------------------
               2       1    3.456789   F        -0.012345   -0.006789   -0.006789
               3       1    4.123456   F        -0.005678   -0.003456   -0.003456
     ## Mulliken charges
"""


HF_SNIPPET = """\
     ## HF-IFIE

                      IJ-PAIR    DIST     DIMER-ES   HF-IFIE        HF-IFIE        HF-IFIE
                                 / A      APPROX.   / Hartree      / kcal/mol     / kJ/mol
        ----------------------------------------------------------------------------------
               2       1    3.456789   F        -0.012345       -7.746450      -32.411490
               3       1    4.123456   F        -0.005678       -3.562903      -14.903495
     ## Mulliken charges
"""


def _write(tmp_path, snippet):
    path = tmp_path / "out.out"
    path.write_text(snippet)
    return str(path)


def test_mp2_output_returns_full_ifie_table(tmp_path):
    """MP2 出力: 従来どおり MP2-IFIE 表をそのまま取る (後方互換)。"""
    io = abinit_io()
    path = _write(tmp_path, MP2_SNIPPET)
    table = io.read_ifie(path)
    assert len(table) == 2
    # cols: [frag_i, frag_j, DIST, F/T, HF/Hartree, MP2/Hartree, PR-TYPE1]
    assert table[0][0] == "2"
    assert table[0][1] == "1"
    # HF/Hartree (col 4) preserved
    assert float(table[0][4]) == -0.012345
    # MP2/Hartree (col 5) non-zero in MP2 output
    assert float(table[0][5]) == -0.006789


def test_hf_output_triggers_fallback_with_zero_mp2(tmp_path):
    """HF-only 出力: HF-IFIE 表を読んで col 5-6 は 0 に上書き。"""
    io = abinit_io()
    path = _write(tmp_path, HF_SNIPPET)
    table = io.read_ifie(path)
    assert len(table) == 2
    # HF/Hartree in col 4 preserved
    assert float(table[0][4]) == -0.012345
    # MP2/Hartree (col 5) must be 0 so downstream HF+MP2 sum == HF
    assert float(table[0][5]) == 0.0
    assert float(table[0][6]) == 0.0


def test_hf_output_downstream_sum_gives_pure_hf():
    """HF fallback 経由で (col 4) + (col 5) = HF/Hartree のみが残ることを確認。

    downstream code (`agg_mbresult.getcontactifie`) は col 4 + col 5 を
    Hartree で累積してから 627.5095× で kcal/mol に変換する。HF-only では
    col 5 = 0 なので、結果は純 HF の kcal/mol sum になるはず。
    """
    io = abinit_io()
    import tempfile
    import os
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "out.out")
        with open(path, "w") as f:
            f.write(HF_SNIPPET)
        table = io.read_ifie(path)

    # HF only: col 4 の和が 2 行分
    hf_sum_hartree = sum(float(row[4]) for row in table)
    mp2_sum_hartree = sum(float(row[5]) for row in table)
    assert abs(hf_sum_hartree - (-0.012345 - 0.005678)) < 1e-10
    assert mp2_sum_hartree == 0.0
    # 下流の変換 (Hartree → kcal/mol)
    hf_kcal = (hf_sum_hartree + mp2_sum_hartree) * 627.5095
    assert abs(hf_kcal - (-11.31)) < 0.01   # ~ -11.31 kcal/mol
