"""
Regression test: IMC BDF at rec=900 should give a stable classification.

Baseline (Luzar-Chandler defaults, T=450 K IMC amorphous):
  dual=10, single=73, free=42, hbonds_cc=31, hbonds_ca=50
"""
import os

import pytest

IMC_BDF = "/home/okuwaki/llm-project/SI/IMC_result450.0_out_rec900.bdf"

pytestmark = pytest.mark.skipif(
    not os.path.exists(IMC_BDF), reason="IMC BDF not available"
)


def test_imc_counts_luzar_chandler():
    from abmptools.hbond import Analyzer, AnalyzerConfig
    cfg = AnalyzerConfig(
        bdf_path=IMC_BDF,
        out_prefix="/tmp/_test_imc_hbond",
        criteria_mode="luzar-chandler",
        do_colorize=False,
        do_plot=False,
        verbose=False,
    )
    a = Analyzer(cfg)
    a.load()
    results = a.run()
    assert len(results) == 1
    fr = results[0]
    # ±2 tolerance to allow for benign numeric variations
    assert abs(fr.n_dual_mols - 10) <= 2, f"dual={fr.n_dual_mols}"
    assert abs(fr.n_single_mols - 73) <= 5, f"single={fr.n_single_mols}"
    assert abs(fr.n_free_mols - 42) <= 5, f"free={fr.n_free_mols}"
    assert abs(fr.n_hbonds_cc - 31) <= 3, f"cc={fr.n_hbonds_cc}"
    assert abs(fr.n_hbonds_ca - 50) <= 5, f"ca={fr.n_hbonds_ca}"
    # invariant: dual + single + free == 125
    assert fr.n_dual_mols + fr.n_single_mols + fr.n_free_mols == 125


def test_imc_strict_yields_fewer_hbonds():
    from abmptools.hbond import Analyzer, AnalyzerConfig
    cfg_lc = AnalyzerConfig(
        bdf_path=IMC_BDF, out_prefix="/tmp/_test_imc_lc",
        criteria_mode="luzar-chandler",
        do_colorize=False, do_plot=False, verbose=False
    )
    cfg_st = AnalyzerConfig(
        bdf_path=IMC_BDF, out_prefix="/tmp/_test_imc_st",
        criteria_mode="strict",
        do_colorize=False, do_plot=False, verbose=False
    )
    n_cc_lc = Analyzer(cfg_lc).load() or Analyzer(cfg_lc).run()  # noqa: E501

    a_lc = Analyzer(cfg_lc)
    a_lc.load()
    r_lc = a_lc.run()[0]

    a_st = Analyzer(cfg_st)
    a_st.load()
    r_st = a_st.run()[0]

    assert r_st.n_hbonds_cc <= r_lc.n_hbonds_cc
    assert r_st.n_hbonds_ca <= r_lc.n_hbonds_ca


def test_colorize_round_trip(tmp_path):
    """colorize_udf renames mols and adds Draw_Attributes entries."""
    from abmptools.hbond import (
        Analyzer, AnalyzerConfig, colorize_udf
    )
    out_path = str(tmp_path / "imc_colored.bdf")
    cfg = AnalyzerConfig(
        bdf_path=IMC_BDF,
        out_prefix=str(tmp_path / "imc"),
        do_colorize=False, do_plot=False, verbose=False
    )
    a = Analyzer(cfg)
    a.load()
    res = a.run()
    cls = res[0].classification
    colorize_udf(IMC_BDF, out_path, cls, base_mol_name="IMC")

    from UDFManager import UDFManager
    u = UDFManager(out_path)
    u.jump(0)
    n_da = u.size("Draw_Attributes.Molecule[]")
    assert n_da >= 3
    names = set()
    for i in range(n_da):
        names.add(u.get(f"Draw_Attributes.Molecule[{i}].Mol_Name"))
    assert {"IMC_DUAL", "IMC_SINGLE", "IMC_FREE"}.issubset(names)
