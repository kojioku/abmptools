"""
Regression test: IMC BDF at rec=900 should give a stable classification.

Baseline (Luzar-Chandler defaults, T=450 K IMC amorphous, per-COOH state):
  dual=10, chain=41, single=38, free=36  (sum=125, one mol = one COOH)
  amide accept=49, amide free=76
  n_hbonds_cc=31, n_hbonds_ca=50

The 4-species split mirrors Yuan et al. (2015) Mol. Pharm. 12, 4518 NMR
deconvolution:
  cyclic dimer (~179 ppm) → dual
  carboxylic acid chain (~176 ppm, disordered chains) → chain
  COOH-amide (~172 ppm) → single
  free COOH (~170 ppm) → free
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
        do_copy_uncolored=False,
        do_plot=False,
        verbose=False,
    )
    a = Analyzer(cfg)
    a.load()
    results = a.run()
    assert len(results) == 1
    fr = results[0]
    cls = fr.classification
    # Mol-level fields (representative role used by colorizer)
    assert abs(fr.n_dual_mols - 10) <= 2, f"dual={fr.n_dual_mols}"
    assert abs(fr.n_chain_mols - 41) <= 5, f"chain={fr.n_chain_mols}"
    assert abs(fr.n_single_mols - 38) <= 5, f"single={fr.n_single_mols}"
    assert abs(fr.n_free_mols - 36) <= 5, f"free={fr.n_free_mols}"
    assert abs(fr.n_hbonds_cc - 31) <= 3, f"cc={fr.n_hbonds_cc}"
    assert abs(fr.n_hbonds_ca - 50) <= 5, f"ca={fr.n_hbonds_ca}"
    assert (fr.n_dual_mols + fr.n_chain_mols
            + fr.n_single_mols + fr.n_free_mols) == 125

    # Per-functional-group fields (new in v1.27 candidate)
    assert cls.n_carboxyls == 125
    assert cls.n_amides == 125
    # IMC: 1 mol = 1 COOH + 1 amide, so per-COOH counts == per-mol counts
    assert cls.n_carboxyls_dual == fr.n_dual_mols
    assert cls.n_carboxyls_chain == fr.n_chain_mols
    assert cls.n_carboxyls_single == fr.n_single_mols
    assert cls.n_carboxyls_free == fr.n_free_mols
    assert (cls.n_carboxyls_dual + cls.n_carboxyls_chain
            + cls.n_carboxyls_single + cls.n_carboxyls_free) == 125
    # Amide accept count ≤ ca H-bond count (one amide can accept multiple)
    assert cls.n_amides_accept <= fr.n_hbonds_ca
    assert cls.n_amides_accept + cls.n_amides_free == 125
    # Ratios in [0, 1] and sum to 1.0
    for ratio in (cls.ratio_carboxyl_dual, cls.ratio_carboxyl_chain,
                  cls.ratio_carboxyl_single, cls.ratio_carboxyl_free,
                  cls.ratio_amide_accept):
        assert 0.0 <= ratio <= 1.0
    assert abs(cls.ratio_carboxyl_dual + cls.ratio_carboxyl_chain
               + cls.ratio_carboxyl_single + cls.ratio_carboxyl_free
               - 1.0) < 1e-9


def test_imc_strict_yields_fewer_hbonds():
    from abmptools.hbond import Analyzer, AnalyzerConfig
    cfg_lc = AnalyzerConfig(
        bdf_path=IMC_BDF, out_prefix="/tmp/_test_imc_lc",
        criteria_mode="luzar-chandler",
        do_colorize=False, do_copy_uncolored=False,
        do_plot=False, verbose=False,
    )
    cfg_st = AnalyzerConfig(
        bdf_path=IMC_BDF, out_prefix="/tmp/_test_imc_st",
        criteria_mode="strict",
        do_colorize=False, do_copy_uncolored=False,
        do_plot=False, verbose=False,
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
        do_colorize=False, do_copy_uncolored=False,
        do_plot=False, verbose=False,
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


def test_uncolored_copy_keeps_mol_name(tmp_path):
    """`<prefix>.bdf` is a plain copy with Mol_Name preserved (J-OCTA pre-render)."""
    from abmptools.hbond import Analyzer, AnalyzerConfig
    cfg = AnalyzerConfig(
        bdf_path=IMC_BDF,
        out_prefix=str(tmp_path / "imc"),
        do_colorize=False, do_copy_uncolored=True,
        do_plot=False, verbose=False,
    )
    a = Analyzer(cfg)
    a.load()
    a.run()
    paths = a.write_outputs()
    uncolored = paths.get("uncolored")
    assert uncolored is not None
    assert uncolored.endswith(".bdf")

    from UDFManager import UDFManager
    u = UDFManager(uncolored)
    u.jump(0)
    n_mol = u.size("Set_of_Molecules.molecule[]")
    assert n_mol == 125
    # Mol_Name on every molecule MUST NOT have been renamed
    for i in range(n_mol):
        name = u.get(f"Set_of_Molecules.molecule[{i}].Mol_Name")
        assert "_DUAL" not in name
        assert "_SINGLE" not in name
        assert "_FREE" not in name


def test_classification_csv_per_group(tmp_path):
    """`<prefix>_classification.csv` lists per-functional-group roles."""
    from abmptools.hbond import Analyzer, AnalyzerConfig
    cfg = AnalyzerConfig(
        bdf_path=IMC_BDF,
        out_prefix=str(tmp_path / "imc"),
        do_colorize=False, do_copy_uncolored=False,
        do_plot=False, verbose=False,
    )
    a = Analyzer(cfg)
    a.load()
    a.run()
    paths = a.write_outputs()
    assert "classification" in paths
    import csv
    with open(paths["classification"]) as f:
        rows = list(csv.DictReader(f))
    # 125 carboxyl + 125 amide rows per record
    carb_rows = [r for r in rows if r["group_type"] == "carboxyl"]
    amide_rows = [r for r in rows if r["group_type"] == "amide"]
    assert len(carb_rows) == 125
    assert len(amide_rows) == 125
    # roles drawn from expected sets (v1.27 4-species)
    assert all(
        r["role"] in {"dual", "chain", "single", "free"} for r in carb_rows
    )
    assert all(r["role"] in {"accept", "free"} for r in amide_rows)
