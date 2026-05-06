# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.mmgbsa.analysis."""
from __future__ import annotations

from pathlib import Path

import pytest

# Force non-interactive backend for matplotlib.
import matplotlib
matplotlib.use("Agg")

from abmptools.genesis.mmgbsa.analysis import (
    CSV_HEADER,
    SYSTEM_NAMES,
    aggregate_target,
    analyze,
    compute_dg_bind,
    compute_dg_components,
    parse_step4_log,
    parse_step4_log_file,
    plot_dg_bind,
    write_results_csv,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Real GENESIS atdyn STEP4 output excerpt (from POC sample-genesis logs).
COMPLEX_STEP4 = """\
[STEP4] Compute Single Point Energy for Molecules

            STEP          ENERGY            BOND           ANGLE        DIHEDRAL        IMPROPER         VDWAALS           ELECT       SOLVATION
 --------------- --------------- --------------- --------------- ---------------
               0      -8855.2124        147.8618        539.1255       3479.7236         23.5583      -1269.9594      -9746.0196      -2029.5026
"""

# 3772L-rename gold values (from POC analysis_results.csv).
GOLD_3772L = {
    "complex_e": -8855.2124, "complex_s": -2029.5026,
    "ligand_e":      8.8654, "ligand_s":     -57.5072,
    "receptor_e": -8825.8078, "receptor_s": -2093.8123,
}

#: Standard MM/GBSA ΔG_bind for 3772L gold values.
#:
#: Per GENESIS doc 05_Energy.rst:564, the ENERGY column is the total
#: potential ``U = U_FF + ΔG_solv`` (SOLVATION is already folded in).
#: So ΔG_bind = E_complex - E_ligand - E_receptor
#:            = -8855.2124 - 8.8654 - (-8825.8078)
#:            = -38.2700 kcal/mol  (negative = favourable binding,
#:                                  single-frame minimisation)
GOLD_3772L_DG = -38.2700


# Path to fixture logs (populated from /home/okuwaki/repos/abmptools-sample/...)
FIXTURE_DIR = Path(__file__).parent / "fixtures" / "gbsa_logs" / "3772L-rename"


# ---------------------------------------------------------------------------
# parse_step4_log
# ---------------------------------------------------------------------------

class TestParseStep4Log:
    def test_basic_parse(self):
        e = parse_step4_log(COMPLEX_STEP4)
        assert e.energy == pytest.approx(-8855.2124)
        assert e.solvation == pytest.approx(-2029.5026)

    def test_missing_step4_raises(self):
        with pytest.raises(ValueError, match="STEP4"):
            parse_step4_log("Some output without the marker.\n")

    def test_truncated_after_marker_raises(self):
        text = "[STEP4] Compute Single Point Energy for Molecules\n"
        with pytest.raises(ValueError, match="data row"):
            parse_step4_log(text)

    def test_too_few_tokens_raises(self):
        # 4 lines below the header but only 5 numbers.
        text = (
            "[STEP4] Compute Single Point Energy for Molecules\n"
            "header1\nheader2\nheader3\n"
            " 0 1.0 2.0 3.0 4.0\n"
        )
        with pytest.raises(ValueError, match=">= 9"):
            parse_step4_log(text)


# ---------------------------------------------------------------------------
# parse_step4_log_file (real fixture)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not FIXTURE_DIR.is_dir(),
    reason=f"Log fixtures missing at {FIXTURE_DIR}",
)
class TestParseStep4LogFileFixture:
    def test_complex_log_3772L(self):
        e = parse_step4_log_file(FIXTURE_DIR / "complex.log")
        assert e.energy == pytest.approx(GOLD_3772L["complex_e"], abs=1e-3)
        assert e.solvation == pytest.approx(GOLD_3772L["complex_s"], abs=1e-3)

    def test_ligand_log_3772L(self):
        e = parse_step4_log_file(FIXTURE_DIR / "ligand.log")
        assert e.energy == pytest.approx(GOLD_3772L["ligand_e"], abs=1e-3)
        assert e.solvation == pytest.approx(GOLD_3772L["ligand_s"], abs=1e-3)

    def test_receptor_log_3772L(self):
        e = parse_step4_log_file(FIXTURE_DIR / "receptor.log")
        assert e.energy == pytest.approx(GOLD_3772L["receptor_e"], abs=1e-3)
        assert e.solvation == pytest.approx(GOLD_3772L["receptor_s"], abs=1e-3)


# ---------------------------------------------------------------------------
# compute_dg_bind
# ---------------------------------------------------------------------------

class TestComputeDgBind:
    def test_3772L_gold(self):
        dg = compute_dg_bind(
            complex_e=GOLD_3772L["complex_e"],
            complex_s=GOLD_3772L["complex_s"],
            ligand_e=GOLD_3772L["ligand_e"],
            ligand_s=GOLD_3772L["ligand_s"],
            receptor_e=GOLD_3772L["receptor_e"],
            receptor_s=GOLD_3772L["receptor_s"],
        )
        assert dg == pytest.approx(GOLD_3772L_DG, abs=1e-3)

    def test_zero_complex_zero_components(self):
        # ΔG = 0 - 0 - 0 = 0
        assert compute_dg_bind(0, 0, 0, 0, 0, 0) == 0.0

    def test_no_solvation_double_counting(self):
        # ENERGY column already includes SOLVATION (per GENESIS doc),
        # so ΔG_bind = E_c - E_l - E_r (no need to add S on top).
        dg = compute_dg_bind(
            complex_e=-100.0, complex_s=-50.0,
            ligand_e=-30.0, ligand_s=-15.0,
            receptor_e=-60.0, receptor_s=-30.0,
        )
        # = -100 - (-30) - (-60) = -100 + 30 + 60 = -10
        assert dg == pytest.approx(-10.0)

    def test_components_sum_equals_dg_bind(self):
        # ΔG_bind = ΔE_MM + ΔG_solv (algebraically identical decomposition).
        comps = compute_dg_components(
            complex_e=-100.0, complex_s=-50.0,
            ligand_e=-30.0, ligand_s=-15.0,
            receptor_e=-60.0, receptor_s=-30.0,
        )
        dg = compute_dg_bind(
            complex_e=-100.0, complex_s=-50.0,
            ligand_e=-30.0, ligand_s=-15.0,
            receptor_e=-60.0, receptor_s=-30.0,
        )
        assert comps["dg_mm"] + comps["dg_solv"] == pytest.approx(dg)
        assert comps["dg_bind"] == pytest.approx(dg)
        # ΔG_solv = -50 - (-15) - (-30) = -50 + 15 + 30 = -5
        assert comps["dg_solv"] == pytest.approx(-5.0)
        # ΔE_MM = (E - S) differences
        # egas_c = -100 - (-50) = -50
        # egas_l = -30 - (-15) = -15
        # egas_r = -60 - (-30) = -30
        # ΔE_MM = -50 - (-15) - (-30) = -5
        assert comps["dg_mm"] == pytest.approx(-5.0)


# ---------------------------------------------------------------------------
# aggregate_target (real fixture)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not FIXTURE_DIR.is_dir(),
    reason=f"Log fixtures missing at {FIXTURE_DIR}",
)
class TestAggregateTarget:
    def test_3772L_gold(self):
        r = aggregate_target(FIXTURE_DIR)
        assert r.name == "3772L-rename"
        assert r.complex_e == pytest.approx(GOLD_3772L["complex_e"], abs=1e-3)
        assert r.dg_bind == pytest.approx(GOLD_3772L_DG, abs=1e-3)

    def test_missing_log_raises(self, tmp_path):
        # tmp_path has no logs.
        with pytest.raises(FileNotFoundError, match="Missing"):
            aggregate_target(tmp_path)

    def test_custom_name(self):
        r = aggregate_target(FIXTURE_DIR, name="custom_label")
        assert r.name == "custom_label"


# ---------------------------------------------------------------------------
# write_results_csv
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not FIXTURE_DIR.is_dir(),
    reason=f"Log fixtures missing at {FIXTURE_DIR}",
)
class TestWriteResultsCsv:
    def test_writes_header_and_one_row(self, tmp_path):
        r = aggregate_target(FIXTURE_DIR)
        out_csv = tmp_path / "results.csv"
        write_results_csv([r], out_csv)
        assert out_csv.exists()
        lines = out_csv.read_text().splitlines()
        assert lines[0] == ",".join(CSV_HEADER)
        assert "3772L-rename" in lines[1]
        # ΔG_bind near gold.
        cells = lines[1].split(",")
        assert float(cells[-1]) == pytest.approx(GOLD_3772L_DG, abs=1e-3)


# ---------------------------------------------------------------------------
# plot_dg_bind
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not FIXTURE_DIR.is_dir(),
    reason=f"Log fixtures missing at {FIXTURE_DIR}",
)
class TestPlotDgBind:
    def test_writes_png(self, tmp_path):
        r = aggregate_target(FIXTURE_DIR)
        out_png = tmp_path / "dg.png"
        plot_dg_bind([r], out_png)
        assert out_png.exists()
        assert out_png.stat().st_size > 0  # non-empty

    def test_empty_results_raises(self, tmp_path):
        with pytest.raises(ValueError, match="empty"):
            plot_dg_bind([], tmp_path / "x.png")


# ---------------------------------------------------------------------------
# analyze (top-level)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not FIXTURE_DIR.is_dir(),
    reason=f"Log fixtures missing at {FIXTURE_DIR}",
)
class TestAnalyze:
    def test_full_flow(self, tmp_path):
        result = analyze(
            target_dirs=[FIXTURE_DIR],
            out_csv=tmp_path / "results.csv",
            out_png=tmp_path / "dg.png",
        )
        assert result.csv_path.exists()
        assert result.png_path.exists()
        assert len(result.targets) == 1
        assert result.targets[0].dg_bind == pytest.approx(GOLD_3772L_DG, abs=1e-3)
