"""Unit tests for abmptools.hbond.distance_dist."""
from __future__ import annotations

import csv
import os
import tempfile

import numpy as np
import pytest

from abmptools.hbond.analyzer import FrameResult
from abmptools.hbond.classifier import (
    CarboxylRole, FunctionalGroupClassification,
)
from abmptools.hbond.distance_dist import (
    aggregate_distance_angle, aggregate_distances_generic,
    aggregate_distances_imc, compute_distance_stats, default_bin_edges,
    write_distance_histogram_csv, write_distance_stats_csv,
)
from abmptools.hbond.hbond_detector import HBond


def _hb(donor_mol, acceptor_mol, d_da, angle=160.0):
    return HBond(
        donor_mol=donor_mol, donor_d=3, donor_h=40,
        acceptor_mol=acceptor_mol, acceptor_a=4,
        d_da=d_da, d_ha=d_da - 1.0, angle=angle,
    )


def test_default_bin_edges_count():
    edges = default_bin_edges(d_min=2.0, d_max=3.6, bin_width=0.05)
    assert edges[0] == pytest.approx(2.0)
    assert edges[-1] == pytest.approx(3.6)
    assert edges.size == 33  # 32 bins between 2.0 and 3.6


def test_compute_distance_stats_basic():
    edges = np.linspace(2.0, 3.0, 11)   # 0.1 Å bins
    s = compute_distance_stats([2.5, 2.6, 2.5, 2.7], "x", edges)
    assert s is not None
    assert s.n == 4
    assert s.mean == pytest.approx(2.575)
    assert s.median == pytest.approx(2.55)
    assert s.peak == pytest.approx(2.55, abs=0.05)
    assert s.std == pytest.approx(np.std([2.5, 2.6, 2.5, 2.7]))
    assert s.counts.sum() == 4


def test_compute_distance_stats_empty_returns_none():
    edges = default_bin_edges()
    assert compute_distance_stats([], "empty", edges) is None


def _imc_frame_with_dual_and_single():
    """Mol 0 ↔ 1 dual COOH-COOH, mol 0 → mol 2 single COOH (no amide)."""
    role_0 = CarboxylRole(
        mol_index=0, carboxyl_index=0, role="dual",
        dual_partners={(1, 0)},
    )
    role_1 = CarboxylRole(
        mol_index=1, carboxyl_index=0, role="dual",
        dual_partners={(0, 0)},
    )
    role_2 = CarboxylRole(
        mol_index=2, carboxyl_index=0, role="free",
    )
    cls = FunctionalGroupClassification(
        n_molecules=3, roles=[], dual_pairs=[(0, 1)], single_pairs=[],
        n_dual_mols=2, n_single_mols=0, n_free_mols=1,
        n_carboxyls=3, n_carboxyls_dual=2, n_carboxyls_free=1,
        carboxyl_roles=[role_0, role_1, role_2], amide_roles=[],
    )
    hb_01 = _hb(0, 1, 2.70, angle=170.0)
    hb_10 = _hb(1, 0, 2.80, angle=165.0)
    hb_02 = _hb(0, 2, 3.10, angle=140.0)
    return FrameResult(
        record=0, n_dual_mols=2, n_single_mols=0, n_free_mols=1,
        n_hbonds_cc=3, n_hbonds_ca=0, classification=cls,
        hbonds_cc=[hb_01, hb_10, hb_02], hbonds_ca=[],
    )


def test_aggregate_imc_dual_subset():
    fr = _imc_frame_with_dual_and_single()
    edges = default_bin_edges()
    out = aggregate_distances_imc([fr], edges)

    assert out["all"].n == 3
    assert out["COOH-COOH"].n == 3
    assert out["COOH-COOH (dual)"].n == 2
    assert out["COOH-COOH (chain/single)"].n == 1
    # No amide acceptor in this frame → empty bucket should be filtered out
    assert "COOH-amide" not in out

    dual_distances = np.array([2.70, 2.80])
    assert out["COOH-COOH (dual)"].mean == pytest.approx(dual_distances.mean())


def test_aggregate_imc_with_amide():
    fr = _imc_frame_with_dual_and_single()
    # Add a COOH→amide H-bond at 2.9 Å
    fr.hbonds_ca = [_hb(0, 2, 2.90, angle=155.0)]
    edges = default_bin_edges()
    out = aggregate_distances_imc([fr], edges)
    assert out["COOH-amide"].n == 1
    assert out["COOH-amide"].mean == pytest.approx(2.90)
    assert out["all"].n == 4   # cc(3) + ca(1)


def test_aggregate_generic_per_pair_type():
    hb_aa = _hb(0, 1, 2.65)
    hb_ab = _hb(1, 2, 3.00)
    fr = FrameResult(
        record=0, n_dual_mols=0, n_single_mols=0, n_free_mols=0,
        n_hbonds_cc=0, n_hbonds_ca=0, classification=None,
        hbonds_cc=[], hbonds_ca=[],
        hbonds_by_pair_type={
            ("hydroxyl", "hydroxyl_O"): [hb_aa],
            ("hydroxyl", "amide_O"): [hb_ab],
        },
    )
    edges = default_bin_edges()
    out = aggregate_distances_generic([fr], edges)
    assert out["all"].n == 2
    assert out["hydroxyl->hydroxyl_O"].n == 1
    assert out["hydroxyl->amide_O"].n == 1
    assert out["hydroxyl->hydroxyl_O"].mean == pytest.approx(2.65)
    assert out["hydroxyl->amide_O"].mean == pytest.approx(3.00)


def test_aggregate_distance_angle_imc():
    fr = _imc_frame_with_dual_and_single()
    fr.hbonds_ca = [_hb(0, 2, 2.90, angle=155.0)]
    d, a = aggregate_distance_angle([fr], "imc")
    assert d.size == 4
    assert a.size == 4
    assert set(np.round(d, 2)) == {2.70, 2.80, 3.10, 2.90}
    assert set(np.round(a, 1)) == {170.0, 165.0, 140.0, 155.0}


def test_csv_writers_round_trip():
    fr = _imc_frame_with_dual_and_single()
    edges = default_bin_edges()
    out = aggregate_distances_imc([fr], edges)
    with tempfile.TemporaryDirectory() as td:
        stats_path = os.path.join(td, "stats.csv")
        hist_path = os.path.join(td, "hist.csv")
        write_distance_stats_csv(out, stats_path)
        write_distance_histogram_csv(out, hist_path)

        with open(stats_path) as f:
            rows = list(csv.DictReader(f))
        assert {r["label"] for r in rows} == set(out.keys())
        dual = next(r for r in rows if r["label"] == "COOH-COOH (dual)")
        assert int(dual["n"]) == 2
        assert float(dual["mean_DA"]) == pytest.approx(2.75)

        with open(hist_path) as f:
            hist_rows = list(csv.DictReader(f))
        # Each label contributes (n_bins) rows
        n_bins = edges.size - 1
        assert len(hist_rows) == n_bins * len(out)
        total_count = sum(
            int(r["count"]) for r in hist_rows
            if r["label"] == "COOH-COOH (dual)"
        )
        assert total_count == 2
