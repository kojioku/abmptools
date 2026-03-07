# -*- coding: utf-8 -*-
"""Tests for abmptools.amorphous.density module."""
import math

import pytest

from abmptools.amorphous.density import (
    AVOGADRO,
    estimate_box_size_nm,
    weight_fractions_to_counts,
)


# ---------------------------------------------------------------------------
# weight_fractions_to_counts
# ---------------------------------------------------------------------------

class TestWeightFractionsToCounts:
    """Tests for weight_fractions_to_counts."""

    def test_single_component_full_fraction(self):
        """Single component with 100% weight fraction returns all molecules."""
        counts = weight_fractions_to_counts([18.0], [1.0], 50)
        assert counts == [50]

    def test_two_equal_components_same_mw(self):
        """Two components with equal MW and equal weight fractions get equal counts."""
        counts = weight_fractions_to_counts([100.0, 100.0], [0.5, 0.5], 100)
        assert counts == [50, 50]

    def test_two_components_different_mw(self):
        """Heavier molecules get fewer counts for the same weight fraction."""
        # wf=0.5 each, MW 100 vs 200 -> raw proportions: 0.5/100=0.005, 0.5/200=0.0025
        # fracs: 0.005/0.0075*60=40, 0.0025/0.0075*60=20
        counts = weight_fractions_to_counts([100.0, 200.0], [0.5, 0.5], 60)
        assert counts == [40, 20]

    def test_counts_sum_to_total(self):
        """Returned counts must always sum to total_molecules."""
        mws = [18.015, 46.07, 180.16]
        wfs = [0.3, 0.3, 0.4]
        total = 200
        counts = weight_fractions_to_counts(mws, wfs, total)
        assert sum(counts) == total

    def test_raises_on_length_mismatch(self):
        """Raises ValueError when list lengths differ."""
        with pytest.raises(ValueError, match="same length"):
            weight_fractions_to_counts([18.0, 46.0], [1.0], 10)

    def test_raises_on_non_positive_mw(self):
        """Raises ValueError if any molecular weight is zero or negative."""
        with pytest.raises(ValueError, match="positive"):
            weight_fractions_to_counts([0.0], [1.0], 10)
        with pytest.raises(ValueError, match="positive"):
            weight_fractions_to_counts([-5.0], [1.0], 10)

    def test_raises_on_bad_weight_fraction_sum(self):
        """Raises ValueError if weight fractions deviate from 1.0 by more than 0.05."""
        with pytest.raises(ValueError, match="sum to"):
            weight_fractions_to_counts([18.0, 46.0], [0.3, 0.3], 10)


# ---------------------------------------------------------------------------
# estimate_box_size_nm
# ---------------------------------------------------------------------------

class TestEstimateBoxSizeNm:
    """Tests for estimate_box_size_nm."""

    def test_water_known_values(self):
        """100 water molecules at 1.0 g/cm^3 gives a predictable box size."""
        mw = 18.015
        n = 100
        density = 1.0
        expected_mass_g = mw * n / AVOGADRO
        expected_vol_cm3 = expected_mass_g / density
        expected_vol_nm3 = expected_vol_cm3 * 1e21
        expected_edge = expected_vol_nm3 ** (1.0 / 3.0)

        result = estimate_box_size_nm([mw], [n], density)
        assert result == pytest.approx(expected_edge, rel=1e-6)

    def test_higher_density_gives_smaller_box(self):
        """Doubling density should reduce the box edge by factor 2^(1/3)."""
        edge_low = estimate_box_size_nm([100.0], [50], density_g_cm3=0.5)
        edge_high = estimate_box_size_nm([100.0], [50], density_g_cm3=1.0)
        assert edge_high < edge_low
        assert edge_low / edge_high == pytest.approx(2.0 ** (1.0 / 3.0), rel=1e-6)

    def test_avogadro_constant(self):
        """AVOGADRO constant has the expected CODATA 2018 value."""
        assert AVOGADRO == pytest.approx(6.02214076e23, rel=1e-9)
