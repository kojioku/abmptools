# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.replica_temperatures."""
from __future__ import annotations

import math

import pytest

from abmptools.genesis.grest.models import ReplicaTemperatureSpec
from abmptools.genesis.grest.replica_temperatures import (
    format_ladder_line,
    generate_ladder,
    ladder_ratios,
)


class TestGenerateLadderManual:
    def test_passthrough(self):
        spec = ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.0, 318.11, 337.11, 357.10],
        )
        ladder = generate_ladder(spec)
        assert ladder == [300.0, 318.11, 337.11, 357.10]

    def test_two_replicas(self):
        spec = ReplicaTemperatureSpec(
            mode="manual", temperatures=[300.0, 350.0]
        )
        ladder = generate_ladder(spec)
        assert ladder == [300.0, 350.0]

    def test_rounding_to_2dp(self):
        # Inputs at 3+ d.p. should be rounded to 2 d.p.
        spec = ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.123, 318.456, 337.111, 357.099],
        )
        ladder = generate_ladder(spec)
        # 357.099 -> 357.10 after rounding.
        assert ladder == [300.12, 318.46, 337.11, 357.10]


class TestGenerateLadderAuto:
    def test_geometric_4_replicas_matches_poc(self):
        # POC ladder: 300.00 318.11 337.11 357.10 (T_min=300, T_max=357.10).
        # Geometric series approximation should be within 0.5 K.
        spec = ReplicaTemperatureSpec(
            mode="auto",
            method="geometric",
            T_min=300.0,
            T_max=357.10,
            n_replicas=4,
        )
        ladder = generate_ladder(spec)
        assert len(ladder) == 4
        assert ladder[0] == 300.0
        # T_max anchored exactly via final-element forcing.
        assert ladder[-1] == 357.10
        # Intermediate values close to POC.
        assert math.isclose(ladder[1], 318.11, abs_tol=0.5)
        assert math.isclose(ladder[2], 337.11, abs_tol=0.5)

    def test_geometric_2_replicas(self):
        spec = ReplicaTemperatureSpec(
            mode="auto",
            method="geometric",
            T_min=300.0,
            T_max=400.0,
            n_replicas=2,
        )
        ladder = generate_ladder(spec)
        assert ladder == [300.0, 400.0]

    def test_geometric_8_replicas_monotone(self):
        spec = ReplicaTemperatureSpec(
            mode="auto",
            method="geometric",
            T_min=300.0,
            T_max=420.0,
            n_replicas=8,
        )
        ladder = generate_ladder(spec)
        assert len(ladder) == 8
        # Monotonically increasing.
        assert all(b > a for a, b in zip(ladder, ladder[1:]))
        # Endpoints anchored.
        assert ladder[0] == 300.0
        assert ladder[-1] == 420.0

    def test_vanderspoel_not_implemented(self):
        spec = ReplicaTemperatureSpec(
            mode="auto",
            method="vanderspoel",
            T_min=300.0,
            T_max=357.10,
            n_replicas=4,
        )
        with pytest.raises(NotImplementedError, match="vanderspoel|deferred"):
            generate_ladder(spec)


class TestFormatLadderLine:
    def test_basic(self):
        # 3-decimal rendering matches POC's parameters1 line format.
        line = format_ladder_line([300.0, 318.11, 337.11, 357.10])
        assert line == "300.000 318.110 337.110 357.100"

    def test_single_element(self):
        # Edge case: single element still renders.
        assert format_ladder_line([300.0]) == "300.000"

    def test_empty(self):
        assert format_ladder_line([]) == ""


class TestLadderRatios:
    def test_4_replicas(self):
        ladder = [300.0, 318.11, 337.11, 357.10]
        ratios = ladder_ratios(ladder)
        assert len(ratios) == 3
        assert math.isclose(ratios[0], 318.11 / 300.0, rel_tol=1e-9)
        assert math.isclose(ratios[1], 337.11 / 318.11, rel_tol=1e-9)
        assert math.isclose(ratios[2], 357.10 / 337.11, rel_tol=1e-9)

    def test_single(self):
        # Single-element ladder -> no pairs.
        assert ladder_ratios([300.0]) == []
