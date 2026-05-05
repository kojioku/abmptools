# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.replica_temperatures
--------------------------------------------
Temperature-ladder generation for gREST replicas.

Two modes (cf. :class:`RESTSelectionSpec`):

- ``auto``: ``method="geometric"`` -- ``T_i = T_min * (T_max/T_min)**(i/(N-1))``.
  ``method="vanderspoel"`` (Patriksson-van der Spoel target-acceptance
  formula) raises :class:`NotImplementedError` in v1.20.0.
- ``manual``: passthrough of the user-supplied list.

The returned list is rounded to two decimal places to match the
GENESIS ``parameters1 = `` printout (e.g. POC ``300.00 318.11 337.11
357.10``). Use :func:`format_ladder_line` to render the
space-separated ladder string for the ``[REMD]`` block.
"""
from __future__ import annotations

import logging
from typing import List

from .models import ReplicaTemperatureSpec

logger = logging.getLogger(__name__)


def generate_ladder(spec: ReplicaTemperatureSpec) -> List[float]:
    """Return the explicit temperature ladder for ``[REMD] parameters1``.

    Parameters
    ----------
    spec
        Validated :class:`ReplicaTemperatureSpec`.

    Returns
    -------
    List[float]
        Length ``spec.n_replicas`` (or ``len(spec.temperatures)`` in
        manual mode), rounded to 2 decimals.
    """
    if spec.mode == "manual":
        return [round(float(t), 2) for t in spec.temperatures]

    if spec.method == "vanderspoel":
        raise NotImplementedError(
            "Patriksson-van der Spoel acceptance-ratio formula is "
            "deferred to v1.20.x. Use mode='manual' or "
            "method='geometric' in v1.20.0."
        )

    if spec.method != "geometric":
        # Defensive: __post_init__ already rejects, but guard anyway.
        raise ValueError(f"Unknown method: {spec.method}")

    # Geometric series: T_i = T_min * (T_max / T_min)**(i / (N-1)).
    n = spec.n_replicas
    ratio = spec.T_max / spec.T_min
    ladder = [
        round(spec.T_min * (ratio ** (i / (n - 1))), 2) for i in range(n)
    ]
    # Ensure last entry exactly matches T_max after rounding (avoid drift).
    ladder[-1] = round(spec.T_max, 2)
    return ladder


def format_ladder_line(ladder: List[float]) -> str:
    """Render the ladder for use as ``parameters1 = ...`` in ``.inp``.

    Format mirrors GENESIS POC: each value with three decimals,
    separated by single spaces (e.g. ``"300.000 318.110 337.110 357.100"``).
    """
    return " ".join(f"{t:.3f}" for t in ladder)


def ladder_ratios(ladder: List[float]) -> List[float]:
    """Return the per-pair temperature ratios ``T_{i+1}/T_i``.

    Used by :mod:`forcefield_check` to flag ladders likely to give
    poor acceptance (rule of thumb: ratios in 1.05--1.07 give >= 0.2
    acceptance for proteins-in-water).
    """
    return [ladder[i + 1] / ladder[i] for i in range(len(ladder) - 1)]
