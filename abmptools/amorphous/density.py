# -*- coding: utf-8 -*-
"""
abmptools.amorphous.density
----------------------------
Box-size estimation and weight-fraction → molecule-count conversion.
"""
from __future__ import annotations

import math
from typing import List

AVOGADRO = 6.02214076e23


def weight_fractions_to_counts(
    molecular_weights: List[float],
    weight_fractions: List[float],
    total_molecules: int,
) -> List[int]:
    """Convert weight fractions to integer molecule counts.

    Parameters
    ----------
    molecular_weights : list of float
        Molecular weight of each component [g/mol].
    weight_fractions : list of float
        Target weight fraction for each component (should sum to ~1).
    total_molecules : int
        Desired total number of molecules.

    Returns
    -------
    list of int
        Molecule counts for each component (sums to *total_molecules*).
    """
    if len(molecular_weights) != len(weight_fractions):
        raise ValueError("molecular_weights and weight_fractions must have same length.")
    if any(mw <= 0 for mw in molecular_weights):
        raise ValueError("All molecular weights must be positive.")

    wf_sum = sum(weight_fractions)
    if abs(wf_sum - 1.0) > 0.05:
        raise ValueError(f"Weight fractions should sum to ~1.0, got {wf_sum:.4f}.")
    # normalise
    wf_norm = [w / wf_sum for w in weight_fractions]

    # n_i ∝ wf_i / MW_i
    raw = [wf / mw for wf, mw in zip(wf_norm, molecular_weights)]
    raw_sum = sum(raw)
    fracs = [r / raw_sum * total_molecules for r in raw]

    # round to integers while preserving total
    counts = [int(round(f)) for f in fracs]
    diff = total_molecules - sum(counts)
    if diff != 0:
        # adjust the component with the largest fractional part
        residuals = [(fracs[i] - counts[i], i) for i in range(len(counts))]
        residuals.sort(reverse=(diff > 0))
        for step in range(abs(diff)):
            idx = residuals[step][1]
            counts[idx] += 1 if diff > 0 else -1
    return counts


def estimate_box_size_nm(
    molecular_weights: List[float],
    counts: List[int],
    density_g_cm3: float = 0.8,
) -> float:
    """Estimate cubic box edge length from target density.

    Parameters
    ----------
    molecular_weights : list of float
        Molecular weight of each component [g/mol].
    counts : list of int
        Number of molecules for each component.
    density_g_cm3 : float
        Target mass density [g/cm^3].

    Returns
    -------
    float
        Box edge length [nm].
    """
    total_mass_g = sum(mw * n / AVOGADRO for mw, n in zip(molecular_weights, counts))
    volume_cm3 = total_mass_g / density_g_cm3
    volume_nm3 = volume_cm3 * 1e21  # 1 cm^3 = 1e21 nm^3
    edge_nm = volume_nm3 ** (1.0 / 3.0)
    return edge_nm
