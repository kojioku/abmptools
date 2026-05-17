"""
hbond_detector.py
-----------------
Geometric H-bond detection with orthogonal PBC minimum-image convention.

Criteria (selectable):
- 'luzar-chandler' (default):  d(D...A) <= 3.5 Å  AND  angle(D-H...A) >= 120°
- 'strict':                    d(H...A) <= 2.5 Å  AND  angle(D-H...A) >= 150°
- 'custom':                    user-specified thresholds via HBondCriteria dataclass

References:
- A. Luzar and D. Chandler, Nature 379, 55-57 (1996)
- F.W. Starr et al., J. Chem. Phys. 113, 9727 (2000)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from .bdf_reader import CellBox


@dataclass
class HBondCriteria:
    """Geometric thresholds for H-bond detection."""
    d_da_max: float = 3.5      # Å, donor-acceptor (heavy atom) distance max
    d_ha_max: Optional[float] = None  # Å, H-acceptor distance max (optional)
    angle_min: float = 120.0   # degrees, ∠(D-H...A) min

    @classmethod
    def luzar_chandler(cls) -> "HBondCriteria":
        return cls(d_da_max=3.5, d_ha_max=None, angle_min=120.0)

    @classmethod
    def strict(cls) -> "HBondCriteria":
        return cls(d_da_max=3.5, d_ha_max=2.5, angle_min=150.0)


@dataclass
class HBond:
    """One detected H-bond."""
    donor_mol: int        # molecule index of donor
    donor_d: int          # local atom index of D (heavy atom, e.g. carboxyl oh)
    donor_h: int          # local atom index of H (e.g. carboxyl ho)
    acceptor_mol: int     # molecule index of acceptor
    acceptor_a: int       # local atom index of A (heavy atom)
    d_da: float           # Å
    d_ha: float           # Å
    angle: float          # degrees


def minimum_image_vector(
    dr: np.ndarray, cell: CellBox
) -> np.ndarray:
    """Apply orthogonal PBC minimum-image convention.

    dr: (..., 3) array of raw distance vectors
    cell: CellBox(a, b, c)
    Returns shortest-image vector(s).
    """
    box = cell.as_array()
    return dr - np.round(dr / box) * box


def _angle_deg(d_pos: np.ndarray, h_pos: np.ndarray, a_pos: np.ndarray,
               cell: CellBox) -> float:
    """Compute D-H...A angle in degrees (apex at H)."""
    dh = minimum_image_vector(d_pos - h_pos, cell)
    ah = minimum_image_vector(a_pos - h_pos, cell)
    norm_dh = np.linalg.norm(dh)
    norm_ah = np.linalg.norm(ah)
    if norm_dh < 1e-9 or norm_ah < 1e-9:
        return 0.0
    cos_th = np.clip(np.dot(dh, ah) / (norm_dh * norm_ah), -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_th)))


@dataclass
class DonorSite:
    """A donor: (mol_index, D_local, H_local). E.g. for COOH: (i, oh, ho)."""
    mol_index: int
    d_local: int
    h_local: int


@dataclass
class AcceptorSite:
    """An acceptor: (mol_index, A_local). E.g. for amide C=O: (j, o)."""
    mol_index: int
    a_local: int


def detect_hbonds(
    donors: List[DonorSite],
    acceptors: List[AcceptorSite],
    positions: List[np.ndarray],
    cell: CellBox,
    criteria: HBondCriteria = None,
    allow_self: bool = False,
) -> List[HBond]:
    """Detect H-bonds between donor and acceptor sites.

    Parameters
    ----------
    donors : list of DonorSite
    acceptors : list of AcceptorSite
    positions : list of (n_atoms, 3) arrays, one per molecule
    cell : CellBox (orthogonal)
    criteria : HBondCriteria (default: Luzar-Chandler)
    allow_self : if False (default), skip donor-acceptor pairs on same molecule

    Returns
    -------
    list of HBond
    """
    if criteria is None:
        criteria = HBondCriteria.luzar_chandler()

    out: List[HBond] = []
    d_da_max = criteria.d_da_max
    d_ha_max = criteria.d_ha_max
    angle_min = criteria.angle_min

    for don in donors:
        d_pos = positions[don.mol_index][don.d_local]
        h_pos = positions[don.mol_index][don.h_local]
        for acc in acceptors:
            if not allow_self and acc.mol_index == don.mol_index:
                continue
            a_pos = positions[acc.mol_index][acc.a_local]
            dr_da = minimum_image_vector(a_pos - d_pos, cell)
            d_da = float(np.linalg.norm(dr_da))
            if d_da > d_da_max:
                continue
            dr_ha = minimum_image_vector(a_pos - h_pos, cell)
            d_ha = float(np.linalg.norm(dr_ha))
            if d_ha_max is not None and d_ha > d_ha_max:
                continue
            angle = _angle_deg(d_pos, h_pos, a_pos, cell)
            if angle < angle_min:
                continue
            out.append(HBond(
                donor_mol=don.mol_index,
                donor_d=don.d_local,
                donor_h=don.h_local,
                acceptor_mol=acc.mol_index,
                acceptor_a=acc.a_local,
                d_da=d_da,
                d_ha=d_ha,
                angle=angle,
            ))
    return out
