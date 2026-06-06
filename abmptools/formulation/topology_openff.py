"""Topology merger (Phase 1 OpenFF route — Windows native).

複数 species の :class:`openff.interchange.Interchange` を merge し、
TIP3P solvate + Joung-Cheatham ions を加えて、 GROMACS ``system.top`` /
``system.gro`` を export する。

設計方針 (詳細: `docs/platform_support.md`):

- ``Interchange.combine(...)`` で peptide / enhancer / bile salt の
  per-species Interchange を 1 つの system に統合
- Solvation は ``Interchange.set_box(...)`` + TIP3P で。 packmol で water 充填
  しなくても OpenFF 経路では `PDBFile.populate_water()` 経路もある
- Ions balance は OpenFF 経路では ``Interchange.add_water_and_ions(...)`` (将来 API)
  または手動で Na+/Cl- を packmol/parmed 経由

依存:
- ``openff-interchange`` (BSD-3)

Status: **NOT IMPLEMENTED YET** — Phase 1 skeleton。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

__all__ = [
    "MergedTopologyResult",
    "merge_interchanges",
    "solvate_and_neutralize_openff",
]


@dataclass
class MergedTopologyResult:
    """Output of :func:`merge_interchanges` (+ solvate)."""

    gro_path: Path
    top_path: Path
    n_atoms_total: int
    box_nm: float
    counts: Dict[str, int]


def merge_interchanges(
    *,
    species_interchanges: Sequence[Any],
    species_counts: Sequence[int],
    box_size_nm: float,
) -> Any:
    """Combine per-species Interchanges into a system Interchange.

    Parameters
    ----------
    species_interchanges
        Per-species ``Interchange``。 ``[peptide, enhancer_neu, enhancer_chg,
        bile_salt, ...]`` の順。
    species_counts
        対応する ``n_copies``。
    box_size_nm
        cubic box edge (nm)。

    Returns
    -------
    combined
        Merged ``Interchange`` (solvent / ions 追加前)。

    Raises
    ------
    NotImplementedError
        Phase 1 が未完。
    """
    raise NotImplementedError(
        "merge_interchanges is a Phase 1 skeleton; implementation pending. "
        "See docs/platform_support.md for the Phase 1 plan."
    )


def solvate_and_neutralize_openff(
    *,
    combined: Any,
    salt_concentration_M: float = 0.15,
    water_model_offxml: str = "tip3p.offxml",
) -> MergedTopologyResult:
    """TIP3P solvate + Joung-Cheatham Na/Cl で neutral + 0.15 M に。

    Parameters
    ----------
    combined
        ``merge_interchanges`` の戻り値。
    salt_concentration_M
        NaCl 濃度 (Hossain 2023 = 0.15 M)。
    water_model_offxml
        TIP3P SMIRNOFF。 default は OpenFF 提供の ``tip3p.offxml``。

    Returns
    -------
    MergedTopologyResult
        ``system.gro`` / ``system.top`` の path + atom 数 + 各 species count。

    Raises
    ------
    NotImplementedError
        Phase 1 が未完。
    """
    raise NotImplementedError(
        "solvate_and_neutralize_openff is a Phase 1 skeleton; implementation "
        "pending. See docs/platform_support.md for the Phase 1 plan."
    )
