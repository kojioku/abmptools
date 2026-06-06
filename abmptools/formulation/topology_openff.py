"""Topology builder (Phase 1 OpenFF route — Windows native).

複数 species の OpenFF Molecule template と packmol で組まれた mixed PDB から、
OpenFF Interchange を作って GROMACS ``system.top`` / ``system.gro`` を export。

設計:
- amorphous.parameterizer.create_interchange + export_gromacs を呼ぶ薄い wrapper
- ``Topology.from_pdb(mixture_pdb, unique_molecules=[...])`` で全 species を
  1 system に combine、 ``Interchange.from_smirnoff`` で SMIRNOFF apply
- ``charge_from_molecules=[...]`` で NAGL / Gasteiger 事前 charges を再利用
  (Windows native で ``sqm`` 不要)
- TIP3P SMIRNOFF (``tip3p.offxml``) を後ろに stack、 water 部分を上書き

Phase 1 制限:
- water + ions は **mixture PDB 側に既に packmol で入れてある前提** (= tleap
  solvatebox に代わる「pre-solvated packmol」経路)
- OpenMM Modeller.addSolvent 経由の wet solvate は Phase 2 で追加

依存 (全 OS install 可):
- ``openff-toolkit``、 ``openff-interchange``、 ``openmm``
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

logger = logging.getLogger(__name__)


__all__ = [
    "MergedTopologyResult",
    "merge_to_interchange",
    "export_gromacs_files",
]


@dataclass
class MergedTopologyResult:
    """Output of :func:`merge_to_interchange` + :func:`export_gromacs_files`."""

    gro_path: Path
    top_path: Path
    interchange: Any
    n_atoms_total: int


def merge_to_interchange(
    *,
    species_molecules: Sequence[Any],
    species_counts: Sequence[int],
    mixture_pdb: Path,
    box_size_nm: float,
    forcefield_offxmls: Sequence[str] = (
        "openff_unconstrained-2.1.0.offxml",
        "tip3p.offxml",
    ),
    use_precomputed_charges: bool = True,
) -> Any:
    """Merge per-species OpenFF Molecule templates + mixture PDB into 1 Interchange.

    Parameters
    ----------
    species_molecules
        Per-species ``openff.toolkit.Molecule`` の list。 順序は packmol
        input と一致させる ([peptide, enhancer_neu, enhancer_chg, bile_salt,
        ...] の順)。
    species_counts
        対応する ``n_copies``。 informational only (sanity check 用途)。
    mixture_pdb
        packmol で組まれた全 species の mixed PDB。 water/ions が同 packmol
        run で入っている場合も含む。
    box_size_nm
        cubic box edge (nm)。
    forcefield_offxmls
        SMIRNOFF OFFXML の stack。 default = Sage 2.1 + TIP3P (= amorphous
        と同じ慣習)。
    use_precomputed_charges
        ``True`` (default) なら each Molecule の ``partial_charges`` を Interchange
        に伝えて ``sqm`` を呼ばない (Windows native 必須)。

    Returns
    -------
    interchange
        Parameterized ``openff.interchange.Interchange``。
    """
    from ..amorphous.parameterizer import create_interchange

    logger.info(
        "OpenFF route: merge %d species (counts=%s) → Interchange "
        "(FF stack=%s, precomputed_charges=%s)",
        len(species_molecules), list(species_counts), list(forcefield_offxmls),
        use_precomputed_charges,
    )
    interchange = create_interchange(
        molecules=list(species_molecules),
        counts=list(species_counts),
        box_size_nm=float(box_size_nm),
        mixture_pdb=str(mixture_pdb),
        forcefield_name=list(forcefield_offxmls),
        use_precomputed_charges=use_precomputed_charges,
    )
    return interchange


def export_gromacs_files(
    *,
    interchange: Any,
    gro_path: Path,
    top_path: Path,
) -> MergedTopologyResult:
    """``Interchange`` → GROMACS ``system.gro`` / ``system.top``。

    Parameters
    ----------
    interchange
        :func:`merge_to_interchange` の戻り値。
    gro_path, top_path
        出力 path。

    Returns
    -------
    MergedTopologyResult
    """
    from ..amorphous.parameterizer import export_gromacs

    gro_path = Path(gro_path)
    top_path = Path(top_path)
    out = export_gromacs(
        interchange=interchange,
        gro_path=str(gro_path),
        top_path=str(top_path),
    )
    # atom 数を gro から数える (Interchange API に直接 accessor が無い)
    n_atoms = 0
    try:
        with gro_path.open() as f:
            f.readline()  # title
            n_atoms = int(f.readline().strip())
    except (OSError, ValueError):
        n_atoms = 0

    logger.info("OpenFF route: GROMACS export done (%s, %s, %d atoms)",
                out["gro"], out["top"], n_atoms)
    return MergedTopologyResult(
        gro_path=Path(out["gro"]),
        top_path=Path(out["top"]),
        interchange=interchange,
        n_atoms_total=n_atoms,
    )
