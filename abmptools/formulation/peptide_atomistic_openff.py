"""Peptide builder (Phase 1 OpenFF route — Windows native compatible).

`abmptools.formulation` の `tleap` 依存 peptide stage を、 OpenFF Toolkit +
`openff-amber-ff-ports` (ff14SB SMIRNOFF port) で置換する **Phase 1 開発中**
モジュール。

設計方針 (詳細: `docs/platform_support.md` Phase 1 計画節):

- Linux/macOS/Windows のどの OS でも install できる科学スタックのみを使用
  (AmberTools が install できない Windows native 環境への配布用)
- amorphous の OpenFF 経路 (`amorphous.molecule_prep` + `amorphous.parameterizer`)
  と同じパターン
- 既存 Amber route (`peptide_atomistic.py`) と並列に存在し、
  `FormulationBuildConfig.force_field_route = "openff"` で切替

依存:
- ``openff-toolkit`` (BSD-3)
- ``openff-amber-ff-ports`` (BSD-3、 ff14SB SMIRNOFF 提供)
- ``openff-interchange`` (BSD-3、 GROMACS export)
- ``rdkit`` (BSD-3、 3D 配座生成)

Status: **NOT IMPLEMENTED YET** — Phase 1 skeleton。 関数 signature 凍結、
本体は次フェーズで実装。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple

from .models import PeptideSpec


__all__ = [
    "PeptideOpenFFResult",
    "build_peptide_openff",
    "check_openff_amber_ff_ports_available",
]


@dataclass
class PeptideOpenFFResult:
    """Output of :func:`build_peptide_openff`.

    Mirrors the relevant subset of the Amber-route result so that
    :mod:`abmptools.formulation.topology` can merge regardless of route.

    Attributes
    ----------
    pdb_path
        単一 peptide の 3D 構造 PDB (packmol 入力)。
    interchange
        ``openff.interchange.Interchange`` インスタンス (forcefield + topology
        adopted)。 後段の :func:`abmptools.formulation.topology_openff.merge_interchanges`
        で全 component を combine する。
    n_atoms_per_copy
        単一 peptide の atom 数 (packmol fixed や cluster_pdb 計算用)。
    """

    pdb_path: Path
    interchange: Any
    n_atoms_per_copy: int


def check_openff_amber_ff_ports_available() -> None:
    """Raise RuntimeError if ``openff-amber-ff-ports`` is not installed.

    `openff-amber-ff-ports` は ff14SB を SMIRNOFF 形式で提供する補助
    パッケージ。 これが無いと OpenFF route で peptide build できない。
    """
    try:
        import openff.amber_ff_ports  # noqa: F401
    except ImportError as e:
        raise RuntimeError(
            "openff-amber-ff-ports is required for the OpenFF peptide route.\n"
            "Install: pip install openff-amber-ff-ports  (or conda-forge)"
        ) from e


def build_peptide_openff(
    *,
    spec: PeptideSpec,
    output_dir: Path,
    forcefield_offxml: str = "openff-amber-ff-14SB.offxml",
    water_model: str = "tip3p",
) -> PeptideOpenFFResult:
    """Build a single peptide topology + 3D PDB via OpenFF / amber-ff-ports.

    Parameters
    ----------
    spec
        Peptide specification (sequence, caps, disulfide pairs).
    output_dir
        ``build/input/`` 等の出力 dir。 ``<name>.pdb`` を書き出す。
    forcefield_offxml
        SMIRNOFF OFFXML 名。 default は ff14SB port。
    water_model
        Disulfide / capping 周辺の H 配置で使う water model 名 (将来用)。

    Returns
    -------
    PeptideOpenFFResult
        PDB path + Interchange (force-field-typed topology) + atom 数。

    Raises
    ------
    NotImplementedError
        Phase 1 が未完。 本実装は次フェーズで追加。
    RuntimeError
        ``openff-toolkit`` or ``openff-amber-ff-ports`` 未 install。
    """
    check_openff_amber_ff_ports_available()
    raise NotImplementedError(
        "build_peptide_openff is a Phase 1 skeleton; implementation pending. "
        "See docs/platform_support.md for the Phase 1 plan and progress."
    )
