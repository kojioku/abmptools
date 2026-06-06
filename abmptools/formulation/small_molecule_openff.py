"""Small-molecule parameterizer (Phase 1 OpenFF route — Windows native).

`abmptools.formulation` の `acpype` / `sqm` / `antechamber` 依存を、 OpenFF
Sage 2.x (SMIRNOFF) + Interchange で置換する **Phase 1 開発中** モジュール。
caprate (Na-C10) と taurocholate を対象。

設計方針 (詳細: `docs/platform_support.md`):

- amorphous.molecule_prep / amorphous.parameterizer の SMIRNOFF パターンを
  formulation の small-molecule stage に流用
- AM1-BCC charges は ``openff-toolkit`` の ``assign_partial_charges('am1bcc')``
  を試みる (内部で `xtb` or 公式 OE backend を呼ぶ可能性あり、 OE が無くても
  動く経路を整える)。 fallback として Gasteiger / 事前計算 charges を許可
- 既存 Amber route (`small_molecule.py`) と並列に存在し、
  ``FormulationBuildConfig.force_field_route = "openff"`` で切替

依存:
- ``openff-toolkit`` (BSD-3)
- ``openff-interchange`` (BSD-3)
- ``rdkit`` (BSD-3)
- (任意) ``xtb-python`` (LGPL-3.0) — AM1-BCC 代替の semiempirical charges

Status: **NOT IMPLEMENTED YET** — Phase 1 skeleton。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from .models import BileSaltSpec, EnhancerSpec


__all__ = [
    "SmallMoleculeOpenFFResult",
    "parameterize_small_mol_openff",
]


@dataclass
class SmallMoleculeOpenFFResult:
    """Output for a single small-molecule species.

    Attributes
    ----------
    pdb_path
        単分子 3D PDB (packmol 入力)。
    interchange
        ``openff.interchange.Interchange`` (SMIRNOFF typed)。
    n_atoms_per_copy
        分子の atom 数。
    net_charge
        formal charge (cation = +1、 anion = -1 等)、 ions balance に使う。
    """

    pdb_path: Path
    interchange: Any
    n_atoms_per_copy: int
    net_charge: int


def parameterize_small_mol_openff(
    *,
    smiles: str,
    resname: str,
    net_charge: int,
    output_dir: Path,
    forcefield_offxml: str = "openff_unconstrained-2.1.0.offxml",
    charge_method: str = "am1bcc",
) -> SmallMoleculeOpenFFResult:
    """SMILES → 3D 配座 → SMIRNOFF typed Interchange.

    amorphous.molecule_prep.smiles_to_molecule + amorphous.parameterizer
    の create_interchange パターンを 1 species 用に再構成。

    Parameters
    ----------
    smiles
        構造 SMILES (charge 情報含む、 例: ``"CCCCCCCCCC(=O)[O-]"``)。
    resname
        GROMACS topology の resname (3-4 文字)、 例: ``"CPN"``、 ``"CPC"``。
    net_charge
        formal charge (acpype の ``-n`` と同義)。
    output_dir
        ``<resname>.pdb`` の出力 dir。
    forcefield_offxml
        SMIRNOFF OFFXML 名。 default は Sage 2.1 unconstrained。
    charge_method
        "am1bcc" (推奨) / "gasteiger" (fallback、 軽量) / "zeros" (debug)。

    Returns
    -------
    SmallMoleculeOpenFFResult

    Raises
    ------
    NotImplementedError
        Phase 1 が未完。
    """
    raise NotImplementedError(
        "parameterize_small_mol_openff is a Phase 1 skeleton; implementation "
        "pending. See docs/platform_support.md for the Phase 1 plan."
    )
