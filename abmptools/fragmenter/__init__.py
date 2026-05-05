# -*- coding: utf-8 -*-
"""
abmptools.fragmenter
--------------------
FMO 自動フラグメント分割ツール。

PDB から小分子・脂質・ポリマーのフラグメント分割を提案し、
ユーザーが Jupyter UI (A 経路) または CLI ヘッドレス経路 (C 経路) で
確定後、abmptools.log2config 互換の segment_data.dat を出力する。

タンパク質・DNA は対象外 (既存 abmptools.log2config 経由で扱う)。

主要 API:
    FragmenterConfig    全体設定 (JSON 往復可能)
    CutSite             1 結合に対する切断状態 (suggested / enabled)
    MoleculeGroup       canonical SMILES でグループ化された分子集合
    FragmentResult      最終分割結果 (segment_data 互換構造)
    load_pdb_molecules  PDB → 連結成分単位の RDKit Mol list (要 rdkit)
    group_by_smiles     canonical SMILES でグループ化 (要 rdkit)

CLI:
    python -m abmptools.fragmenter {suggest,apply,example}

依存:
    - rdkit-pypi (extras: [fragmenter])
    - ipywidgets / ipykernel (extras: [jupyter]、A 経路使用時のみ)
"""
from .models import (
    FragmenterConfig,
    CutSite,
    MoleculeGroup,
    FragmentResult,
)

# rdkit が無い環境でも (FragmenterConfig 等の dataclass のみ) 部分使用可能
try:
    from .pdb_loader import load_pdb_molecules, LoadedMolecule
    from .grouping import group_by_smiles
    _CORE_AVAILABLE = True
except ImportError:
    _CORE_AVAILABLE = False

__all__ = [
    "FragmenterConfig",
    "CutSite",
    "MoleculeGroup",
    "FragmentResult",
]
if _CORE_AVAILABLE:
    __all__ += ["load_pdb_molecules", "LoadedMolecule", "group_by_smiles"]
