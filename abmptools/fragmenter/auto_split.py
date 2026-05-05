# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.auto_split
-------------------------------
C-C 単結合切断候補の自動提案。

アルゴリズム:
    1. 全 bond を走査し、両端 C / SINGLE / 環内除外 / 多重結合除外 /
       ヘテロ隣接除外でフィルタした candidate set を作る
    2. heavy atom only グラフでの直径 (graph diameter) を求めて主鎖と定義
    3. 主鎖を端から walk しつつ、各 atom の MW (側鎖 heavy atom MW を加算)
       を累積、accumulator が target_mw を超えたら次の candidate bond で切断

非線形 (分岐) 分子の場合は graph diameter を主鎖として、それ以外の側鎖は
主鎖 atom のローカルな MW として加算する。複数の鎖端を持つ T 字 / Y 字状
分子では、主鎖 = 最長 shortest path を選ぶ。
"""
from __future__ import annotations

import logging
from collections import deque
from typing import Any, List, Set

from .models import CutSite, FragmenterConfig, MoleculeGroup
from .pdb_loader import LoadedMolecule

logger = logging.getLogger(__name__)

# ヘテロ原子 (C, H 以外で隣接除外対象とする原子番号)
HETERO_ELEMENTS = {7, 8, 16, 15, 9, 17, 35, 53}  # N, O, S, P, F, Cl, Br, I


def suggest_cuts(mol: Any, config: FragmenterConfig) -> List[CutSite]:
    """1 分子に対する切断候補を提案する。

    Parameters
    ----------
    mol
        RDKit Mol (H 含む or 含まない、どちらでも動く)。
    config
        FragmenterConfig (target_mw, exclude_* フィルタ参照)。

    Returns
    -------
    List[CutSite]
        suggested=True, enabled=True で生成された切断点リスト。
    """
    candidate = _enumerate_candidate_bonds(mol, config)
    if not candidate:
        return []

    main_chain = _find_main_chain_heavy(mol)
    if len(main_chain) < 2:
        return []

    cut_bond_indices = _select_cuts_along_path(mol, main_chain, candidate, config.target_mw)

    cuts: List[CutSite] = []
    for bond_idx in cut_bond_indices:
        bond = mol.GetBondWithIdx(bond_idx)
        cuts.append(CutSite(
            bond_idx=bond_idx,
            atom1_idx=bond.GetBeginAtomIdx(),
            atom2_idx=bond.GetEndAtomIdx(),
            suggested=True,
            enabled=True,
        ))
    return cuts


def suggest_cuts_for_groups(
    groups: List[MoleculeGroup],
    molecules: List[LoadedMolecule],
    config: FragmenterConfig,
) -> None:
    """各 group の代表分子に対して切断候補を提案し、group.cut_sites を上書きする。

    in-place で更新。
    """
    for g in groups:
        rep = molecules[g.representative_mol_idx]
        g.cut_sites = suggest_cuts(rep.mol, config)
        logger.info(
            "Suggested %d cut(s) for group SMILES=%s (n_copies=%d)",
            len(g.cut_sites), g.smiles, g.n_copies,
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _enumerate_candidate_bonds(mol: Any, config: FragmenterConfig) -> Set[int]:
    """切断対象 C-C 単結合を列挙する。"""
    from rdkit import Chem

    candidates: Set[int] = set()
    for bond in mol.GetBonds():
        if config.exclude_multibond and bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
            continue
        if config.exclude_ring_cc and bond.IsInRing():
            continue
        if config.exclude_heteroneighbor:
            if _has_hetero_neighbor(a1, exclude_idx=a2.GetIdx()) or \
               _has_hetero_neighbor(a2, exclude_idx=a1.GetIdx()):
                continue
        candidates.add(bond.GetIdx())
    return candidates


def _has_hetero_neighbor(atom: Any, exclude_idx: int) -> bool:
    """atom の隣接 (exclude_idx 以外) にヘテロ原子が含まれるか。"""
    for nb in atom.GetNeighbors():
        if nb.GetIdx() == exclude_idx:
            continue
        if nb.GetAtomicNum() in HETERO_ELEMENTS:
            return True
    return False


def _heavy_neighbors(atom: Any) -> List[int]:
    """atom の heavy atom (H 以外) 隣接 atom idx list (ソート済)。"""
    return sorted([
        nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1
    ])


def _bfs_path(mol: Any, src: int, dst: int) -> List[int]:
    """heavy atom only graph での src → dst の shortest path。

    Returns
    -------
    List[int]
        path 上の atom idx list (両端を含む)。到達不可なら空 list。
    """
    if src == dst:
        return [src]
    visited = {src}
    parent = {src: -1}
    queue = deque([src])
    while queue:
        cur = queue.popleft()
        if cur == dst:
            break
        for nb_idx in _heavy_neighbors(mol.GetAtomWithIdx(cur)):
            if nb_idx in visited:
                continue
            visited.add(nb_idx)
            parent[nb_idx] = cur
            queue.append(nb_idx)
    if dst not in parent:
        return []
    # 復元
    path: List[int] = []
    cur = dst
    while cur != -1:
        path.append(cur)
        cur = parent[cur]
    path.reverse()
    return path


def _bfs_farthest(mol: Any, src: int) -> int:
    """heavy atom only graph で src から最遠の atom idx を返す。"""
    visited = {src}
    queue = deque([src])
    far = src
    while queue:
        cur = queue.popleft()
        far = cur  # 最後に visit した atom が最遠
        for nb_idx in _heavy_neighbors(mol.GetAtomWithIdx(cur)):
            if nb_idx in visited:
                continue
            visited.add(nb_idx)
            queue.append(nb_idx)
    return far


def _find_main_chain_heavy(mol: Any) -> List[int]:
    """heavy atom only グラフの直径 (graph diameter) を主鎖として返す。

    2-pass BFS algorithm (tree や近似的な graph で diameter を求める標準手法):
        1. 任意 atom から最遠 far1 を見つける
        2. far1 から最遠 far2 を見つける
        3. far1 ↔ far2 の shortest path が直径
    """
    heavy_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    if len(heavy_atoms) < 2:
        return list(heavy_atoms)
    far1 = _bfs_farthest(mol, heavy_atoms[0])
    far2 = _bfs_farthest(mol, far1)
    return _bfs_path(mol, far1, far2)


def _atom_total_mw(atom: Any) -> float:
    """atom 自身の質量 + 結合している H の質量。"""
    mw = atom.GetMass()
    n_h = atom.GetTotalNumHs(includeNeighbors=False)
    mw += n_h * 1.008
    return mw


def _sidechain_mw(mol: Any, anchor: int, main_chain_set: Set[int]) -> float:
    """anchor atom から主鎖外の側鎖 heavy atom 全部の MW を集計する。"""
    visited = {anchor}
    total = 0.0
    queue: deque = deque()
    for nb_idx in _heavy_neighbors(mol.GetAtomWithIdx(anchor)):
        if nb_idx in main_chain_set:
            continue
        queue.append(nb_idx)
    while queue:
        cur = queue.popleft()
        if cur in visited or cur in main_chain_set:
            continue
        visited.add(cur)
        total += _atom_total_mw(mol.GetAtomWithIdx(cur))
        for nb_idx in _heavy_neighbors(mol.GetAtomWithIdx(cur)):
            if nb_idx in visited or nb_idx in main_chain_set:
                continue
            queue.append(nb_idx)
    return total


def _select_cuts_along_path(
    mol: Any,
    path: List[int],
    candidate_bonds: Set[int],
    target_mw: float,
) -> List[int]:
    """主鎖 path を辿りつつ、累積 MW ≥ target_mw のたびに candidate bond で切断。"""
    if len(path) < 2:
        return []
    selected: List[int] = []
    main_chain_set = set(path)
    cum_mw = 0.0
    for i in range(len(path) - 1):
        atom = mol.GetAtomWithIdx(path[i])
        cum_mw += _atom_total_mw(atom)
        cum_mw += _sidechain_mw(mol, path[i], main_chain_set)

        bond = mol.GetBondBetweenAtoms(path[i], path[i + 1])
        if bond is None:
            continue
        bond_idx = bond.GetIdx()
        if bond_idx in candidate_bonds and cum_mw >= target_mw:
            selected.append(bond_idx)
            cum_mw = 0.0
    return selected
