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
from typing import Any, List, Set, Tuple

from .models import CutSite, FragmenterConfig, MoleculeGroup
from .pdb_loader import LoadedMolecule

logger = logging.getLogger(__name__)

# ヘテロ原子 (C, H 以外で隣接除外対象とする原子番号)
HETERO_ELEMENTS = {7, 8, 16, 15, 9, 17, 35, 53}  # N, O, S, P, F, Cl, Br, I


def _check_remaining_mw_ok(
    mol: Any, source_atom: int, exclude: Set[int], min_mw: float,
) -> bool:
    """source_atom から exclude 以外の heavy atom MW 累積が min_mw 以上に達するか。

    early-exit BFS で `min_mw` に達した時点で True を返す。`min_mw <= 0` なら
    常に True (= 制約無効)。
    """
    if min_mw <= 0:
        return True
    cum = 0.0
    visited = set(exclude) | {source_atom}
    queue = deque([source_atom])
    while queue:
        cur = queue.popleft()
        cum += _atom_total_mw(mol.GetAtomWithIdx(cur))
        if cum >= min_mw:
            return True
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            if nb.GetAtomicNum() == 1:
                continue
            if nb.GetIdx() in visited:
                continue
            visited.add(nb.GetIdx())
            queue.append(nb.GetIdx())
    return False


def _enumerate_side_chain_candidate_bonds(
    mol: Any, config: FragmenterConfig,
) -> Set[int]:
    """walk_side_chains 用の relaxed candidate set。

    主鎖 walk と違って、side chain は通常 ester / amine 隣接のため
    `exclude_heteroneighbor` で全 bond が除外されてしまう。side chain walk では
    その filter を **無視** し、C-X (X = N/O/S/...) 単結合も許可する
    (ABINIT-MP 慣習 C-side BDA に基づき、`_walk_side_chains_from_main` 内で
    BDA/BAA を決定する)。環内 / 多重結合除外は維持。
    """
    from rdkit import Chem

    candidates: Set[int] = set()
    for bond in mol.GetBonds():
        if config.exclude_multibond and bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        z1, z2 = a1.GetAtomicNum(), a2.GetAtomicNum()
        is_cc = z1 == 6 and z2 == 6
        is_cx = (
            (z1 == 6 and z2 in HETERO_ELEMENTS)
            or (z2 == 6 and z1 in HETERO_ELEMENTS)
        )
        if not (is_cc or is_cx):
            continue
        if config.exclude_ring_cc and bond.IsInRing():
            continue
        # ヘテロ隣接フィルタは無視 (side chain walk 用)
        candidates.add(bond.GetIdx())
    return candidates


def _walk_side_chains_from_main(
    mol: Any,
    main_path: List[int],
    candidate_bonds: Set[int],
    target_mw: float,
    min_terminal_ratio: float = 0.5,
) -> List[Tuple[int, int, int]]:
    """主鎖 main_path の各 atom から side chain を BFS walk して cut を提案。

    Returns
    -------
    List[Tuple[bond_idx, bda_atom, baa_atom]]
        side chain 内で target_mw 累積に達した bond を cut として提案する。
        BDA = walk 起点側 (= main atom 側)、BAA = walk 進行方向側 (= side chain
        の terminal 側)。C-X bond なら C 側を BDA に固定。
    """
    main_set = set(main_path)
    selected: List[Tuple[int, int, int]] = []

    for main_atom in main_path:
        # main atom から伸びる各 side chain を独立に BFS walk
        for nb in mol.GetAtomWithIdx(main_atom).GetNeighbors():
            if nb.GetAtomicNum() == 1 or nb.GetIdx() in main_set:
                continue
            start = nb.GetIdx()
            visited = set(main_set) | {start}
            # queue 各要素: (current atom, cum_mw)
            queue = deque([(start, _atom_total_mw(mol.GetAtomWithIdx(start)))])
            while queue:
                cur, cum = queue.popleft()
                cur_atom = mol.GetAtomWithIdx(cur)
                # 次の heavy 隣接 (main / visited を除く) を探索
                nbrs = sorted([
                    n.GetIdx() for n in cur_atom.GetNeighbors()
                    if n.GetAtomicNum() != 1 and n.GetIdx() not in visited
                ])
                for nb_idx in nbrs:
                    if nb_idx in visited:
                        continue
                    visited.add(nb_idx)
                    bond = mol.GetBondBetweenAtoms(cur, nb_idx)
                    if bond is None:
                        continue
                    new_cum = cum + _atom_total_mw(mol.GetAtomWithIdx(nb_idx))
                    bond_idx = bond.GetIdx()
                    cut_here = (bond_idx in candidate_bonds and new_cum >= target_mw)
                    if cut_here and min_terminal_ratio > 0:
                        # 進行方向 (nb_idx 側) の残り fragment が
                        # target_mw * ratio 未満なら cut 抑制
                        ok = _check_remaining_mw_ok(
                            mol, nb_idx, visited - {nb_idx},
                            target_mw * min_terminal_ratio,
                        )
                        if not ok:
                            cut_here = False
                    if cut_here:
                        # BDA / BAA: C-X case は C 側、それ以外は walk 起点側 (cur)
                        a_cur = mol.GetAtomWithIdx(cur)
                        a_nb = mol.GetAtomWithIdx(nb_idx)
                        zc, zn = a_cur.GetAtomicNum(), a_nb.GetAtomicNum()
                        if zc == 6 and zn != 6:
                            bda, baa = cur, nb_idx
                        elif zn == 6 and zc != 6:
                            bda, baa = nb_idx, cur
                        else:
                            bda, baa = cur, nb_idx  # walk 起点側
                        selected.append((bond_idx, bda, baa))
                        queue.append((nb_idx, 0.0))
                    else:
                        queue.append((nb_idx, new_cum))
    return selected


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

    cut_results = _select_cuts_along_path(
        mol, main_chain, candidate, config.target_mw,
        min_terminal_ratio=config.min_terminal_fragment_ratio,
    )

    # walk_side_chains=True なら主鎖 walk 後に各 main atom の side chain も walk
    # side chain 用に candidate を別途構築 (ester/amine 隣接でも C-C や C-X を許可)
    if config.walk_side_chains:
        side_candidate = _enumerate_side_chain_candidate_bonds(mol, config)
        cut_results = cut_results + _walk_side_chains_from_main(
            mol, main_chain, side_candidate, config.target_mw,
            min_terminal_ratio=config.min_terminal_fragment_ratio,
        )
        # 同じ bond が両 walk で追加されないよう dedupe
        seen_bonds: Set[int] = set()
        dedup: List[Tuple[int, int, int]] = []
        for bond_idx, bda, baa in cut_results:
            if bond_idx in seen_bonds:
                continue
            seen_bonds.add(bond_idx)
            dedup.append((bond_idx, bda, baa))
        cut_results = dedup

    cuts: List[CutSite] = []
    for bond_idx, bda_atom, baa_atom in cut_results:
        bond = mol.GetBondWithIdx(bond_idx)
        cuts.append(CutSite(
            bond_idx=bond_idx,
            atom1_idx=bond.GetBeginAtomIdx(),
            atom2_idx=bond.GetEndAtomIdx(),
            suggested=True,
            enabled=True,
            bda_atom_idx=bda_atom,
            baa_atom_idx=baa_atom,
        ))
    return cuts


def decide_bda_baa_for_manual_cut(
    mol: Any,
    atom1_idx: int,
    atom2_idx: int,
) -> Tuple[int, int]:
    """ユーザーが手動で追加した cut bond の BDA/BAA を decide する。

    Rules (auto_split.suggest_cuts と整合):
        - C-X (X = N/O/S/...): C 側を BDA (ABINIT-MP 慣習、ユーザー指定)
        - C-C: atom_idx 若い側を BDA (deterministic な default)

    Returns
    -------
    (bda_atom_idx, baa_atom_idx)
    """
    a1 = mol.GetAtomWithIdx(atom1_idx)
    a2 = mol.GetAtomWithIdx(atom2_idx)
    z1, z2 = a1.GetAtomicNum(), a2.GetAtomicNum()
    if z1 == 6 and z2 != 6:
        return atom1_idx, atom2_idx
    if z2 == 6 and z1 != 6:
        return atom2_idx, atom1_idx
    # C-C もしくは X-X: atom_idx 若い側を BDA
    if atom1_idx <= atom2_idx:
        return atom1_idx, atom2_idx
    return atom2_idx, atom1_idx


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
    """切断対象の単結合を列挙する。

    デフォルトは C-C 単結合のみ。``config.include_c_heteroatom=True`` の場合
    は C-X 単結合 (X = N/O/S/P/F/Cl/Br/I) も candidate に含める。
    ``exclude_heteroneighbor`` フィルタは **C-C bond にのみ適用** され、
    C-X bond は別途許可される。
    """
    from rdkit import Chem

    candidates: Set[int] = set()
    for bond in mol.GetBonds():
        if config.exclude_multibond and bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        z1, z2 = a1.GetAtomicNum(), a2.GetAtomicNum()

        is_cc = z1 == 6 and z2 == 6
        is_cx = (
            (z1 == 6 and z2 in HETERO_ELEMENTS)
            or (z2 == 6 and z1 in HETERO_ELEMENTS)
        )
        if not is_cc and not (is_cx and config.include_c_heteroatom):
            continue

        if config.exclude_ring_cc and bond.IsInRing():
            continue

        # ヘテロ隣接除外は C-C bond にのみ適用
        if is_cc and config.exclude_heteroneighbor:
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
    """atom 自身の質量 + 結合している H の質量 (implicit + explicit 両方)。

    PDB 由来の Mol は explicit H を持つ (RDKit が PDB から読み込むと H atom が
    別 atom として登録される) ため、`GetTotalNumHs(includeNeighbors=False)` は
    implicit H 0 を返してしまう。結果として CH3 の MW が ~12 (本来 ~15) と
    低くなり、累積 walk の cut 位置がずれる。implicit + explicit H neighbors の
    両方を加算することで、heavy-only Mol と H 込み Mol の挙動を一致させる。
    """
    mw = atom.GetMass()
    n_h_implicit = atom.GetTotalNumHs(includeNeighbors=False)
    n_h_explicit = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1)
    mw += (n_h_implicit + n_h_explicit) * 1.008
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
    min_terminal_ratio: float = 0.5,
) -> List[Tuple[int, int, int]]:
    """主鎖 path を辿りつつ、累積 MW ≥ target_mw のたびに candidate bond で切断。

    Returns
    -------
    List[Tuple[bond_idx, bda_atom_idx, baa_atom_idx]]
        BDA / BAA 役割は walk 方向で決定する: walk 起点側 (path[0] 側) の atom
        が BDA、反対側 (path[-1] 側) の atom が BAA。これにより「main chain
        起点側 fragment が electron pair を保持」する形になる (アミノ酸の
        N→C 方向のアナロジーを generic C-C cut に拡張)。
    """
    if len(path) < 2:
        return []
    selected: List[Tuple[int, int, int]] = []
    main_chain_set = set(path)

    # Pre-compute: 各 path 位置以降の残り MW (main + side chain)。
    # cut 判定で「進行方向 fragment が target_mw * ratio 未満なら抑制」に使う。
    remaining_mw = [0.0] * (len(path) + 1)
    for i in range(len(path) - 1, -1, -1):
        remaining_mw[i] = (
            remaining_mw[i + 1]
            + _atom_total_mw(mol.GetAtomWithIdx(path[i]))
            + _sidechain_mw(mol, path[i], main_chain_set)
        )

    cum_mw = 0.0
    min_terminal_mw = max(0.0, target_mw * min_terminal_ratio)
    for i in range(len(path) - 1):
        atom = mol.GetAtomWithIdx(path[i])
        cum_mw += _atom_total_mw(atom)
        cum_mw += _sidechain_mw(mol, path[i], main_chain_set)

        bond = mol.GetBondBetweenAtoms(path[i], path[i + 1])
        if bond is None:
            continue
        bond_idx = bond.GetIdx()
        if (
            bond_idx in candidate_bonds
            and cum_mw >= target_mw
            and remaining_mw[i + 1] >= min_terminal_mw
        ):
            # BDA / BAA 役割の決定:
            #   - C-X (X = N/O/S/...): C 側を BDA、X 側を BAA (ABINIT-MP 慣習)
            #   - C-C: walk 起点側 (path[i]) を BDA、進行方向側 (path[i+1]) を BAA
            a_i = mol.GetAtomWithIdx(path[i])
            a_j = mol.GetAtomWithIdx(path[i + 1])
            zi, zj = a_i.GetAtomicNum(), a_j.GetAtomicNum()
            if zi == 6 and zj != 6:
                bda, baa = path[i], path[i + 1]      # i は C
            elif zj == 6 and zi != 6:
                bda, baa = path[i + 1], path[i]      # j は C
            else:
                bda, baa = path[i], path[i + 1]      # C-C: walk 起点側
            selected.append((bond_idx, bda, baa))
            cum_mw = 0.0
    return selected
