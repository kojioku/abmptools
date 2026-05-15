# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.chain_splitter
------------------------------------------------
ring 外 chain atoms を connected component に分けて、各 component を
target_mw 基準で切断する。

`abmptools.fragmenter.auto_split._atom_total_mw` (private import) を流用し、
MW walk のロジックを CG segment 用に応用する。
"""
from __future__ import annotations

import logging
from collections import deque
from typing import Any, List, Set

from ..auto_split import _atom_total_mw
from .models import Segment

logger = logging.getLogger(__name__)


def split_chain(
    mol: Any,
    used_atoms_in_rings: Set[int],
    target_mw: float = 200.0,
    start_segment_id: int = 0,
) -> List[Segment]:
    """ring 外 chain atoms を connected component に分けて、各 component を
    target_mw 基準で切断する。

    Parameters
    ----------
    mol
        RDKit Mol。
    used_atoms_in_rings
        ring segment に既に含まれている atom idx set (chain から除外)。
    target_mw
        chain 切断の目安 MW (g/mol)。
    start_segment_id
        生成する Segment の segment_id 開始値。

    Returns
    -------
    List[Segment]
        chain segment のリスト (kind="chain")。
    """
    chain_atoms = {
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetAtomicNum() != 1 and atom.GetIdx() not in used_atoms_in_rings
    }

    if not chain_atoms:
        return []

    # connected component (chain 内のみ、ring atom 経由は跨がない)
    components = _connected_components_within(mol, chain_atoms)

    segments: List[Segment] = []
    sid = start_segment_id
    for comp in components:
        for sub in _split_component_by_mw(mol, comp, target_mw):
            segments.append(Segment(
                segment_id=sid,
                atom_indices=sorted(sub),
                kind="chain",
            ))
            sid += 1

    logger.info(
        "split_chain: %d chain segment(s), %d chain atoms total",
        len(segments), len(chain_atoms),
    )
    return segments


def _connected_components_within(mol: Any, atom_set: Set[int]) -> List[List[int]]:
    """atom_set 内の atoms について、その内部のみで連結成分を列挙する。"""
    visited: Set[int] = set()
    components: List[List[int]] = []
    for start in sorted(atom_set):
        if start in visited:
            continue
        comp: List[int] = []
        queue = deque([start])
        while queue:
            cur = queue.popleft()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            atom = mol.GetAtomWithIdx(cur)
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 1:
                    continue
                nb_idx = nb.GetIdx()
                if nb_idx in atom_set and nb_idx not in visited:
                    queue.append(nb_idx)
        components.append(sorted(comp))
    return components


def _split_component_by_mw(
    mol: Any,
    comp_atoms: List[int],
    target_mw: float,
) -> List[List[int]]:
    """1 つの chain component を target_mw 基準で分割する。

    Strategy: 端 atom (heavy_in_comp = 1) から greedy walk、累積 MW ≥
    target_mw のたびに次 atom を別 segment へ。
    """
    if not comp_atoms:
        return []
    if target_mw <= 0:
        return [list(comp_atoms)]

    comp_set = set(comp_atoms)
    endpoints = []
    for a in comp_atoms:
        atom = mol.GetAtomWithIdx(a)
        heavy_in_comp = sum(
            1 for nb in atom.GetNeighbors()
            if nb.GetAtomicNum() != 1 and nb.GetIdx() in comp_set
        )
        if heavy_in_comp <= 1:
            endpoints.append(a)

    if not endpoints:
        # 環状 chain (chain 内に内部 ring; 通常起きない) → 全体 1 segment
        return [list(comp_atoms)]

    start = sorted(endpoints)[0]
    sub_segments: List[List[int]] = [[]]
    cum_mw = 0.0
    visited: Set[int] = set()
    # (current_atom, prev_atom) — DFS で chain を辿る
    stack = [(start, -1)]
    while stack:
        cur, prev = stack.pop()
        if cur in visited:
            continue
        visited.add(cur)
        atom = mol.GetAtomWithIdx(cur)
        sub_segments[-1].append(cur)
        cum_mw += _atom_total_mw(atom)

        heavy_neighbors = sorted(
            [
                nb.GetIdx() for nb in atom.GetNeighbors()
                if nb.GetAtomicNum() != 1
                and nb.GetIdx() in comp_set
                and nb.GetIdx() != prev
                and nb.GetIdx() not in visited
            ],
            reverse=True,  # DFS で若い idx を後で訪問 (sub_segments の順序を deterministic に)
        )

        if cum_mw >= target_mw and heavy_neighbors:
            sub_segments.append([])
            cum_mw = 0.0

        for nb_idx in heavy_neighbors:
            stack.append((nb_idx, cur))

    return [s for s in sub_segments if s]
