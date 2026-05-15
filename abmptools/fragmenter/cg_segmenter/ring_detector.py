# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.ring_detector
-----------------------------------------------
RDKit `GetRingInfo().AtomRings()` を使った ring segment 抽出。

- 各 SSSR ring を 1 Segment とする
- fused ring (共有 atom あり) は共有 atom を **両方の Segment に含める**
  (`allow_atom_sharing=True` の場合) — atom 共有許容が CG 用途では一般的
- ring に直接 attach する 1 heavy atom 置換基 (-OH, -NH2, -F 等) は近接
  ring segment に吸収する (`absorb_single_substituent=True`、
  `kind="ring_with_substituent"`)
"""
from __future__ import annotations

import logging
from typing import Any, List, Set, Tuple

from .models import Segment

logger = logging.getLogger(__name__)


def detect_ring_segments(
    mol: Any,
    allow_atom_sharing: bool = True,
    absorb_single_substituent: bool = True,
    start_segment_id: int = 0,
) -> Tuple[List[Segment], Set[int]]:
    """RDKit Mol から ring を抽出して Segment list として返す。

    Returns
    -------
    segments : List[Segment]
        ring segment のリスト (kind="ring" or "ring_with_substituent")。
    used_atoms : Set[int]
        ring segment に含まれた atom idx set。chain_splitter で除外する用。
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples (each ring is a tuple of atom idx)

    segments: List[Segment] = []
    used_atoms: Set[int] = set()
    sid = start_segment_id

    for ring_idx, ring_atoms in enumerate(atom_rings):
        atom_set = set(ring_atoms)
        if not allow_atom_sharing:
            atom_set = atom_set - used_atoms
            if not atom_set:
                continue
        seg = Segment(
            segment_id=sid,
            atom_indices=sorted(atom_set),
            kind="ring",
            ring_indices=[ring_idx],
        )
        segments.append(seg)
        used_atoms.update(atom_set)
        sid += 1

    if absorb_single_substituent:
        for seg in segments:
            _absorb_substituents(mol, seg, used_atoms)

    logger.info(
        "detect_ring_segments: %d ring segment(s), %d ring atoms (incl. absorbed substituents)",
        len(segments), len(used_atoms),
    )
    return segments, used_atoms


def _absorb_substituents(mol: Any, seg: Segment, used_atoms: Set[int]) -> None:
    """ring atom に attach する 1 heavy atom 置換基を seg に吸収する。

    Rule: ring atom の heavy neighbor で、ring atom 以外に heavy neighbor を
    持たない atom (= 末端 heavy atom = -OH の O、-NH2 の N、-F の F 等) を
    seg.atom_indices に追加。
    """
    seg_set = set(seg.atom_indices)
    to_add: List[int] = []
    for ring_atom in seg.atom_indices:
        atom = mol.GetAtomWithIdx(ring_atom)
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 1:
                continue
            nb_idx = nb.GetIdx()
            if nb_idx in seg_set or nb_idx in used_atoms:
                continue
            # nb の他の heavy neighbor が ring_atom のみなら末端置換基
            other_heavy = [
                n2.GetIdx() for n2 in nb.GetNeighbors()
                if n2.GetAtomicNum() != 1 and n2.GetIdx() != ring_atom
            ]
            if not other_heavy:
                to_add.append(nb_idx)

    if to_add:
        seg.atom_indices = sorted(seg_set | set(to_add))
        used_atoms.update(to_add)
        seg.kind = "ring_with_substituent"
