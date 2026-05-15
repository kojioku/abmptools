# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.cap_attach
--------------------------------------------
各 Segment の切断面に cap atom (H or CH3) を 3D 配置で追加する。

Rule:
    - 切断面の端 atom (= segment 内の atom) の元素が:
      - C / halogen / その他: **H cap** (1 atom 追加)
      - N / O / S / P (`hetero_cap_methyl_elements`): **CH3 cap**
        (central C 1 atom + 周辺 3 H = 計 4 atoms 追加; CapAtom には central C
        のみ記録、3 H は exporter が展開)

位置計算:
    - 端 atom から 切断相手 neighbor 方向の単位ベクトル × 結合長
    - 結合長は `BOND_LEN` table (C-H 1.09, C-C 1.54, etc.)
"""
from __future__ import annotations

import logging
import math
from typing import Any, List, Optional, Sequence, Tuple

from .models import CapAtom, Segment

logger = logging.getLogger(__name__)

# Single bond length defaults (Å). 詳細値は CRC handbook 等の代表値。
BOND_LEN = {
    ("C", "H"): 1.09,
    ("N", "H"): 1.01,
    ("O", "H"): 0.97,
    ("S", "H"): 1.34,
    ("P", "H"): 1.42,
    ("F", "H"): 0.92,   # 通常 F は H と直接 bond しないがフォールバック
    ("Cl", "H"): 1.27,
    ("Br", "H"): 1.41,
    ("I", "H"): 1.60,
    ("C", "C"): 1.54,
    ("N", "C"): 1.47,
    ("O", "C"): 1.43,
    ("S", "C"): 1.82,
    ("P", "C"): 1.84,
}


def attach_caps(
    mol: Any,
    segments: List[Segment],
    hetero_cap_methyl_elements: Optional[Sequence[str]] = None,
) -> None:
    """各 Segment の境界 atom に cap atom を追加する (in-place 更新)。

    Parameters
    ----------
    mol
        RDKit Mol (3D conformer を持っていること)。
    segments
        Segment list (`atom_indices` が埋まっている)。
    hetero_cap_methyl_elements
        切断面の端 atom がこの元素なら CH3 cap、それ以外は H cap。
        Default ["N", "O", "S", "P"] (不要 H 結合スポット発生を避ける)。
    """
    if hetero_cap_methyl_elements is None:
        hetero_cap_methyl_elements = ("N", "O", "S", "P")
    hetero_set = set(hetero_cap_methyl_elements)

    conf = mol.GetConformer()

    for seg in segments:
        seg_set = set(seg.atom_indices)
        seg.cap_atoms = []  # 上書き (前回計算済の cap をクリア)
        for a in seg.atom_indices:
            atom = mol.GetAtomWithIdx(a)
            end_symbol = atom.GetSymbol()
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 1:
                    continue
                nb_idx = nb.GetIdx()
                # nb が同 segment にもあれば cap 不要 (atom 共有 OK の場合)
                if nb_idx in seg_set:
                    continue

                # cap 元素を決定
                if end_symbol in hetero_set:
                    cap_element = "C"      # CH3 cap (central C)
                    is_methyl = True
                else:
                    cap_element = "H"
                    is_methyl = False

                # 3D 位置: end atom → nb 方向の単位ベクトル × 結合長
                pos_end = conf.GetAtomPosition(a)
                pos_nb = conf.GetAtomPosition(nb_idx)
                dx = pos_nb.x - pos_end.x
                dy = pos_nb.y - pos_end.y
                dz = pos_nb.z - pos_end.z
                norm = math.sqrt(dx * dx + dy * dy + dz * dz)
                if norm < 1e-6:
                    direction = (1.0, 0.0, 0.0)
                else:
                    direction = (dx / norm, dy / norm, dz / norm)

                bond_len = BOND_LEN.get((end_symbol, cap_element), 1.50)
                cap_pos = (
                    pos_end.x + bond_len * direction[0],
                    pos_end.y + bond_len * direction[1],
                    pos_end.z + bond_len * direction[2],
                )

                seg.cap_atoms.append(CapAtom(
                    parent_atom_idx=a,
                    element=cap_element,
                    position=cap_pos,
                    is_methyl_cap=is_methyl,
                ))

    n_caps = sum(len(s.cap_atoms) for s in segments)
    logger.info("attach_caps: total %d cap site(s) across %d segment(s)", n_caps, len(segments))


def methyl_hydrogen_positions(
    central_c_pos: Tuple[float, float, float],
    parent_atom_pos: Tuple[float, float, float],
    bond_len_ch: float = 1.09,
) -> List[Tuple[float, float, float]]:
    """CH3 cap の 3 H 原子の 3D 座標を返す (tetrahedral 配置の簡易近似)。

    parent_atom -- C(central) -- H × 3 の tetrahedral。central C から見て
    parent 方向の逆向き = methyl axis、それに垂直な平面上に 120 度ずつ H を配置。
    """
    cx, cy, cz = central_c_pos
    px, py, pz = parent_atom_pos
    # methyl axis: central C から parent と反対方向の単位ベクトル
    ax = cx - px
    ay = cy - py
    az = cz - pz
    norm = math.sqrt(ax * ax + ay * ay + az * az)
    if norm < 1e-6:
        ax, ay, az = 0.0, 0.0, 1.0
    else:
        ax, ay, az = ax / norm, ay / norm, az / norm

    # axis に垂直な任意の単位ベクトル e1 を作る
    if abs(ax) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)
    # e1 = axis × ref / |...|
    e1x = ay * ref[2] - az * ref[1]
    e1y = az * ref[0] - ax * ref[2]
    e1z = ax * ref[1] - ay * ref[0]
    nrm = math.sqrt(e1x * e1x + e1y * e1y + e1z * e1z)
    e1x, e1y, e1z = e1x / nrm, e1y / nrm, e1z / nrm
    # e2 = axis × e1
    e2x = ay * e1z - az * e1y
    e2y = az * e1x - ax * e1z
    e2z = ax * e1y - ay * e1x

    # H positions: tetrahedral angle = 109.47°, axis 方向に cos(70.53°)、
    # 垂直平面上に sin(70.53°)、120° ずつ回転
    cos_t = math.cos(math.radians(70.53))  # 0.3338
    sin_t = math.sin(math.radians(70.53))  # 0.9428
    positions = []
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        rx = sin_t * (math.cos(phi) * e1x + math.sin(phi) * e2x) + cos_t * ax
        ry = sin_t * (math.cos(phi) * e1y + math.sin(phi) * e2y) + cos_t * ay
        rz = sin_t * (math.cos(phi) * e1z + math.sin(phi) * e2z) + cos_t * az
        positions.append((
            cx + bond_len_ch * rx,
            cy + bond_len_ch * ry,
            cz + bond_len_ch * rz,
        ))
    return positions
