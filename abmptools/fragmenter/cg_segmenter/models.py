# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.models
----------------------------------------
データクラス定義 (rdkit 非依存、JSON 往復可能)。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class CapAtom:
    """切断面に追加する cap 原子の情報。

    Attributes
    ----------
    parent_atom_idx
        cap を付ける元 atom の Mol 内 idx (segment 内の atom)。
    element
        "H" (default cap) or "C" (CH3 cap の central C; 周辺 3 H は別途生成)。
    position
        3D 座標 (Å)。元 bond 方向ベクトルから配置 (C-H=1.09 Å, C-C=1.54 Å)。
    is_methyl_cap
        True なら central C + 3 H で計 4 atoms 追加 (-CH3 group)。
        False なら H 1 atom のみ。
    """

    parent_atom_idx: int
    element: str
    position: Tuple[float, float, float]
    is_methyl_cap: bool = False

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["position"] = list(d["position"])
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CapAtom":
        pos = tuple(d.get("position", (0.0, 0.0, 0.0)))
        return cls(
            parent_atom_idx=d["parent_atom_idx"],
            element=d["element"],
            position=pos,
            is_methyl_cap=d.get("is_methyl_cap", False),
        )


@dataclass
class Segment:
    """1 segment の情報 (CG 粒子に対応する原子グループ)。

    Attributes
    ----------
    segment_id
        系内での通し番号 (0-origin)。
    atom_indices
        元 Mol 内 atom idx の list (heavy atom)。fused ring の共有 atom は
        複数 Segment の `atom_indices` に同時に含まれる場合がある。
    cap_atoms
        切断面に追加された cap atom list。
    smiles
        cap 込み canonical SMILES (検証・表示用)。空文字なら計算未実行。
    kind
        "ring": 環構造、"chain": chain 部分、"ring_with_substituent": 環 +
        小さい置換基 (例: 環の -OH を吸収)。
    ring_indices
        この segment が含む RDKit ring index list (`mol.GetRingInfo()` の出力
        順序)。`kind="chain"` なら空。
    """

    segment_id: int
    atom_indices: List[int] = field(default_factory=list)
    cap_atoms: List[CapAtom] = field(default_factory=list)
    smiles: str = ""
    kind: str = "chain"
    ring_indices: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "segment_id": self.segment_id,
            "atom_indices": list(self.atom_indices),
            "cap_atoms": [c.to_dict() for c in self.cap_atoms],
            "smiles": self.smiles,
            "kind": self.kind,
            "ring_indices": list(self.ring_indices),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Segment":
        return cls(
            segment_id=d["segment_id"],
            atom_indices=list(d.get("atom_indices", [])),
            cap_atoms=[CapAtom.from_dict(c) for c in d.get("cap_atoms", [])],
            smiles=d.get("smiles", ""),
            kind=d.get("kind", "chain"),
            ring_indices=list(d.get("ring_indices", [])),
        )

    @property
    def n_atoms(self) -> int:
        """この segment の原子数 (cap atom 含む。methyl cap は central C + 3 H = 4)。"""
        n = len(self.atom_indices)
        for c in self.cap_atoms:
            n += 4 if c.is_methyl_cap else 1
        return n


@dataclass
class CGSegmenterConfig:
    """CG segmenter の入力設定。

    Attributes
    ----------
    pdb_path
        入力 PDB ファイル。
    output_dir
        出力ディレクトリ (per-segment PDB + XYZ + summary JSON)。
    target_mw
        chain 部分の C-C 切断 (= fragmenter.auto_split) の目安分子量 (g/mol)。
    separate_rings
        True なら各 ring を別 segment に。False なら ring を chain と同様に
        target_mw 基準で切る (実験的)。
    allow_atom_sharing
        True なら fused ring の共有 atom を両 ring segment に含める。False なら
        共有 atom を片方に寄せる (どちらに寄せるかは決定論的に若い atom_idx 側)。
    hetero_cap_methyl_elements
        切断面の atom がこのリストの元素なら CH3 cap、それ以外なら H cap。
        Default: N / O / S / P (不要な水素結合スポット発生を避ける用)。
    absorb_single_substituent
        True なら ring に直接 attach する 1 heavy atom 置換基 (例: -OH, -NH2,
        -F) を近接 ring segment に吸収する (`kind="ring_with_substituent"`)。
        False なら別 chain segment 化。
    """

    pdb_path: str
    output_dir: str = "./cg_segments"
    target_mw: float = 200.0
    separate_rings: bool = True
    allow_atom_sharing: bool = True
    hetero_cap_methyl_elements: List[str] = field(
        default_factory=lambda: ["N", "O", "S", "P"]
    )
    absorb_single_substituent: bool = True

    def to_json(self, path: Optional[str] = None) -> str:
        text = json.dumps(asdict(self), indent=2, ensure_ascii=False)
        if path:
            Path(path).write_text(text)
        return text

    @classmethod
    def from_json(cls, path: str) -> "CGSegmenterConfig":
        data = json.loads(Path(path).read_text())
        return cls(**data)


@dataclass
class SegmentResult:
    """最終分割結果。

    Attributes
    ----------
    segments
        全 Segment。
    total_atoms_with_cap
        cap 込みの総 atom 数 (= sum(s.n_atoms))。
    shared_atom_pairs
        atom が複数 segment に属するペア list。各要素は (atom_idx,
        segment_id_a, segment_id_b)。検証 / レポート用。
    """

    segments: List[Segment]
    total_atoms_with_cap: int = 0
    shared_atom_pairs: List[Tuple[int, int, int]] = field(default_factory=list)
