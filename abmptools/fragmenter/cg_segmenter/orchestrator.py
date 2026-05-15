# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.orchestrator
----------------------------------------------
CGSegmenter — ロード → ring 検出 → chain 切断 → cap 付与 → 出力 を統合する
state 管理クラス。
"""
from __future__ import annotations

import logging
from typing import Any, List, Optional

from ..pdb_loader import load_pdb_molecules
from .cap_attach import attach_caps
from .chain_splitter import split_chain
from .exporter import export_segments
from .models import CGSegmenterConfig, Segment, SegmentResult
from .ring_detector import detect_ring_segments

logger = logging.getLogger(__name__)


class CGSegmenter:
    """CG segmentation の state 管理 orchestrator。

    使用例:
        >>> config = CGSegmenterConfig(pdb_path="cholesterol.pdb", target_mw=200)
        >>> sg = CGSegmenter.from_pdb(config)
        >>> result = sg.export()
        >>> print(f"{len(result.segments)} segments to {config.output_dir}")
    """

    def __init__(self, config: CGSegmenterConfig):
        self.config = config
        self.mol: Optional[Any] = None
        self.segments: List[Segment] = []
        self.result: Optional[SegmentResult] = None

    @classmethod
    def from_pdb(cls, config: CGSegmenterConfig) -> "CGSegmenter":
        """PDB load + segmentation まで一括実行。

        単一分子 PDB を想定 (現状)。複数分子があれば最初の 1 つを使う。
        """
        sg = cls(config)
        loaded = load_pdb_molecules(config.pdb_path)
        if len(loaded) == 0:
            raise RuntimeError(f"No molecule loaded from {config.pdb_path}")
        if len(loaded) > 1:
            logger.warning(
                "cg_segmenter currently uses only the first connected component "
                "(got %d). Multi-molecule support is future work.",
                len(loaded),
            )
        sg.mol = loaded[0].mol
        sg.segment()
        return sg

    def segment(self) -> None:
        """ring 検出 + chain 切断 + cap 付与 (in-place で self.segments を更新)。"""
        if self.mol is None:
            raise RuntimeError("CGSegmenter.mol is None. Call from_pdb() first.")

        ring_segs: List[Segment] = []
        used_in_rings: set = set()
        if self.config.separate_rings:
            ring_segs, used_in_rings = detect_ring_segments(
                self.mol,
                allow_atom_sharing=self.config.allow_atom_sharing,
                absorb_single_substituent=self.config.absorb_single_substituent,
            )

        chain_segs = split_chain(
            self.mol,
            used_in_rings,
            target_mw=self.config.target_mw,
            start_segment_id=len(ring_segs),
        )

        all_segs = ring_segs + chain_segs
        attach_caps(
            self.mol, all_segs,
            hetero_cap_methyl_elements=self.config.hetero_cap_methyl_elements,
        )
        self.segments = all_segs

        logger.info(
            "CGSegmenter.segment: %d ring + %d chain = %d total segment(s)",
            len(ring_segs), len(chain_segs), len(all_segs),
        )

    def export(self, output_dir: Optional[str] = None) -> SegmentResult:
        """PDB + XYZ + summary JSON 出力。

        Parameters
        ----------
        output_dir
            None なら ``self.config.output_dir`` を使う。
        """
        if self.mol is None:
            raise RuntimeError("CGSegmenter.mol is None.")
        out_dir = output_dir or self.config.output_dir
        self.result = export_segments(self.mol, self.segments, out_dir)
        return self.result

    def export_dpdgen(
        self,
        output_dir: Optional[str] = None,
        monomer_name: str = "cg",
        box_size: tuple = (10.0, 10.0, 10.0),
        total_num: int = 100000,
        step: int = 100,
        aij_file: str = "aij.dat",
    ) -> tuple:
        """DPDgen 入力 (monomer + calc_sett) を生成する。

        Returns
        -------
        (monomer_path, calc_sett_path) -- both Path objects
        """
        from .dpdgen_exporter import export_dpdgen
        if self.mol is None:
            raise RuntimeError("CGSegmenter.mol is None.")
        out_dir = output_dir or self.config.output_dir
        return export_dpdgen(
            self.segments, self.mol, out_dir,
            monomer_name=monomer_name,
            box_size=box_size,
            total_num=total_num,
            step=step,
            aij_file=aij_file,
        )

    # ------------------------------------------------------------------
    # Manual edit operations (Jupyter UI から呼ばれる)
    # ------------------------------------------------------------------

    def move_atom(self, atom_idx: int, target_seg_id: int, shared: bool = False) -> None:
        """atom を target segment に移動する (in-place、cap 再計算)。

        Parameters
        ----------
        atom_idx
            移動する atom の Mol 内 idx (heavy)。
        target_seg_id
            移動先 Segment.segment_id。target に既に atom_idx があれば no-op。
        shared
            False (default): 元 segment(s) から atom_idx を削除して exclusive に。
            True: 元 segment(s) にも残して両方に属する shared atom 化。
        """
        target = next((s for s in self.segments if s.segment_id == target_seg_id), None)
        if target is None:
            raise ValueError(f"Target segment {target_seg_id} not found")

        sources = [
            s for s in self.segments
            if atom_idx in s.atom_indices and s.segment_id != target_seg_id
        ]

        if atom_idx not in target.atom_indices:
            target.atom_indices = sorted(set(target.atom_indices) | {atom_idx})

        if not shared:
            for s in sources:
                s.atom_indices = [a for a in s.atom_indices if a != atom_idx]

        self._recompute_caps()
        logger.info(
            "move_atom: atom %d -> seg %d (shared=%s, removed from %d source seg(s))",
            atom_idx, target_seg_id, shared, 0 if shared else len(sources),
        )

    def toggle_cap(self, segment_id: int, cap_index: int) -> None:
        """指定 Segment の cap_index 番目の cap の element / is_methyl_cap を toggle する。

        H ↔ CH3 を切替。位置は parent atom 方向の単位ベクトル × 新結合長で
        再計算する。
        """
        import math
        from .cap_attach import BOND_LEN
        seg = next((s for s in self.segments if s.segment_id == segment_id), None)
        if seg is None:
            raise ValueError(f"Segment {segment_id} not found")
        if cap_index < 0 or cap_index >= len(seg.cap_atoms):
            raise ValueError(
                f"cap_index {cap_index} out of range (0..{len(seg.cap_atoms) - 1})"
            )
        cap = seg.cap_atoms[cap_index]
        cap.is_methyl_cap = not cap.is_methyl_cap
        cap.element = "C" if cap.is_methyl_cap else "H"

        parent = self.mol.GetAtomWithIdx(cap.parent_atom_idx)
        conf = self.mol.GetConformer()
        ppos = conf.GetAtomPosition(cap.parent_atom_idx)
        dx = cap.position[0] - ppos.x
        dy = cap.position[1] - ppos.y
        dz = cap.position[2] - ppos.z
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        if norm < 1e-6:
            return
        new_len = BOND_LEN.get((parent.GetSymbol(), cap.element), 1.50)
        cap.position = (
            ppos.x + new_len * dx / norm,
            ppos.y + new_len * dy / norm,
            ppos.z + new_len * dz / norm,
        )
        logger.info(
            "toggle_cap: seg %d cap %d -> %s",
            segment_id, cap_index, "CH3" if cap.is_methyl_cap else "H",
        )

    def delete_segment(self, segment_id: int) -> None:
        """指定 Segment を削除 (atoms はどの seg にも属さない状態に)。cap 再計算。"""
        before = len(self.segments)
        self.segments = [s for s in self.segments if s.segment_id != segment_id]
        if len(self.segments) == before:
            raise ValueError(f"Segment {segment_id} not found")
        self._recompute_caps()
        logger.info("delete_segment: removed seg %d, %d segments remain",
                    segment_id, len(self.segments))

    def re_segment(
        self,
        target_mw: Optional[float] = None,
        separate_rings: Optional[bool] = None,
        allow_atom_sharing: Optional[bool] = None,
        absorb_single_substituent: Optional[bool] = None,
    ) -> None:
        """config を更新して self.segments を再計算する (全 segment を上書き)。"""
        if target_mw is not None:
            self.config.target_mw = target_mw
        if separate_rings is not None:
            self.config.separate_rings = separate_rings
        if allow_atom_sharing is not None:
            self.config.allow_atom_sharing = allow_atom_sharing
        if absorb_single_substituent is not None:
            self.config.absorb_single_substituent = absorb_single_substituent
        self.segment()
        logger.info(
            "re_segment: %d segment(s) (target_mw=%.1f)",
            len(self.segments), self.config.target_mw,
        )

    def _recompute_caps(self) -> None:
        """全 segment の cap_atoms を attach_caps で再計算 (in-place)。"""
        if self.mol is None:
            return
        attach_caps(
            self.mol, self.segments,
            hetero_cap_methyl_elements=self.config.hetero_cap_methyl_elements,
        )
