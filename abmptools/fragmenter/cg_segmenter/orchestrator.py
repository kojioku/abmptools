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
