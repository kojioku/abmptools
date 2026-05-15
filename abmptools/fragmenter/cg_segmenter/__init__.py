# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter
---------------------------------
粗視化 (CG) セグメント分割ツール。

`abmptools.fragmenter` (FMO 計算用) と姉妹関係:

| 観点 | `fragmenter` (FMO) | `cg_segmenter` (CG) |
|---|---|---|
| BDA / BAA | あり | なし |
| 切断方法 | 擬似 (segment_data.dat に bond マーク) | **物理的に mol を分割、cap 付与** |
| 出力 | `segment_data.dat` (log2config 互換) | per-segment **PDB + XYZ + JSON** |
| atom 共有 | なし | **あり** (fused ring の共有 atom) |
| ring 内切断 | 除外 (環は 1 fragment) | **環ごとに別 segment** |

両者で共通利用するもの:
- `abmptools.fragmenter.pdb_loader.load_pdb_molecules` (連結成分単位)
- `abmptools.fragmenter.auto_split.suggest_cuts` (chain 部分の target_mw 切断)

主要 API:
    CGSegmenterConfig    全体設定 (JSON 往復可能)
    Segment              1 segment の情報 (atom_indices / cap_atoms / smiles)
    CapAtom              切断面に追加する cap (H or CH3)
    CGSegmenter          ロード→ring 検出→chain 切断→cap→export の orchestrator

CLI:
    python -m abmptools.fragmenter.cg_segmenter build --pdb input.pdb \\
        --output-dir ./segs --target-mw 200

依存:
    - rdkit-pypi (`[fragmenter]` extras と共通)
"""
from .models import (
    CGSegmenterConfig,
    CapAtom,
    Segment,
    SegmentResult,
)

try:
    from .ring_detector import detect_ring_segments
    from .chain_splitter import split_chain
    from .cap_attach import attach_caps
    from .exporter import export_segments, render_segments_svg
    from .dpdgen_exporter import export_dpdgen
    from .orchestrator import CGSegmenter
    _CORE_AVAILABLE = True
except ImportError:
    _CORE_AVAILABLE = False

# Optional: ipywidgets-dependent UI (Jupyter extras only)
try:
    from .notebook_ui import open_panel  # noqa: F401
    _UI_AVAILABLE = True
except ImportError:
    _UI_AVAILABLE = False

__all__ = [
    "CGSegmenterConfig",
    "CapAtom",
    "Segment",
    "SegmentResult",
]
if _CORE_AVAILABLE:
    __all__ += [
        "detect_ring_segments",
        "split_chain",
        "attach_caps",
        "export_segments",
        "render_segments_svg",
        "export_dpdgen",
        "CGSegmenter",
    ]
if _UI_AVAILABLE and _CORE_AVAILABLE:
    __all__.append("open_panel")
