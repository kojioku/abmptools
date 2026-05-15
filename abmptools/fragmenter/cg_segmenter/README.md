# abmptools.fragmenter.cg_segmenter — CG セグメント構築 (モジュール構成)

`abmptools.fragmenter` (FMO 用) の姉妹サブモジュール。**粗視化 (CG) 粒子に対応する
segment** を PDB から自動構築する。

詳細は [`docs/cg_segmenter.md`](../../../docs/cg_segmenter.md) を参照。

## モジュール構成

| ファイル | 役割 |
|---|---|
| `models.py` | データクラス (`CGSegmenterConfig` / `Segment` / `CapAtom` / `SegmentResult`) |
| `ring_detector.py` | RDKit `GetRingInfo` → ring segment 抽出 (fused = atom 共有) |
| `chain_splitter.py` | ring 外 chain の target_mw walk 切断 (`auto_split._atom_total_mw` 流用) |
| `cap_attach.py` | 切断面に H or CH3 cap を 3D 配置 |
| `exporter.py` | per-segment PDB + XYZ + summary JSON 出力 |
| `orchestrator.py` | `CGSegmenter` クラス (state 管理 + pipeline) |
| `__main__.py` | CLI (`python -m abmptools.fragmenter.cg_segmenter`) |

## クイックスタート

```bash
# CLI
python -m abmptools.fragmenter.cg_segmenter build \
    --pdb input.pdb --output-dir ./segs --target-mw 200

# Python API
python -c "
from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
config = CGSegmenterConfig(pdb_path='input.pdb', target_mw=200)
sg = CGSegmenter.from_pdb(config)
print(f'{len(sg.segments)} segments')
sg.export()
"
```
