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
| `exporter.py` | per-segment PDB + XYZ + summary JSON 出力 + SVG ハイライト描画 |
| `dpdgen_exporter.py` | DPDgen 入力 (`{name}_monomer` + `{name}_calc_sett`) を出力 (path-based bond hierarchy + angle ポテンシャル) |
| `fcews_export.py` | FCEWS 入力 (`segment_data.dat` `mode='FMO'` + monomer `.xyz`) を出力。connect `[BDA,BAA]` (AJF 並び)、`connect_num`=BAA 数、atom 共有禁止 (partition 必須) |
| `notebook_ui.py` | Jupyter `open_panel()` (ipywidgets) — Move/Toggle/Delete/Re-segment/Export ボタン |
| `orchestrator.py` | `CGSegmenter` クラス (state 管理 + pipeline + edit ops `move_atom` / `toggle_cap` / `delete_segment` / `re_segment` / `export_dpdgen`) |
| `__main__.py` | CLI (`build` / `dpdgen` subcommand) |

## クイックスタート

```bash
# CLI -- segment 構築
python -m abmptools.fragmenter.cg_segmenter build \
    --pdb input.pdb --output-dir ./segs --target-mw 200

# CLI -- DPDgen 入力生成 (segment 構築 + monomer/calc_sett 出力)
python -m abmptools.fragmenter.cg_segmenter dpdgen \
    --pdb input.pdb --output-dir ./dpdgen_input \
    --monomer-name mymol --box 12

# CLI -- FCEWS 入力生成 (segment_data.dat + monomer xyz、mode='FMO')
python -m abmptools.fragmenter.cg_segmenter fcews \
    --pdb input.pdb --output-dir ./fcews_input \
    --target-mw 200 --name mymol

# Python API
python -c "
from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
# FCEWS 用は allow_atom_sharing=False (FMO フラグメントは atom 共有不可)
config = CGSegmenterConfig(pdb_path='input.pdb', target_mw=200, allow_atom_sharing=False)
sg = CGSegmenter.from_pdb(config)
print(f'{len(sg.segments)} segments')
sg.export()                                                       # PDB + XYZ + JSON
sg.export_dpdgen(output_dir='./dpdgen', monomer_name='mymol')    # monomer + calc_sett
sg.export_fcews(output_dir='./fcews_input', name='mymol')        # FCEWS segment_data.dat + xyz

# Jupyter UI
from abmptools.fragmenter.cg_segmenter import open_panel
open_panel(sg)                                                    # ipywidgets パネル
"
```
