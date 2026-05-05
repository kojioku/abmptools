# abmptools.fragmenter — FMO 自動フラグメント分割ツール

PDB から小分子・脂質・ポリマーのフラグメント分割を提案し、
ユーザーが Jupyter UI または CLI ヘッドレス経路で確定後、
`abmptools.log2config` 互換の `segment_data.dat` を出力するサブパッケージ。

タンパク質・DNA は対象外（既存 `abmptools.log2config` 経路で扱う）。

## 設計の要点

- **同一分子判定**: RDKit canonical SMILES。系中に同じ分子が複数あれば
  代表 1 分子だけ表示・編集し、確定時に全コピーへ展開する
- **C-C 切断ルール**: 環内 C-C / 多重結合 / ヘテロ原子隣接を除外、
  目安分子量（default 200 g/mol）で区切る
- **2 系統 UI**:
  - **A**: Jupyter + RDKit + ipywidgets（インタラクティブ）
  - **C**: CLI で SVG + JSON 出力 → 編集 → 再適用（ヘッドレス）
- **ポリマー対応**: 鎖長違いを γ 経路（GUI で「同一パターン」を明示指定）で同一視

## クイックスタート

```bash
# 自動提案 + レビューバンドル出力
python -m abmptools.fragmenter suggest --pdb input.pdb --target-mw 200 --output-dir ./review

# (review/ 内の SVG / JSON を確認・手動編集)

# 確定後、segment_data.dat 出力
python -m abmptools.fragmenter apply --pdb input.pdb --review-dir ./review --output segment_data.dat
```

詳細は [`docs/fragmenter.md`](../../docs/fragmenter.md) を参照（P6 で作成予定）。

## インストール

```bash
pip install -e '.[fragmenter]'           # コア (rdkit-pypi)
pip install -e '.[fragmenter,jupyter]'   # Notebook UI 含む (ipywidgets, ipykernel)
```

## モジュール構成

| ファイル | 役割 | 実装 phase |
|---|---|---|
| `models.py` | データクラス (FragmenterConfig / CutSite / MoleculeGroup / FragmentResult) | P1 |
| `pdb_loader.py` | PDB → RDKit Mol、連結成分分解、bond perception フォールバック | P1 |
| `grouping.py` | canonical SMILES グループ化 | P1 |
| `auto_split.py` | C-C 切断候補の自動提案 (環内/多重/ヘテロ隣接除外、target_mw 基準) | P2 |
| `cut_apply.py` | CutSite list → 各分子内のフラグメント割当 | P2 |
| `expand_to_system.py` | 1 分子の分割を全コピーへ展開、segment_data.dat 出力 | P2 |
| `headless_io.py` | C 経路: SVG + JSON 出力 / 編集 JSON 読込 | P3 |
| `__main__.py` | CLI (suggest / apply / example) | P3 |
| `notebook_ui.py` | A 経路: Jupyter + RDKit + ipywidgets UI | P4 |
| `polymer.py` | γ 経路: 同一視グループのユーザー指定 | P5 |

## 関連サブパッケージ

- `abmptools.log2config` — タンパク質/DNA の残基ごとフラグメント (本ツールの前段／代替)
- `abmptools.pdb2fmo` — `segment_data.dat` + PDB → ABINIT-MP `ajf` 入力 (本ツールの後段)
