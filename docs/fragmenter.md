# abmptools.fragmenter — FMO 自動フラグメント分割ツール

## 概要

`abmptools.fragmenter` は、PDB から **FMO 計算用のフラグメント分割を自動提案**し、
ユーザーがインタラクティブまたはヘッドレスで手動修正後、`abmptools.log2config`
互換の `segment_data.dat` を出力するサブパッケージ (v1.21.0+)。

| 担当範囲 | 経路 |
|---|---|
| **タンパク質 / DNA** | 対象外 — 既存の `abmptools.log2config` 経路で扱う (ABINIT-MP の `nprint=0` log → `segment_data.dat`) |
| **小分子 / 脂質 / ポリマー** | 本サブパッケージ — RDKit ベースの C-C 切断アルゴリズムで自動提案 |

設計の要点:

- **同一分子判定**: RDKit canonical SMILES を group key とし、複数コピーは
  代表 1 分子だけ表示・編集して全コピーへ展開する
- **C-C 切断ルール**: 環内 C-C / 多重結合 / ヘテロ原子隣接を除外、目安分子量
  (default 200 g/mol) で区切る
- **2 系統 UI**:
  - **A**: Jupyter + RDKit + ipywidgets (`open_panel`)
  - **C**: CLI で SVG + JSON 出力 → 編集 → 再 apply (`python -m abmptools.fragmenter`)
- **ポリマー対応**: 鎖長違い (PE N=10 vs N=11 等) は γ 経路 (`declare_same_pattern`)
  で同一視を明示指定する

## インストール

```bash
# コア (RDKit のみ)
pip install -e '.[fragmenter]'

# Notebook UI も使う場合
pip install -e '.[fragmenter,jupyter]'
```

`extras_require`:

| extras | 追加依存 | 用途 |
|---|---|---|
| `[fragmenter]` | `rdkit-pypi >= 2022.09` | コア (PDB ロード / SMILES / SVG / 切断) |
| `[jupyter]` | `ipywidgets >= 8.0`, `ipykernel` | A 経路 (Notebook UI) |

省略可能な外部ツール:

- **Open Babel CLI (`obabel`)**: PDB の bond perception フォールバック (subprocess)。
  AmberTools / amorphous extras と同梱されることが多い。

## クイックスタート

### CLI ヘッドレス経路 (C)

```bash
# 1) 自動提案 + レビューバンドル出力
python -m abmptools.fragmenter suggest \
    --pdb input.pdb --target-mw 200 --output-dir ./review

# 2) 出力された ./review/group_*.json の cut_sites.enabled を手動で toggle
#    ./review/group_*.svg で構造と切断点を可視化

# 3) 編集後 JSON を読み込んで segment_data.dat 出力
python -m abmptools.fragmenter apply \
    --pdb input.pdb --review-dir ./review --output segment_data.dat

# 4) pdb2fmo で ajf に展開
python -m abmptools.pdb2fmo -i input.pdb -p input_param
```

### Jupyter UI (A)

```python
from abmptools.fragmenter import FragmenterConfig, AutoFragmenter, open_panel

config = FragmenterConfig(
    pdb_path="input.pdb",
    target_mw=200.0,
    output_dir="./fragmenter_out",
)
af = AutoFragmenter.from_pdb(config)
af.suggest_cuts()
open_panel(af)   # ipywidgets パネルが notebook 上に表示される
```

操作:
- Dropdown で分子グループを選択
- SVG (heavy atom only) で代表分子の構造、提案 cut sites が赤色ハイライト
- 各 cut site の Checkbox で `enabled` を toggle (re-render)
- "Export segment_data.dat" ボタンで全コピーへ展開して出力

## アルゴリズム詳細

### 1. PDB ロードと連結成分分解 (`pdb_loader.py`)

`Chem.MolFromPDBFile` を 4 段階フォールバックで呼ぶ:

| Strategy | 引数 | 用途 |
|---|---|---|
| 1 | `proximityBonding=True, sanitize=True` | 標準 (CONECT 優先 + 近接補完) |
| 2 | `proximityBonding=False, sanitize=True` | 隣接分子間誤結合の回避 (CONECT のみ信用) |
| 3 | `sanitize=False` | best effort (valence エラー時) |
| 4 | `obabel` 前処理 → 再ロード | bond perception フォールバック |

ロード成功後、`Chem.GetMolFrags` で連結成分に分解し、各 component を 1 つの
`LoadedMolecule` として返す (`mol`, `residue_name`, `atom_indices_in_pdb`)。

### 2. Canonical SMILES グループ化 (`grouping.py`)

各 `LoadedMolecule` を heavy atom only (`Chem.RemoveHs`) にしてから
`Chem.MolToSmiles(canonical=True)` で group key を計算。`MoleculeGroup`
にまとめ、`n_copies`, `member_mol_indices`, `residue_names` を集計する。

`include_residue_name=True` を指定すると `(SMILES, residue_name)` のタプルを
key にして、同 SMILES でも残基名が異なれば別グループに分ける。

### 3. C-C 切断候補の自動提案 (`auto_split.py`)

```
suggest_cuts(mol, config)
├── _enumerate_candidate_bonds   フィルタ: SINGLE / 両端 C / 環内除外 /
│                                 多重結合除外 / ヘテロ隣接除外
├── _find_main_chain_heavy       graph diameter (heavy atom only 2-pass BFS)
└── _select_cuts_along_path      主鎖を walk、累積 MW + 側鎖 MW を加算、
                                  累積 ≥ target_mw のたびに candidate bond で切断
```

**フィルタ詳細** (デフォルト全て ON):

| フィルタ | 説明 | 例 |
|---|---|---|
| 環内 C-C 除外 | `bond.IsInRing()` 真なら除外 | ベンゼン / シクロヘキサン環内は切らない |
| 多重結合除外 | `BondType != SINGLE` を除外 | C=C / C≡C は切らない |
| ヘテロ隣接除外 | 両端 C の隣接にヘテロ原子があれば除外 | アミド C-N 隣接 / エステル C-O 隣接 / カルボニル C=O 隣接 |

**主鎖検出 (graph diameter)**:

1. 任意の heavy atom から BFS で最遠 atom `far1` を見つける
2. `far1` から再度 BFS で最遠 `far2` を見つける
3. `far1` ↔ `far2` の shortest path が直径 (= 主鎖と定義)

**MW walk**:

```
for i in range(len(path) - 1):
    cum_mw += atom_mw(path[i])
    cum_mw += sidechain_mw(path[i], main_chain_set)
    bond = bond_between(path[i], path[i+1])
    if bond is in candidates and cum_mw >= target_mw:
        cut(bond); cum_mw = 0
```

主鎖外の側鎖 heavy atom は `path[i]` のローカル MW として加算する。

### 4. cut 適用と fragment 抽出 (`cut_apply.py`)

`apply_cuts(mol, cut_sites)` は `RWMol` で `enabled=True` の bond を削除し、
`Chem.GetMolFrags(asMols=False)` で連結成分を取得。各 fragment について:

- `atom_indices`: fragment 内 atom の mol 内 index (sorted)
- `charge`: formal charge sum
- `baa_atom_pairs`: BAA (Bond Attachment Atom) の `(this_atom, other_atom)` ペア list

### 5. 系全体への展開と segment_data 出力 (`expand_to_system.py`)

各 group の代表分子で得た `cut_sites` を全コピー (`member_mol_indices`) にも
適用し、log2config 互換の `segment_data` 構造を組み立てる:

```python
seg_data = [
    {
        "name": "<basename>_grp1_<short_smiles>",
        "atom":       [n_atoms_per_fragment, ...],     # n_copies × frag/copy 個
        "charge":     [charge_per_fragment, ...],
        "connect_num":[BAA_count_per_fragment, ...],
        "seg_info":   [[atom_idx_1origin, ...], ...],  # PDB 内 atom index (1-origin)
        "connect":    [[(this_atom, other_atom), ...], ...],
        "nummol_seg": [1],
        "repeat":     [1],
        "pair_file":  [],
        "multi_xyz":  "none",
    },
    # 次の group...
]
```

`write_segment_data` で log2config と同じ Python リテラル形式で書き出すため、
`pdb2fmo` の `ast.literal_eval` がそのまま読み込める。

### 6. ポリマー γ 経路 (`polymer.py`)

```python
from abmptools.fragmenter.polymer import declare_same_pattern

declare_same_pattern(groups, molecules, [
    ["CCCCCCCCCC", "CCCCCCCCCCC"],   # PE N=10 と N=11 を同一視
])
```

実装:
- master = 指定 SMILES のうち `len(cut_sites)` が最大の group
- 他の target group の主鎖を求め、master の cut が主鎖上のどの index 位置か調べる
- target 主鎖の同 index 位置の bond で切断する (target が短くて対応位置が
  なければスキップ)

## API リファレンス

### `FragmenterConfig`

```python
@dataclass
class FragmenterConfig:
    pdb_path: str
    output_dir: str = "./fragmenter_out"
    target_mw: float = 200.0
    exclude_ring_cc: bool = True
    exclude_multibond: bool = True
    exclude_heteroneighbor: bool = True
    skip_protein_dna: bool = True
    polymer_groups: List[List[str]] = []
    include_residue_name: bool = False

    def to_json(self, path: Optional[str] = None) -> str: ...
    @classmethod
    def from_json(cls, path: str) -> "FragmenterConfig": ...
```

### `AutoFragmenter`

```python
af = AutoFragmenter.from_pdb(config)        # ロード + グループ化
af.suggest_cuts()                            # 切断候補を自動提案
af.export_segment_data("segment_data.dat")   # 全コピーへ展開して出力
af.export_review_bundle("./review")          # SVG + JSON 出力
af.import_edited("./review")                 # 編集後 JSON を読み込み再描画
```

### 関数 API (FP スタイル)

```python
from abmptools.fragmenter import (
    load_pdb_molecules, group_by_smiles,
    suggest_cuts_for_groups, export_to_system,
    export_review_bundle, import_edited_review, resync_member_indices,
    declare_same_pattern,
)
```

## 既存ツールとの連携

| 段階 | ツール | 入力 | 出力 |
|---|---|---|---|
| 1 | `abmptools.fragmenter` | PDB | `segment_data.dat` (本ドキュメント) |
| 2 | `abmptools.pdb2fmo` | PDB + `segment_data.dat` + `input_param` | ABINIT-MP `ajf` |
| 3 | ABINIT-MP | `ajf` + PDB | LOG / CPF (FMO 計算結果) |
| 4 | `abmptools.getifiepieda` / `cpfmanager` | LOG / CPF | IFIE / PIEDA 解析 |

タンパク質 / DNA を含む系では、ABINIT-MP の自動分割を使った 1 度のテスト計算 →
`log2config` で `segment_data.dat` を抽出 → 本ツールで小分子部分のみ追加分割、
という併用ワークフローも可能。

## テスト

```bash
pytest tests/fragmenter/ -v
```

14 テスト:

| ファイル | 内容 | テスト数 |
|---|---|---|
| `test_basic.py` | dataclass roundtrip / load / group / suggest / export / CLI | 10 |
| `test_polymer.py` | γ 経路 (PE N=10 / N=11 同一視) | 4 |

PDB fixtures は `conftest.py` で RDKit が tmp_path_factory に毎回生成する
(propane, octane, propane×3 + acetone×2, PE N=10+N=11)。

## ライセンス

| 依存 | ライセンス | バンドル |
|---|---|---|
| abmptools 本体 | Apache-2.0 (v1.23.0+; ≤ v1.22.0 は MIT) | はい |
| RDKit | BSD-3-Clause | いいえ (PyPI 経由) |
| ipywidgets | BSD-3-Clause | いいえ (PyPI 経由) |
| Open Babel CLI | GPL-2.0 | いいえ (subprocess only、optional) |

## 関連リソース

- [`docs/overview.md`](overview.md) — abmptools 全サブパッケージ一覧
- [`docs/peptide_builders.md`](peptide_builders.md) — peptide ビルダー 3 系統
- `abmptools/log2config.py` — タンパク質 / DNA の残基ごとフラグメント
- `abmptools/pdb2fmo.py` — segment_data.dat → ajf
- `abmptools/fragmenter/README.md` — モジュール構成早見表

## バージョン履歴

| version | 日付 | 変更 |
|---|---|---|
| 1.21.0 | (本ドキュメント新設) | `abmptools.fragmenter` サブパッケージ初版 |
| 1.23.0 | 2026-05-10 | abmptools 本体 license を MIT → Apache-2.0 に移行 (機能変更なし、本ドキュメント中の license 表記を更新) |
