# abmptools.fragmenter.cg_segmenter — CG セグメント構築ツール

## 概要

`abmptools.fragmenter.cg_segmenter` は、PDB から **粗視化 (Coarse-Grained / CG) 粒子に
対応するセグメント** を自動構築するサブモジュール (v1.24.0+)。

`abmptools.fragmenter` (FMO 用) の姉妹モジュールで、以下が異なる:

| 観点 | `fragmenter` (FMO) | **`cg_segmenter` (CG)** |
|---|---|---|
| 用途 | FMO 計算用フラグメント分割 | 粗視化粒子のセグメント分割 |
| BDA / BAA | あり (electron pair 帰属) | **なし** |
| 切断方法 | 擬似 (segment_data.dat に bond マーク) | **物理的に mol を分割し、cap atom (H or CH3) を付与** |
| 出力 | `segment_data.dat` (log2config 互換) | per-segment **PDB + XYZ + summary JSON** |
| atom 共有 | なし (各 atom は 1 fragment にのみ) | **あり** (fused ring の共有 atom は複数 segment に含まれる) |
| ring 内切断 | 除外 (環は 1 fragment) | **環ごとに別 segment** |

両者は **共通の `pdb_loader` / `auto_split` (target_mw walk)** を流用し、segment 抽出ロジックと出力部分のみ異なる。

## 公開版のスコープ (open-core)

本モジュール (公開版 abmptools) が担うのは **汎用 (reference) の CG segmentation と
入力生成**まで:

- ✅ ring 検出 / chain split / cap 付与 (汎用 segmentation)
- ✅ per-segment PDB / XYZ / JSON 出力
- ✅ **FCEWS segment_data.dat 生成** (FMO 入力。別系統)
- ✅ **DPD 入力 (`{name}_monomer` + `{name}_calc_sett`) 生成** → `abmptools.cg.dpd` で
  Cognac UDF / dpm に変換 (汎用既定パラメータ)

以下は **公開版の範囲外** (= 研究・production の当て込みノウハウ。別途 private 拡張で対応):

- 系ごとに最適化した bond/angle パラメータ・aij/chi セット
- 命名済み bead-typing / 系特化の分類 (例: `chobond_alkyl12` 等)
- 5 環以上の縮環 chain・anthracene 等への特化 default
- 実系の production workflow

公開版の DPD パラメータ (bond13_150=1.661, bond14_150=2.502, stiffness 等) は
cholesterol 等から得た **汎用既定値**であり、系ごとの最適値ではない。bead 粒子は
index ベース (P0, P1, …) で、命名 typing は範囲外。

## インストール

```bash
pip install -e '.[fragmenter]'   # fragmenter と同じ extras (rdkit-pypi)
```

## クイックスタート

### CLI

```bash
python -m abmptools.fragmenter.cg_segmenter build \
    --pdb cholesterol.pdb \
    --output-dir ./cg_segments \
    --target-mw 200
```

オプション:
- `--no-separate-rings`: ring を chain と同様に MW walk で切断
- `--no-atom-sharing`: fused ring の共有 atom を片方の segment にのみ寄せる
- `--no-absorb-substituent`: ring に attach する 1 heavy atom 置換基 (-OH, -NH2, -F 等) を別 chain segment にする
- `--h-only`: hetero CH3 cap rule を切って全境界 H cap に

#### FCEWS 入力生成 (`fcews` サブコマンド)

同じ segment 分割から **FCEWS** の入力 (`segment_data.dat` + monomer `<name>.xyz`)
を生成する。FMO フラグメントは atom を共有できないため、このサブコマンドは
内部で `allow_atom_sharing=False` を強制する (cap は使わず atom partition のみ利用)。

```bash
python -m abmptools.fragmenter.cg_segmenter fcews \
    --pdb cholesterol.pdb \
    --output-dir ./fcews_input \
    --target-mw 200 \
    --name cholesterol
# -> ./fcews_input/segment_data.dat  (FCEWS form, mode='FMO')
#    ./fcews_input/cholesterol.xyz   (monomer 座標、mol idx 順)
```

詳細は後述「[FCEWS export](#fcews-export-fcews-向け-segment_datadat-生成)」を参照。

### Python API

```python
from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

config = CGSegmenterConfig(
    pdb_path="cholesterol.pdb",
    target_mw=200,
    output_dir="./cg_segments",
)
sg = CGSegmenter.from_pdb(config)
result = sg.export()
print(f"{len(result.segments)} segments, {result.total_atoms_with_cap} atoms (incl. caps)")
for seg in result.segments:
    print(f"  seg {seg.segment_id} ({seg.kind}): {len(seg.atom_indices)} heavy, {len(seg.cap_atoms)} caps")
```

### Jupyter UI (`open_panel`)

`fragmenter.open_panel` と同様の interactive 編集 UI:

```python
from abmptools.fragmenter.cg_segmenter import (
    CGSegmenter, CGSegmenterConfig, open_panel,
)
config = CGSegmenterConfig(pdb_path="cholesterol.pdb", target_mw=200)
sg = CGSegmenter.from_pdb(config)
open_panel(sg)
```

UI 構成 (require `pip install -e '.[fragmenter,jupyter]'`):

| Widget | 機能 |
|---|---|
| **Show atom numbers** | SVG に heavy atom 番号 + shared atom には `*` 注記 |
| **Bold shared atoms** | 複数 segment に属する atom を太い黒縁取りで強調 |
| **SVG output** | 各 segment を別色 (赤/緑/青/オレンジ/紫/シアン/ピンク/ライム/...) で塗る、shared atom 黒縁取り |
| **Segments list** | 各 segment 行に `[Delete]` ボタン + 各 cap の `[cap[n] toggle]` ボタン (H ↔ CH3 切替) |
| **Move atom** | atom dropdown + 移動先 seg dropdown + `Keep in source` (shared option) checkbox + `Move` ボタン |
| **Re-segment** | `Target MW` 入力 + 3 flags (`Separate rings` / `Atom sharing` / `Absorb substituent`) + `Re-segment` ボタン (全 segment 上書き) |
| **Export** | `Export PDB + XYZ + JSON` ボタン |

`CGSegmenter` クラスの edit API:

```python
sg.move_atom(atom_idx, target_seg_id, shared=False)   # exclusive or shared
sg.toggle_cap(segment_id, cap_index)                  # H <-> CH3
sg.delete_segment(segment_id)
sg.re_segment(target_mw=..., separate_rings=..., allow_atom_sharing=...,
              absorb_single_substituent=...)
```

これらのいずれも cap atoms を自動再計算する (`_recompute_caps`)。

## アルゴリズム詳細

### 1. PDB ロードと連結成分分解 (`pdb_loader.py` 共通)

`abmptools.fragmenter.pdb_loader.load_pdb_molecules` を流用。現状は **単一分子 PDB** を
想定 (複数分子なら最初の連結成分のみ処理、`logger.warning` を出す)。

### 2. Ring segment 抽出 (`ring_detector.py`)

```
ring_info = mol.GetRingInfo()
for ring_atoms in ring_info.AtomRings():
    # 各 SSSR ring を 1 Segment に
    # fused ring の共有 atom は両 Segment の atom_indices に含める
    #  (allow_atom_sharing=True の場合; False なら若い ring 側に寄せる)
```

`absorb_single_substituent=True` の場合、ring atom に直接 attach する 1 heavy atom
置換基 (= -OH の O、-NH2 の N、-F の F 等で degree(heavy)=1 のもの) を近接 ring
segment に吸収する (`kind="ring_with_substituent"`)。

### 3. Chain split (`chain_splitter.py`)

ring に含まれない heavy atom を chain として抽出 → connected component に分割 → 各
component を `auto_split._atom_total_mw` を使った MW walk で `target_mw` 基準で切断。

### 4. Cap atom 付与 (`cap_attach.py`)

各 Segment の境界 atom (= segment 外に heavy neighbor を持つ atom) に cap atom を
追加する:

| 境界 atom の元素 | cap |
|---|---|
| C / halogen (F/Cl/Br/I) / その他 | **H** (1 atom) |
| **N / O / S / P** (`hetero_cap_methyl_elements` で指定) | **CH3** (central C + 3 H = 4 atoms) |

cap 位置: 元 bond の方向ベクトル × 結合長 (C-H=1.09 Å, C-C=1.54 Å, N-C=1.47 Å,
O-C=1.43 Å 等)。CH3 cap の 3 H は `methyl_hydrogen_positions` で tetrahedral
(109.47°) 配置。

**設計理由**: O や N が cut 面に出ると未飽和の lone pair が残り、後段の MD で
不要な水素結合スポット (擬似 H-bond donor/acceptor) になりやすいので、CH3 cap で
methyl ether / amine 化して安定化する。

### 5. 出力 (`exporter.py`)

```
output_dir/
├── seg_000.pdb           # PDB format (ATOM records + REMARK)
├── seg_000.xyz           # XYZ format
├── seg_001.pdb
├── seg_001.xyz
├── ...
└── segments.json         # 全 segment のサマリ
```

`segments.json` 構造:

```json
{
  "n_segments": 5,
  "total_atoms_with_cap": 80,
  "shared_atom_pairs": [
    [3, 0, 1],
    [8, 0, 1]
  ],
  "segments": [ {...}, {...}, ... ]
}
```

`shared_atom_pairs` の要素 `[atom_idx, seg_a, seg_b]` は atom が複数 segment に
属することを示す (fused ring の共有 atom 検証用)。

## データクラス

```python
@dataclass
class CapAtom:
    parent_atom_idx: int               # cap が付いた元 atom (segment 内)
    element: str                       # "H" or "C"
    position: Tuple[float, float, float]
    is_methyl_cap: bool = False        # True なら central C + 3 H (= 4 atoms)

@dataclass
class Segment:
    segment_id: int
    atom_indices: List[int]            # heavy atom (共有可能)
    cap_atoms: List[CapAtom]
    smiles: str = ""
    kind: str = "chain"                # "ring" / "chain" / "ring_with_substituent"
    ring_indices: List[int] = []

@dataclass
class CGSegmenterConfig:
    pdb_path: str
    output_dir: str = "./cg_segments"
    target_mw: float = 200.0
    separate_rings: bool = True
    allow_atom_sharing: bool = True
    hetero_cap_methyl_elements: List[str] = ["N", "O", "S", "P"]
    absorb_single_substituent: bool = True
```

## 用例: コレステロール (PubChem CID 5997)

ステロイド骨格は 4 つの fused ring (A, B, C, D) + iso-octyl tail + 3-OH。期待される
segmentation:

| Segment | 内容 | 共有 atom |
|---|---|---|
| ring A | cyclohexane + 3-OH (absorbed) | 2 atoms (with ring B) |
| ring B | cyclohexane | 2 atoms × 2 (with A and C) |
| ring C | cyclohexane | 2 atoms × 2 (with B and D) |
| ring D | cyclopentane | 2 atoms (with C) |
| tail | iso-octyl (target_mw で更に分割) | -- |

cap 配置: 各 ring 間の境界には共有 atom があるので cap 不要、ring D から tail への
boundary は H cap、tail 末端は H cap、3-OH の OH 部分は ring A に吸収済。

## EUD-E (側鎖ありポリマー) サンプル

`sample/fragmenter/eude_block55.pdb` (Eudragit E mimic block 共重合体、PMMA×5 + DMAEMA×5、
heavy=90、MW=1288.7):

```bash
python -m abmptools.fragmenter.cg_segmenter build \
    --pdb sample/fragmenter/eude_block55.pdb \
    --output-dir ./cg_eude \
    --target-mw 200
```

期待 segment 数: 7 程度 (MW 1289 / 200 ≒ 6.4)。CH3 cap は ester O や amine N が
切断面に来た所、H cap は C-C 境界。

## DPDgen export (OCTA COGNAC 入力生成)

CG segments から **[DPDgen](https://github.com/kojioku/dpdgen)** (Koji Okuwaki 作の DPD UDF
生成ツール、ABINIT-MP / FCEWS / COGNAC エコシステム) の入力 2 ファイルを生成する:

```bash
python -m abmptools.fragmenter.cg_segmenter dpdgen \
    --pdb cholesterol.pdb \
    --output-dir ./dpdgen_input \
    --monomer-name chol \
    --box 10 --total-num 100000 --step 100
```

または Python API:

```python
from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

config = CGSegmenterConfig(pdb_path="cholesterol.pdb", target_mw=200)
sg = CGSegmenter.from_pdb(config)
monomer_path, calc_sett_path = sg.export_dpdgen(
    output_dir="./dpdgen_input",
    monomer_name="chol",
    box_size=(10, 10, 10),
)
```

Jupyter UI からも `[Export DPDgen (monomer + calc_sett)]` ボタンで同じ操作が可能。

### 出力ファイル

```
output_dir/
├── {name}_monomer       # bond12 / bond13_150 / bond14_150 / angle13 (Python script)
└── {name}_calc_sett     # total_num_list / phys_param / ratio_list ... (Python script)
```

### bond ポテンシャル (距離制約)

| 種類 | 検出 | distance | stiffness |
|---|---|---|---|
| **bond12** (path 1, 隣接) | 共有 atom (fused boundary) または直接 atom bond | **0.86** / 0.60 (ring-ring) | 50 / 200 |
| **bond13_150** (1-skip ring) | fused ring chain BFS path 2 | **1.661** (cosine law, 150° 仮定) | 200 |
| **bond14_150** (2-skip ring) | fused ring chain BFS path 3 | **2.502** (4-ring chain extension) | 200 |

- bond12 distance 0.86: DPD 均等配置時のセグメント間距離 (default)
- bond12 ring-ring 0.60: cholesterol-like 特例 (ring-ring が密に詰まる場合)
- bond13_150 1.661 = 2 × 0.86 × sin(75°) (隣接 ring 内角 150° の余弦定理)
- bond14_150 2.502 = A-B-C-D 4 環直列で A-B-C 内角 150° + B-C-D 内角 165° による chain extension

> **boundary atom filter**: `_is_unique_atom()` で複数 segment に shared assigned された
> atom を has_direct check から除外。これにより cholesterol B 環内部 bond が A-C / B-D
> の direct bond と誤検出されるのを防ぐ。

### angle ポテンシャル (角度制約)

`bond12` graph で path length 2 のペア `(a, b, c)` を抽出し、 b を中央 segment として
`angle13` / `angle13data` に書き出す。

| triple の kind | eq (cognac 余角) | 実角度 | stiffness |
|---|---|---|---|
| 両端 ring + 中央 ring | **30** | 150° | 5.0 |
| cis 二重結合 周辺 (RDKit `BondType.DOUBLE`) | **60** | 120° | 5.0 |
| その他 (chain / mixed) | **0** | 180° (直線想定) | 5.0 |

> **cognac convention**: 平衡角は **余角** (180° - θ) で指定。150° なら 30、 120° なら 60、
> 直線 (180°) なら 0。 stiffness は全 angle で 5.0 (DPD 用にやや弱め設定)。

### オプション: 1-3 距離制約 (bond ポテンシャル)

`bond13_180` (dist 1.72) / `bond13_120` (dist 1.49) は **angle ポテンシャルで十分**
として default では生成しない (= monomer file の **コメントアウト template** として
出力)。距離制約も併用したい場合は monomer file の該当行を uncomment して値を埋める。

### cholesterol 出力例 (4 fused ring + 1 chain tail)

```python
bond12      = [[0,1], [1,2], [2,3], [0,4]]    # A-B / B-C / C-D / A-tail
bond12h     = [..., 0.60, 200 ..., 0.86, 50 ...]

bond13_150  = [[0,2], [1,3]]                  # A-C / B-D (1-skip)
bond13_150h = [[1, 0, 2, 1.661, 200, 0], ...]

bond14_150  = [[0,3]]                         # A-D (2-skip)
bond14_150h = [[1, 0, 3, 2.502, 200, 0]]

angle13     = [[1,0,4], [0,1,2], [1,2,3]]
angle13data = [
    [1, 0, 4, 0,  5.0],    # B-A-tail (mixed) → eq 0
    [0, 1, 2, 30, 5.0],    # A-B-C   (all ring) → eq 30
    [1, 2, 3, 30, 5.0],    # B-C-D   (all ring) → eq 30
]
```

### 自動生成しないもの

- **aij.dat** (相互作用パラメータ χ) → 本 exporter は生成しない。
  **`fcews-manybody` 等で別途計算済の aij.dat を `calc_sett` の `aij_file`
  パスに配置する。**
- **ratio_list / phys_param 詳細** → calc_sett に TODO コメント付き
  template として書かれるので、ユーザーが系構成に合わせて編集。
- **chobond_alkyl12 / chobond_steroid12 / chobond13_155 / chobond13_105** など
  cholesterol 特有の分類は自動推定しない (geometry を見て手動分類)。

### 出力ファイルの使い方

生成した `{name}_monomer` + `{name}_calc_sett` は、aij.dat と合わせて
**`abmptools.cg.dpd`** で DPD 入力 UDF / dpm に変換する（外部ツール不要・
abmptools 単体で完結）:

```bash
cd ./dpdgen_input
python -m abmptools.cg.dpd build-udf \
    --monomer chol_monomer --aij aij.dat --calc-sett chol_calc_sett \
    --output chol_uin.udf --particle-names P0,P1,P2,P3,P4
# -> chol_uin.udf (Cognac DPD 入力 UDF) を生成 → OCTA COGNAC で DPD MD 実行
```

> **Note**: monomer の bond12 / angle13 等の format は「DPDgen」(Okuwaki 作の
> 外部 DPD UDF 生成ツール) 由来の名残だが、消費側は `abmptools.cg.dpd`
> (自前実装、外部 dpdgen / UDFManager 非依存)。本ツールは外部 dpdgen に
> 依存しない。詳細は [`docs/cg_dpd.md`](cg_dpd.md)。

## FCEWS export (FCEWS 向け segment_data.dat 生成)

DPDgen と並ぶもう一つの下流出力。**同じ segment 分割**を使い、FCEWS が
`abmptools.abinit_io.config_read` 経由で読む `segment_data.dat` (FMO フラグメント
定義) + monomer `.xyz` を生成する。cap は使わず、segment の **atom partition**
だけを利用する。

```python
from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

# FMO フラグメントは atom を共有できないので allow_atom_sharing=False 必須
config = CGSegmenterConfig(pdb_path="cholesterol.pdb", target_mw=200,
                           allow_atom_sharing=False)
sg = CGSegmenter.from_pdb(config)
result = sg.export_fcews(output_dir="./fcews_input", name="cholesterol")
# result['segment_data'], result['xyz'], result['entry']
```

### 出力形式と規約

1 化学種 = 1 entry の FCEWS form (CLI `fcews` / `export_fcews` / `build_fcews_segment_data`):

| field | 内容 |
|---|---|
| `name` | entry 名。**monomer `<name>.xyz` の basename と一致必須** (`mol_io.read_mol_name` がファイル名から名前を取る) |
| `mode` | `'FMO'` |
| `atom` | per-fragment の atom 数 (heavy + H) |
| `charge` | per-fragment の formal charge 合計 |
| `connect_num` | per-fragment の **BAA atom 数** |
| `connect` | 各 cut bond を **`[BDA_atom, BAA_atom]`** で記録した flat list (1-origin local) |
| `seg_info` | per-fragment の atom serial (heavy + H、1-origin local) |

atom 番号は **mol atom idx + 1** (monomer 内 local)。monomer `.xyz` を mol idx 順に
書くので `seg_info` の番号と xyz 行が一致する。

**connect の atom 順序は `[BDA, BAA]`** (ABINIT-MP の AJF fragment 接続情報の並び)。
**`connect_num` は BAA 数** = 各 cut bond の BAA を含む fragment にのみ加算される
(LOGManager の `fbaas` と同じ。AJF の `Frag.` 列 = BAA の所属 fragment)。よって
`connect` の **第2要素 (BAA) が `connect_num>0` の fragment 内**にある (FCEWS 手書き
test-data の nafion: `connect=[[4,7]]`, `connect_num=[0,1]` で BAA=7 が frag2 内、
に一致)。BDA / BAA の **役割** は abmptools の
`auto_split.decide_bda_baa_for_manual_cut` (C-X → C 側 BDA / C-C → 若い idx) を流用。

### 制約

- **atom 共有不可**: fused ring (cholesterol 等) は `allow_atom_sharing=False`
  で partition すること。共有が残っていると `build_fcews_segment_data` が
  `ValueError` を投げる。CLI `fcews` は自動で `False` を強制する
- 第一弾は `mode='FMO'` のみ。rigid cluster (`solv`) / `term`・`repeat`
  (oligomer) は対象外

### cholesterol 出力例

```
'name': 'cholesterol', 'mode': 'FMO',
'atom': [15, 10, 11, 13, 25],          # 計 74 atoms (C27H46O)
'charge': [0, 0, 0, 0, 0],
'connect_num': [1, 2, 2, 2, 0],        # BAA 数 (= 各 connect 第2要素の所属 fragment)
'connect': [[7, 9], [12, 17], [13, 14], [16, 21], [17, 18], [20, 25], [21, 22]],  # [BDA, BAA]
```

## テスト

```bash
pytest tests/cg_segmenter/ -v
```

36 テスト (dataclass roundtrip / e2e 7 / 編集 op 7 / DPDgen 11 [path hierarchy +
angle 検証含む] + CLI 系 + **FCEWS export 7** [`config_read` round-trip /
connect `[BDA,BAA]` 順序 / cholesterol partition / 共有検出エラー])、~0.6s で PASS。

## 関連ドキュメント

- [`docs/fragmenter.md`](fragmenter.md) — 姉妹モジュール `fragmenter` (FMO 用)
- [`docs/cg_dpd.md`](cg_dpd.md) — **下流**: 本モジュール出力 (`{name}_monomer` + `{name}_calc_sett`) を fcews `aij.dat` と組み合わせて Cognac DPD 入力 UDF / OCTA viewer dpm を生成する `cg.dpd` サブパッケージ (v1.26.0 候補)
- [`docs/fcews_segment_data.md`](fcews_segment_data.md) もしくは本ドキュメント「FCEWS export」節 — **下流**: 同じ segment 分割から FCEWS `segment_data.dat` + monomer `.xyz` を生成 (FMO χ パラメータ評価向け)
- [`docs/overview.md`](overview.md) — abmptools 全サブパッケージ一覧
- `abmptools/fragmenter/cg_segmenter/README.md` — モジュール構成早見表

## バージョン履歴

| version | 日付 | 変更 |
|---|---|---|
| 1.24.0 | (本ドキュメント新設) | `abmptools.fragmenter.cg_segmenter` サブモジュール初版 (ring/chain/cap/PDB+XYZ 出力) |
| (Unreleased) | — | **FCEWS export** 追加 (`fcews_export.py` / `export_fcews` / CLI `fcews`): 同じ segment 分割から FCEWS `segment_data.dat` (`mode='FMO'`, connect `[BDA,BAA]`) + monomer `.xyz` を生成 |
