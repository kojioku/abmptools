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

## テスト

```bash
pytest tests/cg_segmenter/ -v
```

11 テスト (dataclass roundtrip 2 + e2e 7 + 出力 1 + CLI 1)、~0.33s で PASS。

## 関連ドキュメント

- [`docs/fragmenter.md`](fragmenter.md) — 姉妹モジュール `fragmenter` (FMO 用)
- [`docs/overview.md`](overview.md) — abmptools 全サブパッケージ一覧
- `abmptools/fragmenter/cg_segmenter/README.md` — モジュール構成早見表

## バージョン履歴

| version | 日付 | 変更 |
|---|---|---|
| 1.24.0 | (本ドキュメント新設) | `abmptools.fragmenter.cg_segmenter` サブモジュール初版 (ring/chain/cap/PDB+XYZ 出力) |
