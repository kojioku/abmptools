# abmptools.cg.peptide — Martini 3 Peptide CG Builder

Martini 3 粗視化 (CG) でペプチドアモルファス系 (peptide in water box) を end-to-end
構築するビルダー。残基配列 (1-letter)・本数・box サイズを JSON / YAML で指定すると、
GROMACS MD 入力一式 (CG `.itp`, `.gro`, `topol.top`, `index.ndx`,
`mdp/{em,nvt,npt,md}.mdp`, `run/run.sh`) を `PeptideCGBuilder.build()` 1 呼び出しで
出力する。v1.18.0 で追加。

姉妹サブパッケージ:

- `abmptools.cg.membrane` (v1.19.0) — ペプチド-膜 PMF。本サブパッケージを Stage 1 で
  sub-call して CG ペプチド単体生成を委譲し、その上に bilayer + umbrella sampling を
  重ねる構成
- `abmptools.membrane` (AA、CHARMM36 / Lipid21 backend) — 全原子版のペプチド-膜 PMF

横断的な選択ガイドは [`docs/peptide_builders.md`](peptide_builders.md) を参照。

## ライセンス上のルール (Martini 3 配布物の取り扱い)

このサブパッケージは **subprocess only** 戦略で、外部ツールのソースを abmptools に
同梱・改変しない設計:

- ✅ **`vermouth-martinize`** (PyPI、Apache-2.0) — `martinize2` CLI を提供。
  abmptools は `subprocess.run(["martinize2", ...])` で呼ぶだけ
- ✅ **`gmx`** (LGPL) — solvate / genion / grompp / make_ndx / mdrun (subprocess)
- ✅ **`tleap`** (AmberTools、LGPL、推奨) — atomistic peptide PDB 生成 (任意、PATH 上に
  あれば自動で使う、無ければ extended-backbone fallback に切替)
- ⚠️ **Martini 3 force field `.itp`** (cgmartini.nl 配布、license 未明記) — abmptools は
  **未同梱**。ユーザーが各自 `martini_v300.zip` を取得して `ff/` に 3 ITP
  (`martini_v3.0.0.itp` / `_solvents_v1.itp` / `_ions_v1.itp`) を unzip する設計。
  `validate` サブコマンドが取得手順を表示
- ⚠️ **Martini 3 water box `.gro`** (`martini_v3.0.0_water.gro`) — cgmartini.nl は
  直接配布せず。`solvent_enabled=True` (default) のとき、`ff/` に置けば使う / 無ければ
  `gmx insert-molecules` で W bead を Martini 標準密度 (8.36 W/nm³、5 nm cubic で約
  1045 beads) で **auto-generate**

このルールにより、abmptools 本体は **MIT のまま、商用利用 OK** を維持できる
(GPL-2.0 ツールへ subprocess する mere aggregation は GPL に感染しない、FSF 公式見解)。

## 全体パイプライン (6 stage)

```
PeptideBuildConfig (JSON / YAML)
  ↓
PeptideCGBuilder.build():
  ↓ (1) atomistic PDB         tleap (推奨) or extended-backbone fallback
  ↓ (2) martinize2 -ff martini3001  -> CG ITP / PDB
  ↓ (3) gmx insert-molecules  -> packed.gro
  ↓ (4) gmx solvate           -> system_solv.gro (Martini W、auto-gen 可)
  ↓ (5) gmx grompp + genion   -> system_ions.gro (NaCl 中和 + 0.15 M)
  ↓ (6) MDP × 4 + index.ndx + run.sh
```

`run.sh` の MD 実行順 (`bash run/run.sh`):

```
em -> nvt -> npt -> md (production)
```

## API (主要 dataclass)

### `PeptideBuildConfig` (top-level)

```python
@dataclass
class PeptideBuildConfig:
    peptides: List[PeptideSpec]                 # 1 つ以上 (混合系 OK)
    box_size_nm: float = 0.0                    # cubic 自動 (>0 で適用)
    box_lengths_nm: Optional[List[float]] = None  # rectangular [x, y, z] (排他)
    solvent_enabled: bool = True
    neutralize: bool = True
    nacl_molar: float = 0.15
    temperature: float = 310.0
    seed: Optional[int] = None
    # Force field
    martini_itp_dir: str = ""                   # 空 -> output_dir/ff
    # External tools
    martinize2_path: str = "martinize2"
    gmx_path: str = "gmx"
    tleap_path: str = "tleap"
    output_dir: str = "."
    # MDP toggles + nsteps
    mdp_em: bool = True
    mdp_nvt: bool = True
    mdp_npt: bool = True
    mdp_md: bool = True
    em_steps: int = 50_000
    nvt_nsteps: int = 250_000                   # 5 ns at dt=20 fs
    npt_nsteps: int = 250_000                   # 5 ns
    md_nsteps: int = 5_000_000                  # 100 ns
    dt_fs: float = 20.0                         # Martini 3 標準
```

### `PeptideSpec` (sub dataclass)

```python
@dataclass
class PeptideSpec:
    name: str = ""                              # 必須、unique
    sequence: str = ""                          # 1-letter code (e.g. "KGG", "AAAAA")
    count: int = 1                              # 系内のコピー数
```

`sequence` は 20 種類の standard amino acid (ACDEFGHIKLMNPQRSTVWY) のみ受理。
非標準残基や cap (ACE/NME) はサポートしないので、必要なら全原子版 `abmptools.membrane`
を使うこと。

JSON / YAML I/O は `to_json` / `from_json` で対応。YAML 入力は PyYAML の有無で
動的に判定 (extras `[cg]` に含まれる)。

## CLI

```bash
# example config 取得
python -m abmptools.cg.peptide example > kgg.json

# config + 依存性検証 (martinize2 / gmx / tleap + 3 ITP)
python -m abmptools.cg.peptide validate --config kgg.json --ff-dir ./ff

# 6 stage build (MDP / topol.top / index.ndx / run.sh 生成)
python -m abmptools.cg.peptide build --config kgg.json --ff-dir ./ff -o ./out

# GROMACS で MD 実行
bash out/run/run.sh
```

## 設計判断とその根拠

### 1. atomistic PDB は tleap 推奨、無ければ extended-backbone fallback

martinize2 は atomistic PDB → CG mapping のため、まず all-atom 構造が必要。

- **tleap (推奨)**: AmberTools の `tleap` で `sequence` から 3D 構造 + sidechain heavy
  atom を生成。sidechain bead 座標が物理的に妥当
- **extended-backbone fallback**: tleap が PATH に無い場合、純 Python で beta-strand 風の
  伸びきり backbone を書き出す。ただし sidechain heavy atom が無いため、**芳香族残基
  (W / F / Y) で martinize2 が CG sidechain bead 座標を NaN にする可能性**あり

研究品質では tleap の使用を強く推奨。`validate` サブコマンドが tleap の有無をチェック
して警告を出す。

### 2. solvate template の auto-generate

Martini 3 water box `.gro` は cgmartini.nl が直接配布していない (チュートリアル
archive 等にしか含まれない)。これを必須にするとユーザーが詰まりやすい。

対処:

- ユーザーが `ff_dir` に `martini_v3.0.0_water.gro` を置いた場合 → そのまま使う
- 無い場合 → `gmx insert-molecules` で W bead を標準密度 (8.36 W/nm³) で **auto-generate**

ユーザー側の負担を減らしつつ、自前 water box を持つ熟練者の override も許す。

### 3. 系混合は 1 つの config で複数 `PeptideSpec` を渡す

```json
{
  "peptides": [
    {"name": "kgg", "sequence": "KGG", "count": 5},
    {"name": "ala5", "sequence": "AAAAA", "count": 3}
  ]
}
```

`gmx insert-molecules` が複数種を順次 packing。`molecules/<name>/` ディレクトリが
peptide ごとに分かれて生成される。

### 4. MDP は 4 種を独立 toggle

`mdp_em` / `mdp_nvt` / `mdp_npt` / `mdp_md` を個別 bool で制御。例えば solvent_enabled=False で
ペプチドのみ build したい場合 (cg.membrane が Stage 1 で sub-call する用途) は
`mdp_*=False` で MD 関連を全部 skip できる。

```python
# cg.membrane が cg.peptide を sub-call する設定例
sub_config = PeptideBuildConfig(
    peptides=[PeptideSpec(name="kgg", sequence="KGG")],
    solvent_enabled=False,    # 膜系では cg.membrane 側で solvate
    neutralize=False,
    mdp_em=False, mdp_nvt=False, mdp_npt=False, mdp_md=False,
    output_dir="_cg_peptide_subbuild",
)
```

### 5. ATM (atomistic + martini) ハイブリッドはサポートしない

`PeptideSpec.sequence` から始める one-shot ビルドのみ。事前に martinize2 を回した
ITP / CG PDB を流用する `cg_pdb_path` のような escape hatch は v1.18 では未提供
(必要になれば cg.membrane 側の `PeptideMembraneSpec.cg_pdb_path` 同等を後付け予定)。

## 既知の caveat / 改善余地

| 問題 | 根本原因 | 対処 / 緩和策 |
|---|---|---|
| **芳香族残基 (W/F/Y) で sidechain bead 座標 NaN** | extended-backbone fallback の構造的限界 | tleap を install (`mamba install -c conda-forge ambertools`) |
| **cap (ACE / NME) 非対応** | Martini 3 の cap 取扱いはユーザー側で `martinize2` 引数調整が必要 | AA 版 `abmptools.membrane` を使う |
| **non-standard residue (HIP / CYM / D-amino acid 等) 非対応** | `PeptideSpec` の validation で reject | AA 版 `abmptools.membrane` または直 `martinize2` |
| **Martini 3 small molecule (CGenFF 不要なので利点) は別 channel** | `cg.peptide` はペプチド専用 | v1.20+ で `cg.smallmol/` (Auto-Martini) を予定 |
| **disulfide bond / GO-model / EN-model 非対応** | martinize2 の追加引数で対応可能だが、本ビルダーはデフォルト elastic network (EN) の有無を切り替えるオプション未提供 | martinize2 を直接呼んで生成した ITP を `cg.membrane.PeptideMembraneSpec.cg_itp_path` に渡す |

## デフォルト数値の理由

| パラメータ | 値 | 由来 |
|---|---|---|
| `dt_fs` | 20.0 | Martini 3 標準 (AA の 10×) |
| `temperature` | 310.0 K | 生体系の慣習 |
| `nacl_molar` | 0.15 | 生理食塩水濃度 |
| `em_steps` | 50,000 | Martini 3 標準 (CG なので AA より steep が必要少ない場合あり) |
| `nvt_nsteps` | 250,000 | 5 ns at dt=20 fs |
| `npt_nsteps` | 250,000 | 5 ns |
| `md_nsteps` | 5,000,000 | 100 ns production の典型値 (CG 換算で AA 約 400-1000 ns 相当) |
| `box_size_nm` | 0.0 | 0 で auto-fit (peptide PDB の minimum BB を見て padding 1 nm 追加) |

## ファイル構成

```
abmptools/cg/peptide/
├── __init__.py                # public API: PeptideCGBuilder, PeptideBuildConfig, PeptideSpec
├── __main__.py                # `python -m abmptools.cg.peptide`
├── _subprocess.py             # subprocess wrapper + logger + ensure_dir
├── builder.py                 # PeptideCGBuilder.build() (6 stage)
├── cli.py                     # build / validate / example
├── forcefield_check.py        # 3 ITP + martinize2 + gmx + tleap 検証
├── martinize_runner.py        # subprocess wrapper for martinize2
├── mdp_templates.py           # CG MDP renderer (em/nvt/npt/md)
├── models.py                  # PeptideBuildConfig + PeptideSpec
├── peptide_atomistic.py       # tleap CLI + extended-backbone fallback
├── system_packer.py           # gmx insert-molecules + genion glue
├── top_writer.py              # topol.top renderer (3 ITP + peptide ITP includes)
├── water_box.py               # martini_v3.0.0_water.gro auto-gen
└── README.md                  # ユーザー向け Quick Start
```

## 出力レイアウト

```
output_dir/
├── config.json                # 入力 config の dump (再現用)
├── topol.top
├── packed.gro
├── system_solv.gro
├── system_ions.gro            # 最終座標 (build() 返値の "gro")
├── ions.{mdp,tpr}
├── index.ndx
├── molecules/
│   └── <name>/
│       ├── <name>_atomistic.pdb
│       ├── <name>_cg.pdb
│       └── <name>.itp
├── mdp/
│   ├── em.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   └── md.mdp
└── run/
    └── run.sh                 # bash run/run.sh で MD 4 段階を順次実行
```

## 関連リソース

- [`abmptools/cg/peptide/README.md`](../abmptools/cg/peptide/README.md) — Quick Start / 入力 schema 概要
- [`docs/cg_membrane.md`](cg_membrane.md) — ペプチド-膜 PMF (本サブパッケージを sub-call)
- [`docs/peptide_builders.md`](peptide_builders.md) — 3 系統 (AA membrane / CG peptide / CG membrane) の選択ガイド
- [`docs/membrane.md`](membrane.md) — AA 版のペプチド-膜 PMF (CHARMM36 / Lipid21)

## バージョン履歴

| version | 日付 | 主な変更 |
|---|---|---|
| 1.18.0 | 2026-05-04 | 新設。Martini 3 ペプチド CG ビルダー、6 stage pipeline、tleap fallback、water box auto-gen |
