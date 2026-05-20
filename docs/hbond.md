# abmptools.hbond — H-bond Analyzer for COGNAC Trajectories

非晶質 MD トラジェクトリ (OCTA COGNAC `.udf` / `.bdf`) から、官能基間 H-bond を
幾何条件で検出し、**官能基単位** (COOH / amide ごと) の役割比率を集計し、
gourmet で可視化できる UDF を生成するサブパッケージ。

v1.26.0+ で **FF 抽象化 (GAFF2/OPLS-AA/CHARMM36/OpenFF) + 任意官能基対選択 +
multi-record lifetime/autocorrelation + secondary amide donor** に対応。
v1.27.0 候補 (per-functional-group classification) で「分子単位 1 役割」から
「官能基単位の役割比率」に統計指標を切り替え。J-OCTA プリ描画用に Mol_Name 維持
コピー (`<prefix>.bdf`) も併出する。

## 概要

ターゲット系: 非晶質インドメタシン (IMC) のような **カルボキシル基 + アミド基**
を持つ薬物分子の amorphous MD。13C CPMAS NMR で観測される 3 種類のカルボニル
環境 (環状二量体 / アミドへの単 H-bond / 自由) を計算側から定量する。

### 入力

- COGNAC UDF/BDF (record 1 つ以上、orthogonal cubic box)
- 分子は GAFF2 atomtype (`c`, `oh`, `ho`, `o`, `n`) を持つこと

### 出力

| ファイル | 内容 |
|---|---|
| `<prefix>_summary.csv` | per-record の **官能基単位** dual/single/free 数 + 比率 + mol-level legacy 数 |
| `<prefix>_classification.csv` | 全 carboxyl / amide ごとの (mol_index, group_index, role, partners) |
| `<prefix>_pairs.csv` | 検出された H-bond ペアの一覧 (距離・角度付き) |
| `<prefix>.bdf` | 入力 BDF の単純コピー (Mol_Name 維持、**J-OCTA プリ描画用**) |
| `<prefix>_colored.bdf` | Mol_Name リネーム + Draw_Attributes (**gourmet / OCTA ポスト描画用**) |
| `<prefix>_count.png` | count vs record プロット |
| `<prefix>_lifetime.csv` | (multi-record のみ) per-pair lifetime |
| `<prefix>_autocorr.{csv,png}` | (multi-record のみ) Luzar-Chandler 自己相関 |

## 検出アルゴリズム

### 1. 官能基検出 (`functional_groups.py`)

GAFF2 atomtype + bond graph をスキャン:

| 官能基 | 判定パターン |
|---|---|
| **carboxyl** (COOH) | `c` が `oh` (+`ho`) と `o` の両方と結合 |
| **amide** (C(=O)N) | `c` が `n` と `o` の両方と結合 |
| **hydroxyl** (free O-H) | `oh` が `ho` と結合 (carboxyl 除外オプション) |

Tertiary amide (N に H なし) は `tert=True` でマーキング。
インドメタシンは tert-amide (N が 2 つの aryl + 1 つの carbonyl で full substituted)。

### 2. H-bond 幾何判定 (`hbond_detector.py`)

直交 PBC 最短像法 (`minimum_image_vector`) で距離・角度を計算。

| 基準 | d(D-A) max [Å] | d(H-A) max [Å] | ∠(D-H-A) min [°] |
|---|---|---|---|
| **luzar-chandler** (default) | 3.5 | (制限なし) | 120 |
| **strict** | 3.5 | 2.5 | 150 |
| **custom** | ユーザ指定 | ユーザ指定 | ユーザ指定 |

参考: A. Luzar, D. Chandler, *Nature* 379, 55-57 (1996).

### 3. 分類 (`classifier.py`) — 官能基単位

**v1.27.0 候補で per-functional-group ベースに変更**。1 分子に複数の COOH や
amide がある場合、それぞれ独立に役割判定される (NMR の COOH-C / amide-C 信号
分離と直接対応)。

| 集計対象 | 役割 | 判定基準 | NMR 対応 (IMC 例) |
|---|---|---|---|
| **COOH (carboxyl)** | **dual** | この COOH と相手 COOH の間で **両方向** に H-bond | ~180 ppm |
|  | **single** | (dual ではない) この COOH の OH-H が amide C=O に donate | ~175 ppm |
|  | **free** | いずれもなし | ~165 ppm |
| **amide C=O** | **accept** | この amide O が COOH から H-bond 受けている | (acceptor 側) |
|  | **free** | なし | (free 側) |

**分子代表 role** (色付け用、`colorize_udf` が使用):

- 分子内に **いずれかの dual COOH** があれば → `dual`
- どれも dual でないが **いずれかの single COOH** があれば → `single`
- すべての COOH が free → `free` (amide が accept しているかは反映しない)

これは「分子の代表的な状態」であって、官能基単位の比率とは別。`summary.csv`
には両方が記録される (前者は `n_carb_dual/single/free`、後者は
`n_dual_mols/single_mols/free_mols`)。

**インドメタシン例** (1 mol = 1 COOH + 1 amide、Luzar-Chandler default、
T=450 K rec=900):

- carboxyl: dual=10 / single=49 / free=66 (8% / 39% / 53%)
- amide: accept=49 / free=76 (39% / 61%)
- mol-level: 上と同じ (1 mol = 1 COOH なので)

**v1.26.0 までは** "分子単位 1 役割" の集計で、COOH→amide H-bond の **両当事者**
(COOH 側 mol + amide 側 mol) を `single` にカウントしていた (上記 IMC で
single=73 になる旧仕様)。新仕様は COOH の状態を直接見るので、amide 側 mol
(COOH 持っていない or COOH free) は `free` に入る。

### 4. UDF 色付け (`colorizer.py`)

色付けは **分子代表 role** ベースで Mol_Name リネーム + Draw_Attributes 追加:

```
Set_of_Molecules.molecule[i].Mol_Name  ← 'IMC_DUAL' / 'IMC_SINGLE' / 'IMC_FREE'
Draw_Attributes.Molecule[].{Mol_Name, color, transparency, radius}
```

| Role | 色 | transparency |
|---|---|---|
| dual | Red | 1.0 (不透明) |
| single | Blue | 1.0 (不透明) |
| free | Gray | 0.3 (薄く) |

**重要**: GOURMET の Draw_Attributes color は `select` 型で 9 色名のみ
(Red/Green/Blue/Magenta/Cyan/Yellow/White/Black/Gray)。
RGBA tuple `[r, g, b, a]` は描画関数 (`sphere()` 等) でのみ使用可能で、
Draw_Attributes には書き込めない。

`transparency` は **1.0 = 不透明, 0.0 = 完全透明** (直観に反するので注意)。

### 5. J-OCTA プリ描画用コピー (`<prefix>.bdf`)

J-OCTA の「プリ描画」 (起動時の自動描画) は Mol_Name の値を内部で参照しており、
`molecular` 等の元名前を別の文字列にリネームすると分子が空表示になる。
そのため、Mol_Name 維持のコピーを `<prefix>.bdf` として併出する:

| ファイル | Mol_Name | 用途 |
|---|---|---|
| `<prefix>.bdf` | `molecular` 維持 | **J-OCTA プリ描画** (色付けなし、起動時自動描画される) |
| `<prefix>_colored.bdf` | `IMC_DUAL/SINGLE/FREE` リネーム | **OCTA ポスト描画** (gourmet または J-OCTA で show アクション実行時に 3 色塗分け) |

## 使い方

### CLI

```bash
python -m abmptools.hbond <bdf_path> \
    --out-prefix hbond_result \
    --criteria luzar-chandler \
    --mol-name IMC
```

オプション:

```
--criteria {luzar-chandler,strict,custom}
--d-da-max FLOAT          (custom mode のみ)
--d-ha-max FLOAT          (custom mode のみ)
--angle-min FLOAT         (custom mode のみ)
--record-start INT        開始 record (default 0)
--record-end INT          終了 record (exclusive, default -1 = all)
--mol-name STR            色付け時の Mol_Name prefix (default "IMC")
--no-colorize             <prefix>_colored.bdf を生成しない
--no-copy-uncolored       <prefix>.bdf (Mol_Name 維持コピー) を生成しない
--no-plot                 count plot を生成しない
--quiet, -q               per-frame 出力を抑制
--force-field FF          GAFF2/OPLS-AA/CHARMM36/OpenFF (default 自動判定)
--donor-groups            carboxyl,amide_donor,amine_donor,hydroxyl (CSV)
--acceptor-groups         carboxyl_O,amide_O,hydroxyl_O,ether_O (CSV)
--no-lifetime             lifetime / autocorrelation を計算しない
--gap-tolerance INT       intermittent lifetime の gap 許容 frame 数
--dt FLOAT                record 間時間 (autocorr のτ_HB 単位、default 1.0)
```

### Python API

```python
from abmptools.hbond import Analyzer, AnalyzerConfig

cfg = AnalyzerConfig(
    bdf_path="/path/to/imc.bdf",
    out_prefix="./output/imc_hbond",
    criteria_mode="luzar-chandler",
    base_mol_name="IMC",
)
analyzer = Analyzer(cfg)
analyzer.load()
results = analyzer.run()      # List[FrameResult]
paths = analyzer.write_outputs()

# Per-functional-group counts (primary)
cls = results[-1].classification
print(f"COOH dual/single/free: "
      f"{cls.n_carboxyls_dual}/{cls.n_carboxyls_single}/{cls.n_carboxyls_free}"
      f" of {cls.n_carboxyls}")
print(f"  ratios: {cls.ratio_carboxyl_dual:.2f}/"
      f"{cls.ratio_carboxyl_single:.2f}/{cls.ratio_carboxyl_free:.2f}")
print(f"amide accept/free: "
      f"{cls.n_amides_accept}/{cls.n_amides_free} of {cls.n_amides}")

# Mol-level representative (for color)
print(f"dual mols: {results[-1].n_dual_mols}")
print(f"colored bdf: {paths['colored']}")
print(f"uncolored bdf (J-OCTA pre-render): {paths['uncolored']}")
```

### Jupyter UI

```python
from abmptools.hbond import open_panel
panel = open_panel("/path/to/imc.bdf")
```

- mol[0] の 2D 構造図 (RDKit) で検出されたカルボキシル (赤) / アミド (青) を表示
- 検出基準・record range・出力 prefix をウィジェットで設定
- Run ボタン → 統計表 + count plot がインラインで表示

## OCTA での可視化

### gourmet (Linux) / J-OCTA ポスト描画

```bash
gourmet output/imc_hbond_colored.bdf
```

開いたら **左パネルの Python タブ**で `show` アクションのスクリプトを確認し、
`show.all("line", "mol", ...)` の 3 番目引数を **"molname"** に書き換えてから
`Run` する → Mol_Name (`IMC_DUAL/SINGLE/FREE`) と Draw_Attributes に基づき
3 色グループに塗り分けられた分子描画が表示される。Picking で分子を選んで
アクション実行も可能。

### J-OCTA プリ描画

`output/imc_hbond_colored.bdf` (Mol_Name リネーム済) は J-OCTA のプリ描画 (起動時の
自動描画) で **分子が空表示**になる (内部の Mol_Name lookup が壊れるため)。

プリ描画用には `output/imc_hbond.bdf` (Mol_Name 維持) を渡す:

```
File → Open Particle UDF → output/imc_hbond.bdf
```

色付けはなしだが、空表示にならず分子形状が確認できる。色付き解析は
`_colored.bdf` を別途 OCTA ポスト描画で開く。

## v1.26.0+ の追加機能

### 1. FF 抽象化 (`func_tags.py`)

GAFF2 / OPLS-AA / CHARMM36 / OpenFF の 4 FF に組み込み対応。FF 名を指定しない場合は
atom type の overlap で自動判定:

```python
from abmptools.hbond import detect_force_field, get_mapping
ff = detect_force_field(atom_types)         # → "GAFF2" / "OPLS-AA" / "CHARMM36"
mapping = get_mapping(ff)
```

| FF | 例 atom type → 機能タグ |
|---|---|
| GAFF2 | `c` → carbonyl_C, `oh` → hydroxyl_O, `hn` → amide_H |
| OPLS-AA | `opls_267` → carbonyl_C, `opls_268` → hydroxyl_O |
| CHARMM36 | `CC` → carbonyl_C, `NH1` → amide_N (sec/peptide), `H` → amide_H |

ユーザ拡張も可能:
```python
GAFF2.add_mapping("custom_type", "carbonyl_C")
```

### 2. 任意官能基対選択

donor/acceptor の組み合わせを CLI/Python API/Jupyter UI から指定:

| Donor group | 対象 |
|---|---|
| `carboxyl` (default) | COOH OH-H |
| `amide_donor` | secondary/primary amide N-H (peptide backbone 等) |
| `amine_donor` | amine N-H |
| `hydroxyl` | alcohol OH-H |

| Acceptor group | 対象 |
|---|---|
| `carboxyl_O` (default) | carboxyl C=O |
| `amide_O` (default) | amide C=O |
| `hydroxyl_O` | OH O |
| `ether_O` | ester/ether O |

```bash
# Secondary amide donor + amide acceptor (peptide H-bond network)
python -m abmptools.hbond traj.bdf --out-prefix peptide \
    --donor-groups amide_donor --acceptor-groups amide_O
```

### 3. Multi-record + Lifetime (`lifetime.py`)

Multi-record trajectory で per-pair lifetime + Luzar-Chandler 自己相関を計算:

| 量 | 関数 | 内容 |
|---|---|---|
| Continuous lifetime | `compute_lifetimes(per_rec, gap_tolerance=0)` | 連続存在区間の strict 集計 |
| Intermittent lifetime | `compute_lifetimes(per_rec, gap_tolerance=k)` | k フレームまでの gap を許容 |
| Autocorrelation C(t) | `compute_autocorrelation(per_rec, max_lag)` | Luzar-Chandler `<h(0)h(t)>/<h(0)>` (unbiased) |
| τ_HB | `integrate_autocorrelation(c, dt)` | ∫C(t)dt (台形則) |

CLI:
```bash
python -m abmptools.hbond traj.bdf --out-prefix prefix \
    --gap-tolerance 2 --dt 0.5     # dt = 0.5 ps/record の場合
```

出力 (multi-record 時のみ追加):
- `<prefix>_lifetime.csv` — per-pair の continuous_max/mean、intermittent_max/mean、occupancy
- `<prefix>_autocorr.csv` — lag × C(t)
- `<prefix>_autocorr.png` — autocorrelation plot

### 4. Secondary amide donor 対応

Tertiary (N に H なし) / secondary/primary (N に H あり) を `AmideGroup.tert` フラグで
区別。secondary amide の N-H は `detect_amine_donors()` 経由で donor 集合に追加可能。
IMC は tertiary なので donor にはならない (acceptor のみ)。

```python
from abmptools.hbond import detect_amine_donors
donors = detect_amine_donors(molecules, include_amide=True, include_amine=False)
# → AmineDonorGroup(mol_index, n_atom, h_atom, from_amide=True) のリスト
```

## 既知の制限

1. **直交 cubic/orthorhombic box のみ対応** (三斜晶/モノクリニックは未対応)
2. FF mapping は組み込み 4 FF のみ。AMBER ff14SB 等は `add_mapping` で拡張可能
3. **多分子内多官能基**: 1 分子に複数の COOH/amide がある場合、官能基単位に
   独立に分類される (v1.27 候補)。**分子代表 role** (色付け用) は
   「分子内のいずれかの COOH が dual なら dual」優先度で決まる
4. 大規模 trajectory (1M+ records) でのメモリ効率は未最適化 (`H` 行列が大きくなる)

## サンプル

`abmptools/sample/hbond/imc_amorphous/`:
- `run_cli.sh`: CLI ワンライナー
- `run_notebook.ipynb`: Jupyter UI デモ
- `README.md`: 期待値 + gourmet 表示方法

## テスト

```bash
pytest abmptools/tests/hbond/ -v
```

- `test_functional_groups.py`: IMC BDF から carboxyl/amide が正しく検出されるか
- `test_hbond_detector.py`: 既知 geometry (8 角度パターン + PBC wrap) の単体テスト
- `test_classifier.py`: per-functional-group + per-mol 分類ロジック (11 件)
- `test_imc_regression.py`: IMC count ベースライン
  (carb dual/single/free = 10/49/66、amide accept/free = 49/76、
  hb_cc=31、hb_ca=50、`<prefix>.bdf` の Mol_Name 維持、
  `<prefix>_classification.csv` の per-group 行数)
