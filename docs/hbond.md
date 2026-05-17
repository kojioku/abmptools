# abmptools.hbond — H-bond Analyzer for COGNAC Trajectories

非晶質 MD トラジェクトリ (OCTA COGNAC `.udf` / `.bdf`) から、官能基間 H-bond を
幾何条件で検出し、各分子を **dual / single / free** の 3 グループに分類して gourmet
で可視化できる UDF を生成するサブパッケージ。

v1.26.0+ で **FF 抽象化 (GAFF2/OPLS-AA/CHARMM36/OpenFF) + 任意官能基対選択 +
multi-record lifetime/autocorrelation + secondary amide donor** に対応。

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
| `<prefix>_summary.csv` | per-record の dual/single/free count |
| `<prefix>_pairs.csv` | 検出された H-bond ペアの一覧 (距離・角度付き) |
| `<prefix>_colored.bdf` | gourmet で 3 色グループ可視化される UDF |
| `<prefix>_count.png` | count vs record プロット |

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

### 3. 分類 (`classifier.py`)

各分子に role を割り当て (優先度: dual > single > free):

- **dual**: 分子ペア (i, j) で **両方向** に COOH→COOH H-bond が成立
- **single**: COOH 供与体 → アミド受容体の H-bond に関与 (donor/acceptor 区別なし)
- **free**: いずれにも関与しない

### 4. UDF 色付け (`colorizer.py`)

GOURMET の `Draw_Attributes.Molecule[]` schema を利用:

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
--no-colorize             colored.bdf を生成しない
--no-plot                 count plot を生成しない
--quiet, -q               per-frame 出力を抑制
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
print("dual:", results[-1].n_dual_mols)
print("colored bdf:", paths["colored"])
```

### Jupyter UI

```python
from abmptools.hbond import open_panel
panel = open_panel("/path/to/imc.bdf")
```

- mol[0] の 2D 構造図 (RDKit) で検出されたカルボキシル (赤) / アミド (青) を表示
- 検出基準・record range・出力 prefix をウィジェットで設定
- Run ボタン → 統計表 + count plot がインラインで表示

## gourmet での可視化

```bash
gourmet output/imc_hbond_colored.bdf
```

開いたら `show` アクションを実行 → 3 色グループに塗り分けられた分子描画が表示
される。Picking で分子を選んでアクション実行も可能。

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
3. **多分子内多官能基**: 1 分子に複数の COOH/amide があると全て処理対象になり、
   分類は最も強い役割が優先される
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
- `test_classifier.py`: dual/single/free 分類ロジック
- `test_imc_regression.py`: IMC count ベースライン (dual=10, single=73, free=42)
