# abmptools.hbond — H-bond Analyzer for COGNAC Trajectories

非晶質 MD トラジェクトリ (OCTA COGNAC `.udf` / `.bdf`) から、官能基間 H-bond を
幾何条件で検出し、**官能基単位** (COOH / amide ごと) の役割比率を集計し、
gourmet で可視化できる UDF を生成するサブパッケージ。

**End-to-end チュートリアル**: IMC を例にした step-by-step は
[`tutorial_hbond_imc.md`](./tutorial_hbond_imc.md) を参照 (環境構築 →
CLI 実行 → 可視化 3 経路 → Jupyter UI → NMR 比較 → トラブルシューティング)。

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
| `<prefix>_summary.csv` (imc mode) | per-record の **官能基単位** dual/chain/single/free 数 + 比率 + mol-level legacy 数 |
| `<prefix>_pair_stats.csv` (generic mode) | per-record の donor_type × acceptor_type 集計 (n_hbonds, n_uniq_donors, n_uniq_acceptors, ratio_*_busy) |
| `<prefix>_classification.csv` (imc mode) | 全 carboxyl / amide ごとの (mol_index, group_index, role, partners) |
| `<prefix>_pairs.csv` | 検出された H-bond ペアの一覧 (距離・角度付き)。`kind` 列は imc では `cc/ca`、generic では `<donor_type>-><acceptor_type>` |
| `<prefix>.bdf` | 入力 BDF を Mol_Name 維持でコピー + 各 functional-group atom の `Attributes[]` に `Name='hbond' Value='Dual'/'Single'/'Free'/'Accept'` を append (**J-OCTA プリ描画用**、Attribute フィルタで官能基カテゴリ可視化可) |
| `<prefix>_colored.bdf` | Mol_Name リネーム + Draw_Attributes (`colorize_mode` が `molname` / `both` の時) |
| `<prefix>_action.bdf` + `<prefix>_show.act` | Mol_Name 維持 + autorun action (OCTA gourmet 用、`colorize_mode` が `action` / `both` の時) |
| `<prefix>_show.py` | autorun ラッパーなしのプレーン script (J-OCTA Python パネル用、`<prefix>.bdf` に対して Run、`colorize_mode` が `action` / `both` の時) |
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

### 3. 分類 (`classifier.py`) — 官能基単位、4 species

**v1.27.0 候補で per-functional-group ベース + 4 species 分類**。1 分子に複数の
COOH や amide がある場合、それぞれ独立に役割判定される。4 species 分類は
Yuan et al. (2015) *Mol. Pharm.* 12, 4518 (DOI 10.1021/acs.molpharmaceut.5b00705)
の amorphous IMC NMR deconvolution に対応:

| 集計対象 | 役割 | 判定基準 | NMR 対応 (Yuan 2015 Fig 5 / Table 1) |
|---|---|---|---|
| **COOH (carboxyl)** | **dual** | この COOH と相手 COOH の間で **両方向** に H-bond | ~179.3 ppm — cyclic dimer (58.5%) |
|  | **chain** | (dual ではない) cc 一方向 H-bond に参加 (donor または acceptor 側) | ~176.3 ppm — disordered chain (15.2%) |
|  | **single** | (dual / chain ではない) この COOH の OH-H が amide C=O に donate | ~172.4 ppm — COOH-amide (18.9%) |
|  | **free** | いずれもなし | ~170.4 ppm — free COOH (7.5%) |
| **amide C=O** | **accept** | この amide O が COOH から H-bond 受けている | (acceptor 側) |
|  | **free** | なし | (free 側) |

**`chain` の定義**: cyclic dimer ではない COOH-COOH 一方向 H-bond の片端
(donor or acceptor)。論文の "disordered chains" は「長さ様々の COOH 連鎖、
chain 末端、二量体より大きな ring」を含むとされ、`cc_relations` で逆方向が
無い arrow をすべて chain として拾う設計と一致。

**分子代表 role** (色付け用、`colorize_udf` が使用、優先度: dual > chain > single > free):

- 分子内に **いずれかの dual COOH** があれば → `dual`
- 上記でなく **いずれかの chain COOH** があれば → `chain`
- 上記でなく **いずれかの single COOH** があれば → `single`
- 上記すべて該当しない → `free`

**インドメタシン例** (1 mol = 1 COOH + 1 amide、Luzar-Chandler default、
T=450 K rec=900):

| species | NMR (Yuan 2015) | MD (本実装) |
|---|---|---|
| dual (cyclic dimer) | 58.5% | 8.0% (10/125) |
| chain (disordered) | 15.2% | 32.8% (41/125) |
| single (COOH-amide) | 18.9% | 30.4% (38/125) |
| free | 7.5% | 28.8% (36/125) |

amide: accept=49 / free=76 (39% / 61%)、mol-level は per-COOH と同じ
(1 mol = 1 COOH なので)。

比較プロットは `sample/hbond/imc_amorphous/plot_nmr_comparison.py` で生成可能
(出力: `output/imc_hbond_nmr_comparison.png`)。

### 4. UDF 色付け (`colorizer.py`)

色付けは **分子代表 role** ベースで Mol_Name リネーム + Draw_Attributes 追加:

```
Set_of_Molecules.molecule[i].Mol_Name  ← 'IMC_DUAL' / 'IMC_SINGLE' / 'IMC_FREE'
Draw_Attributes.Molecule[].{Mol_Name, color, transparency, radius}
```

| Role | 色 | transparency |
|---|---|---|
| dual | Red | 1.0 (不透明) |
| chain | Magenta | 1.0 (不透明) |
| single | Blue | 1.0 (不透明) |
| free | Gray | 0.3 (薄く) |

**重要**: GOURMET の Draw_Attributes color は `select` 型で 9 色名のみ
(Red/Green/Blue/Magenta/Cyan/Yellow/White/Black/Gray)。
RGBA tuple `[r, g, b, a]` は描画関数 (`sphere()` 等) でのみ使用可能で、
Draw_Attributes には書き込めない。

`transparency` は **1.0 = 不透明, 0.0 = 完全透明** (直観に反するので注意)。

### 5. J-OCTA プリ描画用コピー (`<prefix>.bdf`) + Attributes[] タグ

J-OCTA の「プリ描画」 (起動時の自動描画) は Mol_Name の値を内部で参照しており、
`molecular` 等の元名前を別の文字列にリネームすると分子が空表示になる。
そのため、Mol_Name 維持のコピーを `<prefix>.bdf` として併出する:

| ファイル | Mol_Name | 用途 |
|---|---|---|
| `<prefix>.bdf` | `molecular` 維持 | **J-OCTA プリ描画** (色付けなし、起動時自動描画される) |
| `<prefix>_colored.bdf` | `IMC_DUAL/SINGLE/FREE` リネーム | **OCTA ポスト描画** (gourmet または J-OCTA で show アクション実行時に 3 色塗分け) |

さらに `<prefix>.bdf` の **各 functional-group atom** の
``Set_of_Molecules.molecule[].atom[].Attributes[]`` に新規 entry を末尾 append:

| 対象 atom | Attribute Name | Value |
|---|---|---|
| carboxyl `c`/`o`/`oh`/`ho` | `hbond` | `Dual` / `Single` / `Free` |
| amide `c`/`o`/`n` | `hbond` | `Accept` / `Free` |
| その他 atom | (追記しない) | — |

既存の `Attributes[]` (J-OCTA が書く `Name='1' Value='molecular:<id>'` 等) は
**触らず末尾に push** するので round-trip 安全。J-OCTA の Attribute フィルタで
`hbond=Dual` の atom のみ表示など、色分けせずカテゴリ可視化が可能 (RGBA 色付け
は `<prefix>_action.bdf` + `<prefix>_show.{act,py}` の経路で別途行う)。

CLI: `--no-write-attributes` でこの append を抑制、`--attributes-name NAME` で
Name を変更可能 (default `hbond`)。

### 6. Python action 経由の per-functional-group 色付け (`<prefix>_action.bdf` + `<prefix>_show.act` + `<prefix>_show.py`)

GOURMET `Draw_Attributes` の schema は `Atom_Type[]` / `Bond_Potential[]` /
`Molecule[]` の 3 種類のみで、**個別 atom 単位の色付けや、Mol_Name を変えずに分子
役割を反映**する経路は公式には存在しない。これを補うため、`colorize_mode="action"`
(または `"both"`) で Python action ファイルを併出する経路を提供:

| ファイル | 内容 |
|---|---|
| `<prefix>_action.bdf` | 入力 BDF をコピーし、ヘッダの `Action:` フィールドに `<prefix>_show.act` を追記。**Mol_Name は変更しない** (J-OCTA プリ描画でも分子形状が描画される) |
| `<prefix>_show.act` | `autorun: showHbond()` 形式の action ファイル (各 carboxyl / amide の atom indices + role を埋込)。**OCTA gourmet 用** |
| `<prefix>_show.py` | 同じ描画ロジックを autorun ラッパーなしの **プレーン Python script** として出力。**J-OCTA Viewer の Python パネル用** (autorun 形式で落ちるケースに対応) |

仕組み:

- **OCTA gourmet 経路**: `<prefix>_action.bdf` を開くと、Action ヘッダ経由で
  `<prefix>_show.act` を自動 load。`autorun: showHbond()` が file open 時 +
  record 切替時に毎回実行され、`sphere()` で各官能基の atom に role 色を上書き描画
- **J-OCTA Viewer 経路** (autorun が落ちる場合): `<prefix>.bdf` (Mol_Name 維持
  uncolored copy) を開いて、Python パネルから `<prefix>_show.py` を `Load…`
  もしくはコピペ → `Run`。スクリプトは autorun ラッパー無しの直接実行形式
- 色: carboxyl dual=red、single=blue、free=gray、amide accept=cyan、amide free=gray
  (RGBA float タプル、`DEFAULT_ACTION_COLORS` で変更可)

**メリット**:
- Mol_Name 維持 → J-OCTA プリ描画と互換
- 1 分子内に複数 COOH / amide があっても官能基ごとに別色で描画される
- BDF binary section に触らないので UDFManager の round-trip 安全
- gourmet (autorun) と J-OCTA Python パネルの双方をカバー

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
--classify-mode {imc,generic}
                          imc (default): COOH 中心 4-species (dual/chain/single/free)
                          generic: donor-type x acceptor-type pair 統計
                          (PVA/peptide/アルコール等の任意系向け)
--no-element-fallback     element + bond-graph fallback を無効化 (strict mode)。
                          default ON で atom_type 不明な atom も自動 tag
                          (OpenFF SMIRNOFF UDF の `MOL0_X` 対応)
--no-colorize             color 出力を一切生成しない
--colorize-mode {molname,action,both}
                          molname (default): Mol_Name rename + Draw_Attributes;
                          action: Python action 経由 (Mol_Name 維持、J-OCTA 互換);
                          both: 両方併出
--no-copy-uncolored       <prefix>.bdf (Mol_Name 維持コピー) を生成しない
--no-write-attributes     <prefix>.bdf に hbond=Dual/Single/Free/Accept の
                          Attributes append を抑制
--attributes-name NAME    append する Attribute Name (default hbond)
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

## Generic mode (任意系向け) — v1.28 候補

COOH 中心の `imc` mode (dual/chain/single/free) ではなく、**donor type × acceptor
type の pair 統計** を出すモード。PVA (hydroxyl × hydroxyl_O)、peptide (amide_donor
× amide_O)、アルコール / ester 混合系等、COOH を持たない系で使う。

### 切替

```bash
python -m abmptools.hbond traj.bdf \
    --classify-mode generic \
    --donor-groups hydroxyl \
    --acceptor-groups hydroxyl_O
```

Jupyter UI からも `Mode:` dropdown で切替可能 (`open_panel(bdf_path)`)。

### 出力

| ファイル | 内容 |
|---|---|
| `<prefix>_pair_stats.csv` | per-record、donor_type × acceptor_type 集計 |
| `<prefix>_pairs.csv` | H-bond ペア (kind 列が `<donor>-><acceptor>` 形式) |
| `<prefix>.bdf` の Attributes | 各 atom に `hbond=Donor/Acceptor/Both` (active のみ書込、Candidate は skip) |
| `<prefix>_action.bdf` + `_show.act` | autorun overlay (gourmet 用)、atom role 色付け |
| `<prefix>_show.py` | Python パネル用 plain script (J-OCTA 用) |

### Element + bond-graph fallback (v1.28+)

OpenFF SMIRNOFF 等で UDF の ``Atom_Type_Name`` が per-atom unique
(``MOL0_0`` / ``MOL0_1`` / ...) になる場合、FF mapping が hit しない。
v1.28+ では ``fallback_tag_by_element`` が automatically:

- O atom: H と bond → ``hydroxyl_O``、無し → ``carbonyl_O``
- H atom: O と bond → ``hydroxyl_H``、N と bond → ``amide_H``
- N atom: C と bond → ``amide_N``
- C atom: ``carbonyl_O`` neighbor → ``carbonyl_C``

を assign する。default ON、CLI ``--no-element-fallback`` で strict
mode (FF mapping のみ) に切替可能。GAFF2 / OPLS / CHARMM 系では FF
mapping が hit するので fallback は no-op。

実機検証 (PVA 10-mer × 30、Luzar-Chandler default、101 records):
fallback ON で antechamber GAFF patch ありと完全同等 (rec=0 で
188 H-bonds、donor busy 62%、acceptor busy 60%)。

### atom role の色 (`DEFAULT_GENERIC_COLORS`)

| role | 色 (RGBA) | 意味 |
|---|---|---|
| `Donor` | `[1.0, 0.0, 0.0, 1.0]` (red) | atom が H-bond donor として active |
| `Acceptor` | `[0.0, 0.8, 0.8, 1.0]` (cyan) | atom が H-bond acceptor として active |
| `Both` | `[0.85, 0.0, 0.85, 1.0]` (magenta) | donor + acceptor 兼任 (例: alcohol O が donor で別 H-bond の acceptor) |
| `Candidate` | `[0.7, 0.7, 0.7, 0.25]` (faint gray) | donor/acceptor 候補だが active でない (描画 skip) |

### J-OCTA Attribute フィルタとの組合せ

`<prefix>.bdf` を J-OCTA で開いて Attribute フィルタで `hbond=Donor` / `Acceptor` /
`Both` の atom のみ選択表示可能。

## 既知の制限

1. **直交 cubic/orthorhombic box のみ対応** (三斜晶/モノクリニックは未対応)
2. FF mapping は組み込み 4 FF のみ。AMBER ff14SB 等は `add_mapping` で拡張可能
3. **多分子内多官能基**: 1 分子に複数の COOH/amide がある場合、官能基単位に
   独立に分類される (v1.27 候補)。**分子代表 role** (色付け用) は
   「分子内のいずれかの COOH が dual なら dual」優先度で決まる
4. 大規模 trajectory (1M+ records) でのメモリ効率は未最適化 (`H` 行列が大きくなる)

## サンプル

### `abmptools/sample/hbond/imc_amorphous/` (imc mode、COOH 4-species)

- `run_cli.sh`: CLI ワンライナー
- `run_notebook.ipynb`: Jupyter UI デモ
- `plot_nmr_comparison.py`: Yuan 2015 NMR vs MD 比較 plot (3-row layout)
- `README.md`: 期待値 (dual=10 / chain=41 / single=38 / free=36) + gourmet 表示方法

### `abmptools/sample/amorphous/pva_amorphous/` (generic mode、PVA, v1.28+)

- `input/pva_10mer.sdf`: PVA 10-mer atactic (rdkit 生成、75 atoms incl. H)
- `run_sample.sh`: OpenFF Sage build + 5-stage MD + xtc 生成 (AmberTools25 PATH 追加で AM1-BCC 有効化)
- `md/build_bdf.py`: 05_npt_final_pbc.xtc → multi-record UDF (101 frames) 変換
- `output/pva_hbond_pair_stats.csv`: hydroxyl→hydroxyl_O の per-record 集計
- 平均値 (101 records): n_hbonds=198.8, ratio_donor_busy=65.2%, ratio_acceptor_busy=61.7%
- **antechamber GAFF type patch 不要** (v1.28 の element + bond-graph fallback により)

## テスト

```bash
pytest abmptools/tests/hbond/ -v
```

- `test_functional_groups.py`: IMC BDF から carboxyl/amide が正しく検出されるか
- `test_hbond_detector.py`: 既知 geometry (8 角度パターン + PBC wrap) の単体テスト
- `test_classifier.py` (13 件): per-functional-group + per-mol 4-species 分類
  ロジック (no_hbonds / single_only / dual_pair / dual_precedence / one-way が
  chain / open trimer / dual priority over chain / multi_cooh_per_mol /
  partner 記録 / ratio_zero / back_compat alias 等)
- `test_func_tags.py` (9 件): 4 FF tag lookup + auto-detect + case-insensitive +
  validation (`add_mapping` の上書き / 競合検出)
- `test_lifetime.py` (8 件): continuous / intermittent / multiple pairs / gap
  tolerance / autocorr unbiased estimator / τ_HB integral / summary
- `test_amine_donor.py` (5 件): 合成 N-methylacetamide (GAFF2 + CHARMM36) +
  IMC tert amide 否定 + filter (include_amide / include_amine)
- `test_pair_type_stats.py` (5 件, v1.28+): generic mode の集約 (no_hbonds /
  single_pair / both_role / multiple_pair / unique_dedup)
- `test_element_fallback.py` (7 件, v1.28+): alcohol_OH / carboxyl pattern /
  amide_NH / preserves_existing / detect_carboxyls_unknown /
  detect_hydroxyls_unknown / disabled_no_groups
- `test_imc_regression.py` (5 件): IMC count ベースライン **(4-species)**
  carb dual / chain / single / free = 10 / 41 / 38 / 36 (sum=125)、
  amide accept / free = 49 / 76、hb_cc=31、hb_ca=50、
  `<prefix>.bdf` の Mol_Name 維持、`<prefix>_classification.csv` の per-group 行数

合計 **64 tests PASS** (2026-05-22 時点)。
