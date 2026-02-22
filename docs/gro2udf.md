# gro2udf — GROMACS .gro to COGNAC-UDF Converter

## 概要

`abmptools.gro2udf` は GROMACS が生成した `.gro` ファイル（位置・速度・セル情報）を
COGNAC-UDF の Structure レコードに書き戻すパッケージです。

逆方向（UDF → GRO）の変換は `abmptools.udf2gro` が担います。

---

## モジュール構成

### udf-and-gro モード（既存）

```
abmptools/abmptools/gro2udf/
├── __init__.py        # Exporter / TopExporter を公開
├── __main__.py        # python -m abmptools.gro2udf 用
├── cli.py             # CLI エントリ（main()）
├── exporter.py        # 薄いオーケストレータ
├── gro_parser.py      # GRO ファイル → GROFrame（Parser 層）
├── gro_adapter.py     # GROFrame → AtomPosition / CellGeometry（Adapter 層）
└── udf_writer.py      # AtomPosition / CellGeometry → UDF レコード（Writer 層）
```

### --from-top モード（追加）

```
abmptools/abmptools/gro2udf/
├── top_parser.py          # TOP/ITP → TopRawData（Parser 層）
├── top_adapter.py         # TopRawData + GRO → TopModel（Adapter 層）
├── top_model.py           # 中間表現 dataclass 群
├── top_exporter.py        # TopModel → COGNAC UDF（Writer 層）
├── mdp_parser.py          # MDP → MdpParams（シミュレーション条件）
└── default_template.udf   # テンプレート未指定時のデフォルト UDF
```

### 設計方針

変換パイプラインは以下の 3 段階:

```
GROParser  →  GROFrame
GROAdapter →  List[AtomPosition], CellGeometry   ← SystemModel の dataclass を再利用
UDFWriter  →  UDF レコード（Steps, Time, Position, Velocity, Unit_Cell）
Exporter   →  上記の薄い協調役（フィルタ・ループ管理）
```

| 層 | クラス/モジュール | 責務 |
|---|---|---|
| Parser  | `GROParser`  | GRO ファイルをフレーム単位で読み込む |
| Adapter | `GROAdapter` | GROFrame → `AtomPosition` / `CellGeometry` |
| Writer  | `UDFWriter`  | UDF に Steps/Time/Position/Velocity/Cell を書く |
| Exporter| `Exporter`   | Parser→Adapter→Writer を協調させる |
| CLI     | `cli.py`     | 引数解析・`Exporter` 呼び出し |

**中間表現**: `abmptools.udf2gro.system_model` の `AtomPosition` / `CellGeometry`
を再利用します（新たな中間表現クラスは作りません）。

`gro2udf` 文脈での `AtomPosition` フィールドの解釈:
- `mol_id`  → 0-based UDF mol インデックス（`Structure.Position.mol[m]` の m）
- `atom_id` → 0-based atom インデックス（`atom[a]` の a）
- `vx/vy/vz` → [nm/ps]（UDFWriter が × 1000 して [m/s] で put）

---

## 使い方

### udf-and-gro モード（後方互換）

```bash
# gro2udf リポジトリの既存スクリプト経由（後方互換）
python gro2udf.py test.udf output.gro
# → test_groout.udf が生成される（カレントディレクトリ）

# abmptools パッケージ経由
python -m abmptools.gro2udf test.udf output.gro
```

### --from-top モード

GROMACS TOP + GRO ファイルから COGNAC UDF を新規生成します。

```bash
# 最小構成（デフォルトテンプレートを使用）
python -m abmptools.gro2udf --from-top system.top output.gro

# テンプレートを明示指定
python -m abmptools.gro2udf --from-top system.top output.gro \
    --template template.udf

# MDP ファイルを与えて NH-Q・Ewald カットオフを自動計算
python -m abmptools.gro2udf --from-top system.top output.gro \
    --mdp system.mdp

# 全オプション
python -m abmptools.gro2udf --from-top system.top output.gro \
    --template template.udf --mdp system.mdp --out result.udf
```

#### テンプレート解決順序

`--template` 未指定時、以下の優先順位でテンプレートを選択します。

1. `<top_stem>.udf`（TOP ファイルと同じディレクトリ・同じ stem）
2. パッケージ同梱の `default_template.udf`
   （`abmptools/abmptools/gro2udf/default_template.udf`）

### Python API

```python
from abmptools.gro2udf import Exporter, TopExporter

# udf-and-gro モード
Exporter().export("test.udf", "output.gro")

# --from-top モード（MDP なし）
TopExporter().export("system.top", "output.gro",
                     template_path="template.udf",
                     out_path="result.udf")

# --from-top モード（MDP あり → Q・Ewald カットオフ自動計算）
TopExporter().export("system.top", "output.gro",
                     template_path="template.udf",
                     out_path="result.udf",
                     mdp_path="system.mdp")
```

---

## テスト手順

### udf-and-gro モード

```bash
cd gro2udf/test/input
python -m abmptools.gro2udf test.udf output.gro
diff test_groout.udf ../output/test_groout.udf
```

### --from-top モード

```bash
cd abmptools/sample/gro2udf/gro_top_mode/input

# テンプレートなし（デフォルトテンプレート使用）
python -m abmptools.gro2udf --from-top test.top output.gro \
    --mdp test.mdp --out output_fromtop.udf

diff <(sed 's/\r//g' output_fromtop.udf | grep -v "^$") \
     <(sed 's/\r//g' ../output/output_fromtop_ref.udf | grep -v "^$")
```

---

## 変換仕様（--from-top モード）

### Nose-Hoover 熱浴質量 Q

UDF フィールド: `Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q`

**単位**: amu·Å²

**計算式**:

```
Q = g · k_B · T · τ²
```

| 変数 | 意味 | 値・出典 |
|---|---|---|
| `g` | 自由度 | `3·N − 3`（N = 全原子数、並進 COM 3 自由度を除く） |
| `k_B` | ボルツマン定数（COGNAC 内部単位） | `0.83144626 amu·Å²/(ps²·K)` |
| `T` | 設定温度 [K] | MDP の `ref_t`（未指定時 300.0 K） |
| `τ` | 緩和時間 [ps] | MDP の `tau_t`（未指定時 0.1 ps） |

**k_B の単位変換**:

```
k_B [amu·Å²/(ps²·K)]
    = 1.38064852×10⁻²³ J/K
      / (1.66054×10⁻²³ J/(amu·Å²/ps²))
    = 0.83144626
```

> MDP を指定しない場合は `T=300.0 K, τ=0.1 ps` のデフォルト値が使用されます。
> 正確な Q を得るには `--mdp` オプションで MDP を渡してください。

---

### Ewald 静電相互作用デフォルト値

`--from-top` モードは以下の Ewald パラメータを UDF に自動設定します。

| UDF パス | 設定値 | 備考 |
|---|---|---|
| `Interactions.Electrostatic_Interaction[0].Name` | `"POINT_CHARGE"` | 固定 |
| `Interactions.Electrostatic_Interaction[0].Algorithm` | `"Ewald"` | 固定 |
| `Interactions.Electrostatic_Interaction[0].Scale_1_4_Pair` | `0.83333333333333` | 5/6（AMBER 慣例） |
| `Interactions.Electrostatic_Interaction[0].Ewald.Dielectric_Constant` | `0.0` | 固定 |
| `Interactions.Electrostatic_Interaction[0].Ewald.R_cutoff` | 下記フォーミュラで計算 **[Å]** | COGNAC 内部単位 |
| `Interactions.Electrostatic_Interaction[0].Ewald.Ewald_Parameters` | `"Auto"` | 固定 |

> `Simulation_Conditions.Calc_Potential_Flags.Electrostatic` は `1` に設定されます。

#### Ewald R_cutoff の計算式

Deserno & Holm (J. Chem. Phys. **109**, 7678, 1998) の最適化条件（実空間誤差 ≈ 逆空間誤差）に基づく式:

```
R_cutoff = sqrt(11.5) / α × 10   [Å]   (COGNAC 内部長さ単位)

α = sqrt(π) × (5.5 N / V²)^(1/6) [nm⁻¹]   (V は nm³)
```

| 変数 | 意味 | 出典 |
|---|---|---|
| N | 全原子数 | TopModel.n_atoms_total |
| V | シミュレーションボックス体積 [nm³] | GRO 第1フレームのセル a×b×c |

GRO ファイルがない場合（空フレーム）は定数フォールバック値 `11.303883305209` Å を使用します。

---

## 変換仕様（udf-and-gro モード）

### フレームフィルタ

`Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps` が N の場合、
GRO ファイル内の `step % N == 0` のフレームのみ UDF レコードとして書き込まれます。

### 単位変換

| 量 | GRO単位 | UDF put 単位 |
|---|---|---|
| 位置 x/y/z | nm | `[nm]` |
| 速度 vx/vy/vz | nm/ps | `[m/s]`（×1000 変換） |
| セルサイズ a/b/c | nm | `[nm]` |
| 時刻 time | ps（`t=` フィールド） | `[ps]` |

### セル形式

- 3値（直交系）: a=box[0], b=box[1], c=box[2], α=β=γ=90°
- 9値（三斜晶）: GROMACS ベクトル表現から a, b, c, α, β, γ を計算

---

## 出力ファイル名

| モード | デフォルト出力ファイル名 |
|---|---|
| udf-and-gro | `<udf_stem>_groout.udf`（カレントディレクトリ） |
| --from-top | `<gro_stem>_fromtop.udf`（カレントディレクトリ） |
