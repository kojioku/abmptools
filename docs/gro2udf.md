# gro2udf — GROMACS .gro to COGNAC-UDF Converter

## 概要

`abmptools.gro2udf` は GROMACS が生成した `.gro` ファイル（位置・速度・セル情報）を
COGNAC-UDF の Structure レコードに書き戻すパッケージです。

逆方向（UDF → GRO）の変換は `abmptools.udf2gro` が担います。

---

## モジュール構成

```
abmptools/abmptools/gro2udf/
├── __init__.py        # Exporter を公開
├── __main__.py        # python -m abmptools.gro2udf 用
├── cli.py             # CLI エントリ（main()）
├── exporter.py        # 薄いオーケストレータ
├── gro_parser.py      # GRO ファイル → GROFrame（Parser 層）
├── gro_adapter.py     # GROFrame → AtomPosition / CellGeometry（Adapter 層）
└── udf_writer.py      # AtomPosition / CellGeometry → UDF レコード（Writer 層）
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

### CLI (後方互換ラッパ経由)

```bash
# gro2udf リポジトリの既存スクリプト経由（後方互換）
python gro2udf.py test.udf output.gro
# → test_groout.udf が生成される（カレントディレクトリ）
```

### CLI (パッケージ経由)

```bash
python -m abmptools.gro2udf test.udf output.gro
```

### Python API

```python
from abmptools.gro2udf import Exporter

Exporter().export("test.udf", "output.gro")
# → test_groout.udf が生成される（カレントディレクトリ）
```

---

## テスト手順

テスト入出力は `gro2udf/test/` にあります。

```bash
cd gro2udf/test/input
python ../../gro2udf.py test.udf output.gro
diff test_groout.udf ../output/test_groout.udf
# 差分なしであればテスト成功
```

または `abmptools` パッケージ経由:

```bash
cd gro2udf/test/input
python -m abmptools.gro2udf test.udf output.gro
diff test_groout.udf ../output/test_groout.udf
```

---

## 変換仕様

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

入力 UDF のベース名（拡張子除く）+ `_groout.udf` がカレントディレクトリに生成されます。

例: `test.udf` → `test_groout.udf`
