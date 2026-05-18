# abmptools.cg.dpd — CG → DPD 系入力 (Cognac UDF / J-OCTA dpm) ビルダー

`abmptools.cg.dpd` (v1.26.0 候補) は、 [`cg_segmenter`](cg_segmenter.md) で生成した
CG segment + [fcews](https://github.com/okuwaki/fcews) の `aij.dat` (相互作用パラメータ
Python 辞書 script) から、 OCTA COGNAC で実行可能な DPD 入力ファイルを生成する。

## 3 ルート構成

```
 [cg_segmenter R0 (既存)]              [fcews]
   monomer + calc_sett                  aij.dat (Python 辞書)
        |                                  |
        +-----------+----------------------+
                    |    |
              [R1]  |    | [R2]
                    ↓    ↓
            *_uin.udf   *.dpm + monomer-lib/<seg>/Virtual.mom + #Message.txt
            (Cognac      (J-OCTA GUI 用、 編集後 UDF 出力)
             直接実行)
```

| Route | 入力 | 出力 | 機能 |
|---|---|---|---|
| **R0** (既存) | CG segment (PDB / RDKit Mol) | `{name}_monomer` + `{name}_calc_sett` | [`cg_segmenter`](cg_segmenter.md) で生成 |
| **R1** (新) | `{name}_monomer` + `{name}_calc_sett` + `aij.dat` | `*_uin.udf` (Cognac DPD 入力) | plain text writer、 UDFManager 非依存 |
| **R2** (新) | `{name}_monomer` + `aij.dat` (+ user template dpm + Virtual.mom) | `*.dpm` + `<dir>/monomer-lib/<seg>/Virtual.mom` + `#Message.txt` | **B 案**: user template の `\begin{data}` 5 ブロックを patch |

## 権利配慮の設計

J-OCTA は商用ソフトであり、 dpm の class 定義 (`\begin{def}` 800+ 行) は OCTA 公開
仕様には含まれない可能性があるため、 abmptools は以下の方針:

| File 種別 | abmptools の出力範囲 | 解決方法 |
|---|---|---|
| **Cognac UDF** (R1) | `\begin{header}` + `\begin{data}` のみ。 冒頭 `\include{"cognac112.udf"}` 1 行 | J-OCTA install dir 経由で UDF Manager が class 定義を resolve |
| **DPM** (R2) | `\begin{data}` 内の 5 top-level block のみを **patch** (template の def section は触らない) | **B 案** = user が J-OCTA で空 dpm を 1 回作成して abmptools に渡す |
| **Virtual.mom** | 全 segment dir に user 提供の同一ファイルを copy するのみ | abmptools は内容を生成しない |

## API

### Python (推奨)

```python
from abmptools.cg.dpd import CGDpdBuilder

# 1) cg_segmenter で先に生成
# 2) fcews で aij.dat を生成 (もしくは hand-craft)
builder = CGDpdBuilder.from_files(
    monomer="chol_monomer",
    aij="aij.dat",
    calc_sett="chol_calc_sett",   # R2 のみなら省略可
    particle_names=None,           # cg_segmenter の汎用ラベル "P0..Pn-1" を上書きしたい場合
    project_name="chol-DPD",
)

# R1: Cognac DPD 入力 UDF
udf = builder.build_udf("chol_uin.udf")

# R2: J-OCTA dpm + monomer-lib (要 user 提供 template)
dpm = builder.build_dpm(
    template="empty.dpm",                  # user が J-OCTA で 1 回作成
    output_dir="./chol_project",
    virtual_mom_template="Virtual.mom",    # user が J-OCTA Monomer Modeler で 1 回作成
    dpm_filename="chol.dpm",
)
```

### CLI

```bash
# R1: UDF 生成 (single monomer)
python -m abmptools.cg.dpd build-udf \
    --monomer chol_monomer --aij aij.dat --calc-sett chol_calc_sett \
    --output chol_uin.udf

# R2: DPM 生成 (B 案、 single monomer)
python -m abmptools.cg.dpd build-dpm \
    --monomer chol_monomer --aij aij.dat \
    --template empty.dpm --virtual-mom Virtual.mom \
    --output-dir ./chol_proj --dpm-filename chol.dpm
```

## 多 monomer 混合系 (C1)

cholesterol + water、 lipid + protein など複数 monomer を含む系は ``from_multi_files``
で構築する。 各 monomer の particle 名はそれぞれ独立に指定可能。

### Python API

```python
from abmptools.cg.dpd import CGDpdBuilder

builder = CGDpdBuilder.from_multi_files(
    monomer_specs=[
        {"monomer_file": "chol_monomer",
         "particle_names": ["A_Chol","B_Chol","C_Chol","D_Chol","Tail_Chol"]},
        {"monomer_file": "wat_monomer",
         "particle_names": ["W"]},
    ],
    aij="aij_mixed.dat",          # 全 6 segment (A_Chol..Tail_Chol + W) を含む
    calc_sett="mixed_calc_sett",
)
udf = builder.build_udf("chol_wat_uin.udf")
```

### CLI (JSON で渡す)

`monomers.json`:
```json
[
  {"monomer_file": "chol_monomer",
   "particle_names": ["A_Chol", "B_Chol", "C_Chol", "D_Chol", "Tail_Chol"]},
  {"monomer_file": "wat_monomer",
   "particle_names": ["W"]}
]
```

```bash
python -m abmptools.cg.dpd build-udf \
    --multi-monomer monomers.json \
    --aij aij_mixed.dat --calc-sett mixed_calc_sett \
    --output chol_wat_uin.udf
```

> ``--monomer`` (single) と ``--multi-monomer`` (multi) は **排他**、 どちらか必須。

### 整合性 warning

``from_multi_files`` は monomer の segment 名が aij.dat に含まれていない場合、
**logging.WARNING** で警告する (R1/R2 出力では `aii=25.0` default にフォールバック):

```
WARNING [...] CGDpdBuilder.from_multi_files: 2 segment(s) not in aij.dat:
['X_Chol', 'Y_Chol'] (R1/R2 出力で default aii=25 が fallback されます)
```

## データクラス (`models.py`)

| Class | 概要 |
|---|---|
| `AijMatrix` | fcews aij.dat の構造化表現。 `segments` / `pairs` / `mode` (a or chi) / `aii`。 `to_a_values()` で chi → a 自動変換 (Groot-Warren: `a = aii + chi/0.306`) |
| `MonomerSpec` | cg_segmenter `{name}_monomer` の構造化表現 (`bond12` / `bond12h` / `bond13_150*` / `bond14_150*` / `angle13` / `angle13data`) |
| `CalcSett` | `{name}_calc_sett` の構造化表現 (total_num / step / box / phys_param 等) |
| `DpdSystemSpec` | R1/R2 共通の system 全体仕様 (monomers + aij + calc_sett + project_name) |

## fcews aij.dat 形式 (fcews 統一)

```python
# mode 'a' (Groot-Warren a パラメータ直接)
aij = [['W','A', 55.0], ['W','W', 50.0], ['A','A', 50.0], ...]

# mode 'chi' (Flory-Huggins χ → a に変換)
chi = [['A_Asp','B_Glu', -4.108], ['A_Asp','C_Thr', -6.528], ...]
```

各要素 `[seg_i, seg_j, value]` の 3 要素 Python list。 `exec(open(path).read())` で読み込み、
`aij_dic['aij']` または `aij_dic['chi']` を抽出。 両方ある場合は `aij` を優先。

## R1 (UDF) の構造

`\include{"cognac112.udf"}` 1 行 + `\begin{header}` + `\begin{data}` の plain text 構造。

| Section | 内容 | 由来 |
|---|---|---|
| `Simulation_Conditions` | DPD dynamics (γ=20, λ=0.65, cutoff=1.0) + dt / total_steps / output_interval | calc_sett |
| `Initial_Structure` | cell (a, 90°×3) + Restart で粒子位置は空 | calc_sett.box_size |
| `Atom_Type[]` | 全 segment {name, mass=1.0} | spec.segment_names() |
| `Bond_Potential[]` | Harmonic、 `"seg_i-seg_j"` 名 | cg_segmenter bond12h |
| `Angle_Potential[]` | Theta、 cognac 余角 convention (eq=30/60/0) | cg_segmenter angle13data |
| `Pair_Interaction[]` | DPD type、 `{a, 4.5}` | fcews aij (chi 自動変換含む) |
| `Set_of_Molecules` | 各 monomer の atom_list / bond_list / angle_list | monomer + particle_names |

**注**: 粒子位置は **空** で出力し、 COGNAC の `Position_Generation_From_Cognac=1` で実行時に
random 配置する設計。 manual で初期位置を埋めたい場合は monomer file 後処理 or dpdgen
`camuslib.gen_init_pos_file` を別途利用。

## R2 (DPM、 B 案) の patch 対象 5 ブロック

`patch_dpm()` は user 提供 dpm template の `\begin{data}` 内 (最後の `\end{def}` 以降)
で以下 5 つの top-level field を **brace counting で識別 + 値を新規生成内容で置換**:

| Field | 内容 |
|---|---|
| `SegmentModel[]` | 各 segment 1 entry (radius=4.5, mass=aii=25.0, 色 default) |
| `SegmentPairModel[]` | 各 (i ≠ j) pair の DPD パラメータ (a 値は aij から) |
| `PolymerModel[]` | 各 monomer 1 entry (bond12 を polymer 内 bond として) |
| `DpdInput` | dt / total_steps / DPDBond[] (bond12) / DPDAngle[] (angle13、 cognac 余角) / solvent[] |
| `FcewsParam` | 空テンプレ (FCEWS 未連携、 user 後追加可) |

template の class 定義 (J-OCTA 商用 spec) は **触らずに温存**。

### dpm の field 識別ロジック (`_replace_dpm_field`)

dpm は入れ子 `\begin{def}` (header 内 def + main def) を持つので、 **最後の** `\end{def}`
以降の `\begin{data}` から field 検索を始める。 これで `\begin{def}` 内の型宣言
(例 `DpdInput:DPDInput`) と data section の actual value (`DpdInput:{...}`) を正しく区別。

## テスト

```bash
pytest tests/cg_dpd/ -v
```

30 テスト (~1.2s) — models + I/O round-trip / R1 UDF e2e / R2 DPM e2e / CLI / chi→a 変換。

`tests/cg_dpd/conftest.py` の fixture:
- `cholesterol_cg`: cholesterol_rdkit.pdb → cg_segmenter で chol_monomer + chol_calc_sett 生成
- `sample_aij_a_mode`: 5-segment 'a' mode aij (15 pair)
- `sample_aij_chi_mode`: 3-segment 'chi' mode aij (chi → a 変換テスト用)
- `dpm_template_path`: `man/octa/dpdfile-test/dpm-sample.dpm` (abmptools repo 外、 不在なら skip)
- `virtual_mom_template_path`: 同上 `Virtual.mom`

## サンプル

`sample/cg_dpd/cholesterol/` に cholesterol → R1 UDF の生成例一式 (再現コマンド含む)。

## Particle 名のマッピング (cg_segmenter → fcews aij)

cg_segmenter `dpdgen_exporter` は particle に `P0..Pn-1` の汎用ラベルを与えるだけで、
fcews の aij.dat 側で実際に使う segment 名 (例 `A_Asp` / `segA` / `W`) との **対応は
user が明示的に行う**。

### Python API

```python
from abmptools.cg.dpd import CGDpdBuilder

builder = CGDpdBuilder.from_files(
    monomer="chol_monomer", aij="aij.dat", calc_sett="chol_calc_sett",
    particle_names=["A_Chol", "B_Chol", "C_Chol", "D_Chol", "Tail_Chol"],
)
```

または直接 `assign_particle_names`:

```python
from abmptools.cg.dpd import read_monomer, assign_particle_names
mono = assign_particle_names(read_monomer("chol_monomer"),
                              ["A_Chol", "B_Chol", "C_Chol", "D_Chol", "Tail_Chol"])
```

### CLI

```bash
python -m abmptools.cg.dpd build-udf \
    --monomer chol_monomer --aij aij.dat --calc-sett chol_calc_sett \
    --output chol_uin.udf \
    --particle-names "A_Chol,B_Chol,C_Chol,D_Chol,Tail_Chol"
```

### 命名規約 (推奨)

| Segment 種別 | 例 | 備考 |
|---|---|---|
| Cholesterol 環 (4 環 fused) | `A_Chol` / `B_Chol` / `C_Chol` / `D_Chol` | 構造的に区別、 dpdgen sample 流儀 |
| Cholesterol tail (acyl chain) | `Tail_Chol` | 環と区別 |
| Water | `W` / `WAT` | 標準慣習 (dpdgen template も `W`) |
| Amino acid (peptide CG) | `A_Asp` / `B_Glu` / ... | fcews ssPalm/protein sample に対応 |

### 名前不一致の症状

`particle_names` と aij.dat の `segments` が完全一致しない場合:
- **R1 UDF**: `Pair_Interaction[]` の DPD pair 名前が aij と異なる → COGNAC 実行時
  lookup 失敗 (no-force or segfault)
- **R2 DPM**: `SegmentPairModel[]` で aij の値が反映されず、 default `aii=25.0` で fallback

整合性チェック (実装側で warning する用):

```python
seg_names = set(builder.spec.segment_names())
aij_segs = set(builder.spec.aij.segments)
missing = seg_names - aij_segs
if missing:
    print(f"WARNING: aij に未定義の segment: {missing}")
```

## 関連

- [`docs/cg_segmenter.md`](cg_segmenter.md) — 上流の CG segment builder (本パッケージの R0)
- `abmptools/cg/dpd/` — 実装本体 (~1435 行、 9 ファイル)
- `tests/cg_dpd/` — pytest 回帰テスト 30 件
- [DPDgen](https://github.com/kojioku/dpdgen) — Koji Okuwaki 本人作の DPD UDF generator (本実装の参考、 logic を移植 — 詳細 attribution は abmptools NOTICE 参照)
