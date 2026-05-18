# sample/cg_dpd/cholesterol — `cg.dpd` 動作検証用サンプル

`abmptools.cg.dpd` (v1.26.0 候補) の検証用サンプル。 cholesterol を CG segment 化し、
そこから Cognac DPD 入力 UDF (R1) を生成した一式。

## 内容

| ファイル | 用途 |
|---|---|
| `aij.dat` | 5-segment fcews 統一形式 (mode='a', aii=25.0) の dummy aij、 P0-P4 間に hydrophobic-mismatch を仮定 |
| `chol_monomer` | cg_segmenter (`dpdgen` subcommand) で生成した monomer ファイル (Python script: bond12 / angle13 等) |
| `chol_calc_sett` | cg_segmenter で生成した calc_sett ファイル (box=(12,12,12), total_num=100000, step=100) |
| `chol_uin.udf` | **R1 出力**: `CGDpdBuilder.build_udf` で生成した Cognac DPD 入力 UDF (10.5 KB, 615 行) |
| `seg_*.pdb` / `seg_*.xyz` | cg_segmenter の per-segment 構造 (cap atom 込み) |
| `segments.json` | cg_segmenter のメタデータ |

## 再現方法

```bash
# 1) cg_segmenter で cholesterol を 5 segment + monomer + calc_sett に分解
python -m abmptools.fragmenter.cg_segmenter dpdgen \
    --pdb ../../cg_segmenter/cholesterol_rdkit.pdb \
    --output-dir . --monomer-name chol --box 12

# 2) aij.dat を fcews 統一形式 (Python 辞書 script) で用意 (本サンプルでは hand-craft)
#    実用時は fcews 出力をそのまま流用

# 3) R1: Cognac DPD 入力 UDF を生成
python -m abmptools.cg.dpd build-udf \
    --monomer chol_monomer --aij aij.dat --calc-sett chol_calc_sett \
    --output chol_uin.udf
```

## R1 (UDF) の構造 (生成された `chol_uin.udf` の section 概要)

| Section | 内容 |
|---|---|
| `\include{"cognac112.udf"}` | class 定義は J-OCTA install dir 経由で解決 (権利配慮、 abmptools は OCTA spec を持たない) |
| `\begin{header}` | `EngineType:"COGNAC"` / Ver112 / ProjectName / Comment |
| `Simulation_Conditions` | DPD dynamics (γ=20, λ=0.65)、 dt=0.05, total_steps=100, output_interval=100 |
| `Initial_Structure` | cell (12.0, 90°×3)、 粒子位置は空 (COGNAC の Position_Generation_From_Cognac で発生) |
| `Molecular_Attributes.Atom_Type[]` | P0..P4 (mass=1.0 each) |
| `Molecular_Attributes.Bond_Potential[]` | 4 entries (cholesterol bond12 = 4 pair の Harmonic) |
| `Molecular_Attributes.Angle_Potential[]` | 3 entries (cholesterol angle13 = 3 entry の Theta、 cognac 余角 convention) |
| `Molecular_Attributes.Pair_Interaction[]` | 15 entries (5 segment の対称 a パラメータ、 DPD type) |
| `Set_of_Molecules` | `"chol"` 1 entry (5 particle、 atom_list + bond_list + angle_list) |

## R2 (DPM、 B 案) について

R2 は user 提供の **空 dpm template** が必要なので、 sample dir には同梱しない。
J-OCTA の DPD Modeler で 1 回だけ空 dpm を作成し、 `CGDpdBuilder.build_dpm` に渡す:

```bash
python -m abmptools.cg.dpd build-dpm \
    --monomer chol_monomer --aij aij.dat \
    --template <user-provided>.dpm \
    --virtual-mom <user-provided>/Virtual.mom \
    --output-dir ./chol_proj --dpm-filename chol.dpm
```

詳細は [`docs/cg_dpd.md`](../../../docs/cg_dpd.md) 参照。

## 関連

- [`docs/cg_dpd.md`](../../../docs/cg_dpd.md) — `cg.dpd` リファレンス
- [`docs/cg_segmenter.md`](../../../docs/cg_segmenter.md) — 上流の CG segment builder
- `tests/cg_dpd/` — pytest 回帰テスト 30 件
