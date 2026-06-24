# `abmptools.udfcharge` — OCTA/COGNAC UDF への per-atom 電荷割り当て

単分子 UDF (電荷あり) から per-atom partial charge を抽出し、 バルク系 UDF の
**同名分子すべて**へ同じ電荷を転写する小モジュール。 力場の電荷だけを別途
用意した 1 分子 UDF から、 大規模系へ一括反映したいとき (例: 量子化学/FMO で
求めた電荷を MD バルク系へ流し込む) に使う。

## なぜ必要か

OCTA/COGNAC のバルク UDF は同じ分子を多数含むが、 電荷が 0 のまま (あるいは
力場 default) になっていることがある。 一方、 単分子の参照 UDF には正しい
per-atom 電荷が入っている。 `udfcharge` はこの 1 分子→バルクの**電荷転写**を
1 コマンドで行う。

## COGNAC UDF の電荷規約

電荷は `Set_of_Molecules.molecule[i].electrostatic_Site[j]` に格納される
(gro2udf / udf2gro と共通):

| フィールド | 意味 |
|---|---|
| `electrostatic_Site[].Type_Name` | `"POINT_CHARGE"` |
| `electrostatic_Site[].ES_Element` | partial charge [e] × `18.224159264` (COGNAC 内部スケール) |
| `electrostatic_Site[].atom[0]` | この site が属する atom の分子内 0-origin index |

読み戻し: `charge[e] = ES_Element / 18.224159264`。
`Set_of_Molecules` は static record (`UDFManager.jump(-1)`) にある。

> **注意**: `UDFManager.put` は numpy float32 をサイレントに 0 化するため、
> ES_Element は必ず Python `float()` で渡す。

## CLI

2 サブコマンド: **`transfer`** (電荷転写) と **`restore`** (形式電荷の復元)。

```bash
# transfer: 単分子 UDF の電荷をバルク UDF の同名分子へ転写
python -m abmptools.udfcharge transfer \
    --template mol.udf \        # 電荷を持つ単分子 UDF
    --bulk     bulk.udf \       # 割り当て先 (同名分子が複数、 電荷なし)
    --out      bulk_charged.udf \   # 省略時 <bulk>_charged.udf
    [--mol-name NAME] [--mol-index I] \
    [--no-verify-types] [--non-strict] [-v]

# restore: 中和 (Σq≈0) された 1 分子 UDF の電荷を指定形式電荷に復元
python -m abmptools.udfcharge restore \
    --udf mol.udf --formal-charge 12 \   # 目標形式電荷 (整数)
    --out mol_q+12.udf \                  # 省略時 <udf>_q<±S>.udf
    [--mol-index I] [--mol-name NAME] [-v]
```

> サブコマンドを省略した旧形式 (`python -m abmptools.udfcharge --template ...`) は
> 後方互換で `transfer` として解釈される。

| オプション (transfer) | 説明 |
|---|---|
| `--mol-name` / `--mol-index` | template から抽出する分子の指定 (省略時 先頭分子) |
| `--no-verify-types` | `Atom_Type_Name` 列の一致検証を無効化 |
| `--non-strict` | 検証不一致を例外でなく warning + skip |

## Python API

```python
from abmptools.udfcharge import (
    read_molecule_charges, assign_charges_to_bulk, restore_formal_charge,
)

# (1) 転写: 単分子 → バルク
tmpl = read_molecule_charges("mol.udf", mol_name="MeOH")
# MoleculeChargeTemplate: mol_name / n_atoms / charges([e]) /
#                         atom_type_names / atom_names / net_charge
res = assign_charges_to_bulk("bulk.udf", tmpl, "bulk_charged.udf")
# AssignResult: out_path / n_molecules_assigned / n_molecules_total / ...

# (2) 復元: 中和 UDF (Σq≈0) → 形式電荷 S
r = restore_formal_charge("mol.udf", 12, "mol_q+12.udf")
# RestoreResult: out_path / mol_name / n_atoms / formal_charge /
#                input_total / output_total / lam
```

## 形式電荷の復元 (`restore`)

MD では系を中性 (Σq=0) にして走らせるため、 元々整数の形式電荷を持つ分子の
電荷が `|q|` 比例で分散・中和されて UDF に入っていることがある。 `restore` は
その中和電荷 (`B`) と目標形式電荷 (`S`、 整数) から、 中和前の電荷 (`A`、 Σ=S)
を復元する。

- **forward (中和)**: `B_i = A_i − S·|A_i| / Σ|A|`。 λ = S/Σ|A| とおくと
  `B_i = (1−λ)A_i` (A_i>0) / `(1+λ)A_i` (A_i<0)
- **reverse (本機能)**: `A_i = B_i/(1−λ)` (B_i>0) / `B_i/(1+λ)` (B_i<0)
- λ は `Σ A = S` の制約から `S·λ² + (P−N)·λ + (P+N−S) = 0`
  (P=Σ_{B>0}B, N=Σ_{B<0}B) の **|λ|<1 の根**として、 **B と S だけ**から求まる

符号が保存される (λ<1) 範囲でのみ一意復元できる。 補正が大きすぎて正電荷が
反転するケース (|λ|≥1) は `ValueError`。 詳細は `SI/A列再現方法.md`。

## 割り当ての契約

- **atom index 対応**: template の atom *i* → 対象分子の atom *i*。 分子の atom
  順序が同じであることが前提。
- 割り当て前に **atom 数 (必須)** と **`Atom_Type_Name` 列** (`verify_atom_types`)
  を検証。 不一致時は `strict=True` で例外、 `False` で warning + skip。
- 出力は別ファイル (入力 bulk は無改変)。
- **座標・トポロジーは無改変**: 更新するのは `electrostatic_Site` のみで、
  `Structure.Position` (座標) / `Unit_Cell` / `bond`・`angle`・`torsion` (結合) は
  そのまま保持される。 したがって座標 + 結合付きバルク UDF に電荷を載せても、
  OCTA viewer で**分子の形ごと**そのまま開ける (`sample/udfcharge/` は座標 +
  bond/angle/dihedral 付き)。

## 制限

- point charge (`electrostatic_Site`) のみ対象。 多極子等は非対応。
- 同名分子は同一構造 (同じ atom 順序) を仮定。 構造が異なる同名分子は
  `verify_atom_types` で弾かれる。

## 関連

- `abmptools.gro2udf` / `abmptools.udf2gro` — 同じ電荷規約で UDF を読み書き
- `sample/udfcharge/` — methanol を題材にした end-to-end デモ
