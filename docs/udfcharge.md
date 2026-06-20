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

```bash
python -m abmptools.udfcharge \
    --template mol.udf \        # 電荷を持つ単分子 UDF
    --bulk     bulk.udf \       # 割り当て先 (同名分子が複数、 電荷なし)
    --out      bulk_charged.udf \   # 省略時 <bulk>_charged.udf
    [--mol-name NAME] [--mol-index I] \
    [--no-verify-types] [--non-strict] [-v]
```

| オプション | 説明 |
|---|---|
| `--mol-name` / `--mol-index` | template から抽出する分子の指定 (省略時 先頭分子) |
| `--no-verify-types` | `Atom_Type_Name` 列の一致検証を無効化 |
| `--non-strict` | 検証不一致を例外でなく warning + skip |

## Python API

```python
from abmptools.udfcharge import read_molecule_charges, assign_charges_to_bulk

tmpl = read_molecule_charges("mol.udf", mol_name="MeOH")
# MoleculeChargeTemplate: mol_name / n_atoms / charges([e]) /
#                         atom_type_names / atom_names / net_charge

res = assign_charges_to_bulk("bulk.udf", tmpl, "bulk_charged.udf")
# AssignResult: out_path / n_molecules_assigned / n_molecules_total /
#               n_atoms_per_mol / skipped_indices
```

## 割り当ての契約

- **atom index 対応**: template の atom *i* → 対象分子の atom *i*。 分子の atom
  順序が同じであることが前提。
- 割り当て前に **atom 数 (必須)** と **`Atom_Type_Name` 列** (`verify_atom_types`)
  を検証。 不一致時は `strict=True` で例外、 `False` で warning + skip。
- 出力は別ファイル (入力 bulk は無改変)。

## 制限

- point charge (`electrostatic_Site`) のみ対象。 多極子等は非対応。
- 同名分子は同一構造 (同じ atom 順序) を仮定。 構造が異なる同名分子は
  `verify_atom_types` で弾かれる。

## 関連

- `abmptools.gro2udf` / `abmptools.udf2gro` — 同じ電荷規約で UDF を読み書き
- `sample/udfcharge/` — methanol を題材にした end-to-end デモ
