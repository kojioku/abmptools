# `abmptools.udfcharge` サンプル — 単分子 UDF → バルク UDF への電荷転写

単分子 UDF (電荷あり) の per-atom partial charge を、 バルク系 UDF の
**同名分子すべて**へ割り当てる。

## ファイル

| ファイル | 内容 |
|---|---|
| `make_example_udfs.py` | デモ UDF 生成スクリプト (methanol CH3OH) |
| `methanol_single_charged.udf` | 1 分子 `MeOH`、 per-atom 電荷あり (template) |
| `methanol_bulk_uncharged.udf` | `MeOH` × 8、 電荷なし (割り当て先) |
| `methanol_bulk_charged.udf` | 上記に電荷を転写した結果 |

methanol の電荷 (和 = 0): C `+0.145` / O `-0.683` / Ho `+0.418` / H×3 `+0.040`

## 実行

```bash
# 1) デモ UDF を生成 (single + bulk)
python make_example_udfs.py

# 2) 電荷を転写 (single → bulk の MeOH 全 8 分子)
python -m abmptools.udfcharge \
    --template methanol_single_charged.udf \
    --bulk     methanol_bulk_uncharged.udf \
    --out      methanol_bulk_charged.udf -v
```

出力:
```
template : MeOH (n_atoms=6, net_charge=+0.0000)
assigned : 8/8 molecules
output   : methanol_bulk_charged.udf
```

## Python API

```python
from abmptools.udfcharge import read_molecule_charges, assign_charges_to_bulk

tmpl = read_molecule_charges("methanol_single_charged.udf")   # 先頭分子から抽出
#   tmpl.mol_name / tmpl.charges ([e]) / tmpl.atom_type_names / tmpl.net_charge
assign_charges_to_bulk("methanol_bulk_uncharged.udf", tmpl,
                       "methanol_bulk_charged.udf")
```

## 仕組み (COGNAC UDF の電荷規約)

電荷は `Set_of_Molecules.molecule[i].electrostatic_Site[j]` に格納:

- `.Type_Name`  = `"POINT_CHARGE"`
- `.ES_Element` = partial charge [e] × **18.224159264** (COGNAC 内部スケール)
- `.atom[0]`    = この site が属する atom の (分子内 0-origin) index

`udfcharge` は template から `ES_Element / 18.224159264` で [e] を読み、
対象分子の各 atom へ `charge[e] × 18.224159264` を書き戻す。 atom 数と
`Atom_Type_Name` 列の一致を検証してから割り当てる (`--no-verify-types` /
`--non-strict` で緩和可能)。
