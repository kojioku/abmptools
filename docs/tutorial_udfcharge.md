# Tutorial: 単分子 UDF → バルク UDF への電荷転写 (`abmptools.udfcharge`)

OCTA/COGNAC のバルク系 UDF は同じ分子を多数含むが、 電荷が 0 (または力場
default) のまま、ということがある。 一方で正しい per-atom 電荷を持つ **単分子
UDF** を別途用意できる場合 (例: 量子化学 / FMO で求めた電荷を入れた 1 分子 UDF)。

このチュートリアルでは、 その単分子 UDF の電荷を **バルク UDF の同名分子すべて**
へ一括転写する手順を、 methanol を題材に示す。

リファレンス: [`docs/udfcharge.md`](udfcharge.md) / サンプル: `sample/udfcharge/`

---

## 1. 環境

`abmptools` 本体に加えて **OCTA 同梱の `UDFManager`** が import できること
(PyPI には無い。 OCTA / J-OCTA インストールに付属)。

```bash
python -c "from UDFManager import UDFManager; print('UDFManager OK')"
python -c "import abmptools.udfcharge; print('udfcharge OK')"
```

`udfcharge` 自体は `UDFManager` 以外に追加依存はない (numpy のみ)。

---

## 2. サンプルを動かす (methanol、 ~数秒)

```bash
cd sample/udfcharge

# 2.1 デモ UDF を生成
#   methanol_single_charged.udf  : MeOH × 1、 per-atom 電荷あり (template)
#   methanol_bulk_uncharged.udf  : MeOH × 8、 電荷なし (割り当て先)
python make_example_udfs.py

# 2.2 単分子 → バルクへ電荷転写
python -m abmptools.udfcharge \
    --template methanol_single_charged.udf \
    --bulk     methanol_bulk_uncharged.udf \
    --out      methanol_bulk_charged.udf -v
```

期待出力:

```
template : MeOH (n_atoms=6, net_charge=+0.0000)
assigned : 8/8 molecules
output   : methanol_bulk_charged.udf
```

methanol の電荷 (和 = 0): C `+0.145` / O `-0.683` / Ho `+0.418` / H×3 `+0.040`。
転写後、 バルク 8 分子すべての各 atom にこの電荷が入る。

---

## 3. ステップ解説

### 3.1 template から電荷を読む

```python
from abmptools.udfcharge import read_molecule_charges

tmpl = read_molecule_charges("methanol_single_charged.udf")  # 既定: 先頭分子
print(tmpl.mol_name)          # 'MeOH'
print(tmpl.n_atoms)           # 6
print(tmpl.charges)           # [0.145, -0.683, 0.418, 0.04, 0.04, 0.04]  [e]
print(tmpl.atom_type_names)   # ['c3','oh','ho','h1','h1','h1']
print(tmpl.net_charge)        # 0.0
```

単分子 UDF に複数分子があるときは `mol_name=` か `mol_index=` で選ぶ:

```python
tmpl = read_molecule_charges("ref.udf", mol_name="MeOH")
```

読み出しは `electrostatic_Site[].ES_Element / 18.224159264` で [e] に戻し、
`electrostatic_Site[].atom[0]` で原子に対応づける (ES site が無い atom は 0)。

### 3.2 バルクへ割り当て

```python
from abmptools.udfcharge import assign_charges_to_bulk

res = assign_charges_to_bulk(
    "methanol_bulk_uncharged.udf", tmpl,
    "methanol_bulk_charged.udf",       # 省略時 <bulk>_charged.udf
)
print(res.n_molecules_assigned, "/", res.n_molecules_total)   # 8 / 8
print(res.skipped_indices)                                    # []
```

`Mol_Name` が template と一致する分子だけが対象。 割り当て前に **atom 数**
(必須) と **`Atom_Type_Name` 列** (`verify_atom_types`、 既定 True) を検証し、
一致したものに `electrostatic_Site` (`Type_Name="POINT_CHARGE"`、
`ES_Element=charge[e]×18.224159264`、 `atom[0]=index`) を書き込む。

出力は別ファイルなので、 入力バルク UDF は無改変。

### 3.3 結果を検証

```python
from UDFManager import UDFManager
E2Q = 18.224159264
u = UDFManager("methanol_bulk_charged.udf"); u.jump(-1)
for i in range(u.size("Set_of_Molecules.molecule[]")):
    na = u.size("Set_of_Molecules.molecule[].atom[]", [i])
    q = [u.get("Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element", [i, j]) / E2Q
         for j in range(na)]
    print(i, [round(c, 3) for c in q])
```

---

## 4. 実ユースケース: FMO/QM 電荷を MD バルク系へ

1. 量子化学 / FMO で 1 分子の partial charge を求め、 それを `electrostatic_Site`
   に入れた **単分子 UDF** を作る (または既存の電荷付き単分子 UDF を使う)。
2. MD で組んだ **バルク UDF** (同じ分子が数百〜数千、 電荷 0) を用意。
3. `udfcharge` で 1 → バルクへ転写:

```bash
python -m abmptools.udfcharge \
    --template qm_charges_single.udf \
    --bulk     md_bulk.udf \
    --out      md_bulk_qmcharge.udf -v
```

複数種の分子があるバルクには、 種ごとに template を変えて繰り返す
(`--mol-name` で対象を切り替え、 出力を次の `--bulk` に渡す):

```bash
python -m abmptools.udfcharge --template solute.udf  --bulk md_bulk.udf        --out step1.udf
python -m abmptools.udfcharge --template solvent.udf --bulk step1.udf          --out md_bulk_charged.udf
```

---

## 5. 失敗モードと対処

| 症状 | 原因 | 対処 |
|---|---|---|
| `atom 数 ... 不一致` で例外 | template とバルクの同名分子の atom 数が違う | 本当に同じ分子か確認。 別物なら `Mol_Name` を分ける |
| `Atom_Type_Name 列 ... 不一致` で例外 | atom 順序や型名が違う | 同一構造か確認。 型名だけ違う等で意図的に進めるなら `--no-verify-types` |
| `一致する分子がありません` | template の `Mol_Name` がバルクに無い | `--mol-name` でバルク側の名前に合わせる |
| 割り当てたのに電荷が 0 | numpy 値を直接 put した (silent 0 化) | 本モジュールは `float()` cast 済。 自前で put する場合は必ず `float()` |
| `UDFManager` が import できない | OCTA 未インストール | OCTA/J-OCTA を入れる (PyPI には無い) |

> **契約**: 電荷は **atom index 対応** (template atom *i* → 対象分子 atom *i*) で
> 割り当てる。 同名分子は同一の atom 順序であることが前提。 point charge
> (`electrostatic_Site`) のみ対象 (多極子は非対応)。
