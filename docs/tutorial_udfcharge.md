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
python -m abmptools.udfcharge transfer \
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
python -m abmptools.udfcharge transfer \
    --template qm_charges_single.udf \
    --bulk     md_bulk.udf \
    --out      md_bulk_qmcharge.udf -v
```

複数種の分子があるバルクには、 種ごとに template を変えて繰り返す
(`--mol-name` で対象を切り替え、 出力を次の `--bulk` に渡す):

```bash
python -m abmptools.udfcharge transfer --template solute.udf  --bulk md_bulk.udf --out step1.udf
python -m abmptools.udfcharge transfer --template solvent.udf --bulk step1.udf   --out md_bulk_charged.udf
```

---

## 5. 形式電荷の復元 (`restore`)

MD は系を中性 (Σq=0) にして走らせるため、 元々整数の形式電荷を持つ分子の
電荷が `|q|` 比例で分散・中和されていることがある。 `restore` はその中和電荷
(B) と目標形式電荷 (S、 整数) から中和前の電荷 (A、 Σ=S) を復元する。

```bash
# 中和済み (Σq≈0) の 1 分子 UDF を、 形式電荷 +12 に復元
python -m abmptools.udfcharge restore \
    --udf system_neutral.udf --formal-charge 12 \
    --out system_q+12.udf -v
```

出力例:
```
molecule     : SYS (n_atoms=864)
input total  : +0.000000  (中和済み想定)
formal charge: +12  (λ=0.06593239)
output total : +12.000000
output       : system_q+12.udf
```

Python:
```python
from abmptools.udfcharge import restore_formal_charge
r = restore_formal_charge("system_neutral.udf", 12, "system_q+12.udf")
print(r.input_total, "->", r.output_total, "λ=", r.lam)
```

復元式 (B と S だけから一意):
`A_i = B_i/(1−λ)` (B_i>0) / `B_i/(1+λ)` (B_i<0)、
λ は `S·λ² + (P−N)·λ + (P+N−S) = 0` (P=Σ_{B>0}B, N=Σ_{B<0}B) の |λ|<1 の根。
`sample/udfcharge/restore_example.py` (methylammonium +1) と
`SI/A列再現方法.md` を参照。 補正が大きく正電荷が反転する (|λ|≥1) ケースは
一意復元できず `ValueError`。

---

## 6. 失敗モードと対処

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
