# acetaminophen amorphous — PubChem 由来の薬物非晶質系

公開薬物 **acetaminophen (paracetamol, C8H9NO2)** の amorphous MD を作る例。
分子は **PubChem から名前で 3D SDF を取得**する経路を使う (SMILES 入力の
`ketoprofen` / `pentane_benzene` に対する対照)。

## 系

| 項目 | 値 |
|---|---|
| 分子 | acetaminophen (PubChem name 取得、3D SDF) |
| 規模 | 64 分子 × 20 atoms = 1280 atoms、box 2.37 nm、density 1.2 g/cm³ |
| FF / 電荷 | OpenFF Sage 2.1.0 / AM1-BCC |
| MD | GROMACS 5-stage anneal、T_high=400 K、production 500 ps (5 ps stride → 101 records) |
| 官能基 | 2 級アミド (N-H / C=O)、フェノール -OH |

## 再現手順

```bash
cd <abmptools>/sample/amorphous/acetaminophen_amorphous
bash run_sample.sh                                          # build (~1 min, AM1-BCC)
export OMP_NUM_THREADS=4
cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 4' bash run_all.sh     # 5-stage MD
python wrap_pbc.py && python build_bdf.py                   # xtc -> UDF (101 rec)
```

`T_high` は 400 K。OpenFF 既定の 600 K だと小分子は気体化しうる。

## 備考

- この系は水素結合解析の題材としても使っている。解析側のワークフローは
  本パッケージには含まれない (README の "Beyond this package" を参照)。
- PubChem 取得の詳細は [`docs/amorphous.md`](../../../docs/amorphous.md)。
