# PVA amorphous — ポリマー非晶質系の構築

ポリビニルアルコール (PVA) 10-mer atactic を amorphous box に詰めて MD を流す
sample。`sample/amorphous/` の他の例が低分子中心なのに対し、**オリゴマー鎖を
packmol で詰める**ケースの実例になっている。

## System

| 項目 | 値 |
|---|---|
| 分子 | PVA 10-mer atactic (`CH3-(CH(OH)-CH2)9-CH(OH)-CH3`) |
| SMILES | `CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)C` |
| atoms / mol | 75 (incl. H), 31 heavy |
| MW / mol | 456.57 g/mol |
| 分子数 | 30 (= 約 2,250 atoms) |
| 目標密度 | 1.2 g/cm³ (PVA バルク ~1.25) |
| 力場 | OpenFF Sage 2.1.0 (`openff_unconstrained-2.1.0.offxml`) |
| 電荷 | AM1-BCC |
| T_high | **400 K** |
| T_prod | 300 K |

**`T_high` を 400 K にしている理由**: OpenFF 既定の 600 K は小さめの有機分子だと
超臨界側に飛んで系が気体化することがある。ポリマーでも安全側に倒して 400 K で
アニールしている。他の系を組む時もここは既定値をそのまま使わずに見直すこと。

## Build + MD

```bash
cd sample/amorphous/pva_amorphous
bash run_sample.sh
cd md && bash run_all.sh && python wrap_pbc.py
```

5-stage MD (~10-20 min on CPU 8 cores): EM → NVT (T_high) → NPT (T_high)
→ annealing → NPT (T_prod)。出力: `md/test_05_output_rec*.bdf`。

`run_sample.sh` は AM1-BCC のために AmberTools の `antechamber` / `sqm` を要求する。

## COGNAC UDF へ

```bash
cd md && python build_bdf.py     # xtc -> multi-record UDF (101 records)
```

変換の詳細は [`docs/gro2udf.md`](../../../docs/gro2udf.md)。

## 備考

- この系は水素結合解析の題材としても使っている。解析側のワークフローは
  本パッケージには含まれない (README の "Beyond this package" を参照)。
- OpenFF Sage が書く UDF の `Atom_Type_Name` は per-atom unique な `MOL0_X` になる
  (SMIRNOFF が atom type の概念を持たないため)。下流で atom type を前提にする処理を
  書く場合はここに注意。
