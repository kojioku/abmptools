# `insulin_aggregation_10mM` — 論文 insulin の論文準拠系

論文 Hossain 2023 の **insulin** (51 残基 = A 鎖 21 + B 鎖 30、 disulfide × 3、
natural L-AA only) を 10 mM 濃度 (= peptide × 6 in 12 nm cubic box) で再現。
RCSB **PDB 2G4M** を入力 (論文も同 PDB を使用)。

## 論文との対応

| 項目 | 論文 insulin | 我々 |
|---|---|---|
| PDB | 2G4M (RCSB) | **2G4M** ✅ |
| 残基 | A21 + B30 = 51、 natural AA | 同 ✅ |
| disulfide bond | A6-A11、 A7-B7、 A20-B19 (3 本) | 同 (OpenFF Topology で declare) |
| n_copies | 6 (10 mM) | **6** ✅ |
| box | 10 nm cubic | **12 nm** (insulin は大きい → 12 nm) |
| 力場 | CHARMM36 + CGenFF | **ff14SB + GAFF2 (商用 OK)** |

box を 10 → 12 nm に拡大した理由: insulin × 6 = ~770 atoms × 6 = ~4,600 peptide
heavy atoms、 caprate × 32 + TCH × 2 で box 余裕 + 0.15 M NaCl 調整。

## 期待 atom 数 + wall

| 項目 | 値 |
|---|---|
| atom 数 (推定) | ~110,000 (water + peptide + small mol) |
| wall (RTX 4070 Ti) | ~13-15 時間 / 100 ns |

## Build + run

```bash
cd sample/formulation/insulin_aggregation_10mM
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/insulin_10mM

cd /tmp/insulin_10mM
PATH=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256:$PATH \
NT=8 MDRUN_OPTS="-nb gpu -pin on" \
/usr/bin/time -v bash run.sh > run_100ns.log 2>&1
```

## 期待される観察 (論文 Fig 1-2 より)

| 項目 | 論文 insulin 結果 |
|---|---|
| % aggregated | "rapidly aggregated within 50 ns、 100% throughout simulation" |
| max aggregate size | hexamer (irreversible) |
| dynamic coexistence | **無** (octreotide と対照的) |

→ insulin は monotonic aggregation (degarelix と同様)。 我々の D-octreotide も
monotonic だったが、 insulin は天然で aggregation tendency が強い (β-sheet 形成、
fibril precursor)。

## Post-process + 解析

```bash
GMX=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256/gmx
python -m abmptools.trajectory thin_nojump --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10 --gmx $GMX
python -m abmptools.formulation analyze \
    --traj prod/prod.xtc --top prod/prod.gro --out analysis/ \
    --n-peptides 6 --enhancer-resnames "CPN,CPC" --bile-salt-resnames "TCH"
```

## 保管先

artifact: OneDrive `abmptools-dump/formulation-octreotide-pmf-smoke/aggregation-10mM/insulin/`
