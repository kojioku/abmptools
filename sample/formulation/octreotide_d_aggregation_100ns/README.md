# `octreotide_d_aggregation_100ns` — D-octreotide natural aggregation (Hossain Config A)

D-octreotide (RCSB 1SOC、 D-Phe1 + D-Trp4 + Thr-ol を含む実構造) × 2 +
Na-caprate × 32 (16 neutral + 16 charged) + taurocholate × 2 を **完全
random** に配置し、 **100 ns NPT で bias 一切なしで natural aggregation**
を観察。 `octreotide_l_aggregation_100ns/` の D 体ペア。

- peptide は whole-peptide GAFF + Gasteiger charges (acpype `-c gas`)
- 論文 (Hossain 2023) の Configuration A 準拠 (box 10 nm + random init +
  production 中 restraint なし)
- L 体 (ff14SB) と比べ D 体は water-peptide 親和性が GAFF Gasteiger で
  弱く、 production 中 NPT で box が大きく collapse する可能性あり (cluster
  smoke で踏んだ問題)。 100 ns の十分な長さがあれば water equilibrium 達成
  + collapse 緩和の見込み

## Build + run

```bash
cd sample/formulation/octreotide_d_aggregation_100ns
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/oct_d_agg100

cd /tmp/oct_d_agg100
GMX=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256/gmx \
PATH=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256:$PATH \
NT=8 MDRUN_OPTS="-nb gpu -pin on" \
/usr/bin/time -v bash run.sh > run_100ns.log 2>&1
```

## Post-process (基本セット)

prod.xtc を可視化・解析サイズに変換 (`abmptools.trajectory`、 Windows 互換):

```bash
cd /tmp/oct_d_agg100
python -m abmptools.trajectory thin_nojump \
    --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10
# → prod/prod_nojump_skip10.xtc (~300 MB、 100 frame、 1 ns stride、 -pbc nojump 済)

# VMD:
vmd prod/prod.tpr prod/prod_nojump_skip10.xtc
```

Python API:

```python
from abmptools.trajectory import thin_and_nojump
thin_and_nojump(trajectory="prod/prod.xtc", tpr="prod/prod.tpr", skip=10)
```

詳細は `sample/formulation/_postprocess/README.md` と
`python -m abmptools.trajectory --help` 参照。

## 期待 wall (RTX 4070 Ti)

L 体実機 (87,647 atoms, 9h56m) と同等の atom 数なら **~10 時間**。
D 体は acpype で生成された topology が ff14SB と微妙に異なる atom 数になる
可能性 (peptide 全体を 1 ligand 扱いだが water 充填量は box 体積で律速)。

### 実測 (2026-05-30、RTX 4070 Ti + GROMACS 2026.0 + 8 OpenMP thread)

| 項目 | 値 |
|---|---|
| atoms (実際の構成) | **83,337** |
| production performance | **246 ns/day** |
| **100 ns wall** | **9 時間 47 分** |
| メモリピーク | 286 MB |
| prod.xtc | 2.9 GB (1000 frame、100 ps stride) |

L 体 (87,647 atoms, 9h56m, 242 ns/day) とほぼ同等。

### aggregation 結果 (100 ns 平均)

| 観察項目 | D (GAFF Gasteiger) | L (ff14SB、 比較用) |
|---|---|---|
| Peptide↔Enhancer COM 距離 | 1.97 nm | 2.26 nm |
| Peptide↔Enhancer min heavy-atom 距離 | 0.31 nm (contact) | 0.29 nm |
| Peptide↔Enhancer heavy-atom contact pairs (<0.5 nm) | **452 (約 1.8×)** | 251 |
| Peptide↔BileSalt COM 距離 | 3.73 nm | 3.14 nm |
| Peptide↔BileSalt contact pairs | 22 | 46 |

→ **D 体は L 体より peptide-caprate aggregation が顕著に強い** (contact 数約 1.8×)。
D 体 GAFF + Gasteiger は water-peptide 親和性が低い分、 相対的に
peptide-enhancer 親和性が強く出る (cluster cohesion 強)。 これは 5/25 の
D 体 cluster smoke で box collapse が顕著だった事実とも整合 (compact
cluster + water 不足)。

artifact は OneDrive `abmptools-dump/formulation-octreotide-pmf-smoke/aggregation-100ns/`
に保管。L vs D 比較 plot (`oct_l_vs_d_agg100ns.png`) も同所。

## 比較対象

- `octreotide_l_aggregation_100ns/` — L 体、 ff14SB tleap
- 同 box 10 nm + 同 random init + 同 100 ns NPT で **力場差 (ff14SB vs
  whole-peptide GAFF Gasteiger)** が aggregation kinetics に与える影響を
  評価可能
