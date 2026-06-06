# `octreotide_d_aggregation_10mM` — D 体真の論文準拠 (Hossain 2023 主要系)

D-octreotide (RCSB 1SOC + ACE/NMe cap、 D-Phe1 + D-Trp4 + Thr-ol 含む実 octreotide)
× **6** + Na-caprate × 32 + taurocholate × 2 を **完全 random** に配置し、
**100 ns NPT** で natural aggregation を観察する系。

既存 `octreotide_d_aggregation_100ns/` (peptide × 2 = 3.3 mM hybrid 構成) の
**peptide × 6 = 10 mM 版**。 L 体 `octreotide_l_aggregation_10mM/` と並列。

## 既存 100ns 系との違い

| 項目 | 既存 (100ns/) | 本系 (10mM/) |
|---|---|---|
| D-octreotide n_copies | 2 (= 3.3 mM) | **6 (= 10 mM)** |
| 力場 | whole-peptide GAFF + Gasteiger | 同 |
| caprate (16+16=32)、 TCH × 2、 box 10 nm | 同 | 同 |
| atom 数 (推定) | ~83k | **~76k** (peptide 増分で water 減) |
| Fig 1-2 transition network 比較 | 不可 (max dimer) | **可能** (hexamer まで) |

## Build + run

```bash
cd sample/formulation/octreotide_d_aggregation_10mM
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/oct_d_10mM

cd /tmp/oct_d_10mM
PATH=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256:$PATH \
NT=8 MDRUN_OPTS="-nb gpu -pin on" \
/usr/bin/time -v bash run.sh > run_100ns.log 2>&1
```

## 期待 wall (RTX 4070 Ti、 GROMACS 2026)

| 系 | atoms | performance | wall |
|---|---|---|---|
| L 体 10mM (実測) | 81,944 | 243 ns/day | 9h52m |
| D 体 10mM (予測) | ~76,000 | ~250 ns/day | **~9h30m** |

## Post-process (基本セット)

```bash
cd /tmp/oct_d_10mM
GMX=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256/gmx
python -m abmptools.trajectory thin_nojump \
    --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10 --gmx $GMX
```

## 期待される観察項目 (論文 Fig 対応)

| 解析 | 論文 Fig | L 体 10mM (参考) |
|---|---|---|
| Aggregate size transition (hexamer 観察) | Fig 1c / Fig 2b | hexamer 51.5%、 trimer 19.8% |
| Pep↔Enh contacts (per peptide) | Fig 4 | **120.7** (Pep↔Pep の 3.1×) |
| Secondary structure (β-turn 等) | Fig 5 | turn 9.4% + 3₁₀-helix 4.4% + β-bridge 4.7% |

D 体は L 体と比較して mirror image + chiral discrimination + β-turn 安定性差
(memory `project_abmptools_formulation.md` 参照) で違いが期待される。

## 保管先

- 入力: 本 dir (`config.json` + `README.md`)
- artifact: OneDrive `abmptools-dump/formulation-octreotide-pmf-smoke/aggregation-10mM/octreotide_d/`
  (build 完了後に backup)
