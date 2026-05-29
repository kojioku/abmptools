# `octreotide_l_aggregation_100ns` — Natural aggregation run (Hossain Configuration A)

L-octreotide (FCFWKTCT + Cys2-Cys7 disulfide + ACE/NMe cap) × 2 + Na-caprate ×
32 (16 neutral + 16 charged) + taurocholate × 2 を **完全 random** に配置し、
**100 ns NPT で bias 一切なしで natural aggregation** を観察する系。
Hossain 2023 の **Configuration A** (Methods page 3 を参照) に準拠。

- 論文と box size 一致 (10 nm cubic)
- 論文と initial layout 一致 (random uniform)
- 論文と production 流儀一致 (NPT、 restraint なし)
- 違いは production 時間のみ (論文 500 ns vs 100 ns)

## Build + run

build は abmptoolsenv (tleap + acpype + packmol)、 mdrun は gmxcudaenv
(GROMACS 2026 GPU 対応) で行う:

```bash
# Build (~30 秒)
cd sample/formulation/octreotide_l_aggregation_100ns
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/oct_l_agg100

# Run (em + nvt + npt + 100 ns prod) — GPU 利用、 wall 測定込み
cd /tmp/oct_l_agg100
# em + nvt + npt は CPU で軽い、 production だけ GPU 指定
NT=8 MDRUN_OPTS="-nb gpu -pin on" /usr/bin/time -v bash run.sh
```

## 期待 wall time

| Atom 数 | system size | GPU (RTX 4070 Ti) | CPU 4-core |
|---|---|---|---|
| ~47,000 | 10 nm cubic | ~8-12 時間 / 100 ns | ~3-7 日 / 100 ns |

`gmx mdrun` の log に `Performance: X ns/day` が出るので、 実際の throughput
で 100 ns wall を予測できる。

### 実測 (2026-05-27/28、RTX 4070 Ti + GROMACS 2026.0 + 8 OpenMP thread)

| 項目 | 値 |
|---|---|
| atoms (実際の構成) | **87,647** (box 10 nm cubic、TIP3P 充填) |
| production performance | **242 ns/day** |
| **100 ns wall** | **9 時間 56 分** |
| メモリピーク | 296 MB |
| prod.xtc | 3.0 GB (1000 frame、100 ps stride) |

GROMACS 2026 は `-nb gpu` 指定だけで PME / update / bonded も自動で
GPU offload する (auto-detection)。 mdrun log に
`PME tasks will do all aspects on the GPU` /
`PP task will update and constrain coordinates on the GPU` と出る。

### aggregation 結果 (100 ns 平均)

| 観察項目 | Peptide↔Enhancer | Peptide↔BileSalt |
|---|---|---|
| Mean COM 距離 | 2.26 nm | 3.14 nm |
| Mean min heavy-atom 距離 | 0.29 nm (contact) | 0.35 nm |
| Mean heavy-atom contact pairs (<0.5 nm) | 251 | 46 |

→ peptide-enhancer aggregate は 100 ns 全期間維持 (Hossain の octreotide
"dynamic coexistence" と整合)、 bile salt は 2 mol で間欠的接触。
artifact は OneDrive `abmptools-dump/formulation-octreotide-pmf-smoke/aggregation-100ns/`
に保管。

## 観察項目

- 完了後 `prod/prod.xtc` を MDAnalysis で時系列解析:
  - peptide-enhancer COM 距離の時間変化
  - heavy-atom contact count
  - aggregation onset time (peptide-enhancer min distance が <0.5 nm に下がる時刻)
- Hossain Figure 1b の "% aggregated molecules vs time" と比較

## Configuration B (cluster-init) との違い

このサンプルは **配置 bias 一切なし** で自然 aggregation を観察。
`octreotide_l_cluster_smoke/` (Configuration B) は peptide を中心に
pre-form して production restraint で保持する shortcut で、 PMF を smoke
で取得するための便宜上の手段。
