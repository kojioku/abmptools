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

### Release PMF (100 ns 最終 frame から release_us)

config の `release_us` section (n_windows=6, spacing=0.3 nm,
k=2000 kJ/mol/nm², window 200 ps, pull 400 ps) で 100 ns prod の最終
構造を入力に pull + 6 windows + WHAM を実行 (`/tmp/oct_l_agg100/us/run_us.sh`、
GPU wall 9 分):

| 指標 | L (ff14SB) | D (GAFF Gasteiger、 参考) |
|---|---|---|
| pull range | [1.65, 1.94] nm | [0.82, 1.13] nm |
| PMF min | -5.5 kJ/mol @ z=1.74 | -1.1 kJ/mol @ z=0.84 |
| PMF max (barrier) | +60.4 kJ/mol @ z=3.18 | +183.5 kJ/mol @ z=2.33 |
| **Release barrier** | **~66 kJ/mol** | **~185 kJ/mol** |

D 体の barrier が L 体の約 3 倍。 これは D 体 GAFF Gasteiger の cluster
cohesion が強い (contact pairs 452 vs L 251、 約 1.8×) ことと整合
(aggregation の cohesion が強い分、 release が熱力学的に難しい)。 ただし
D 体絶対値は GAFF Gasteiger force field の過剰評価の可能性あり、
Hossain CHARMM36 + CGenFF と直接比較は不可。

artifact + L vs D release PMF 比較 plot (`octreotide_release_pmf_l_vs_d.png`)
も OneDrive 保管。

### Release PMF (改善 sampling 版、 6 win × 5 ns × k=4000)

config の release_us を以下に拡張:
- `n_windows=6`、 `window_spacing_nm=0.2`
- `window_nsteps=2500000` (5 ns、 200 ps smoke の 25×)
- `force_constant_kj_mol_nm2=4000` (=2× smoke 2000、 thermal 抑制で box-half 制約回避)
- `pull_nsteps=200000` (400 ps、 smoke と同じ)

| 指標 | L (ff14SB) | D (GAFF Gasteiger) |
|---|---|---|
| wall (RTX 4070 Ti) | 2h45m | 2h40m |
| z range | [1.48, 2.93] | [0.64, 1.98] |
| **Cluster-bound z (E=0 基準)** | **z=1.48** | **z=0.64** |
| **Release barrier peak** | **~97 kJ/mol @ z=2.6 nm** | **~90 kJ/mol @ z=1.85 nm** |

200 ps smoke (L:60 / D:185 kJ/mol) と比べ:
- D 体は **185 → 90 kJ/mol に下落** (sampling 厚で overestimate 解消)
- L 体は **60 → 97 kJ/mol に上昇** (sampling 厚で barrier 顕在化)
- 両系で release barrier ~90 kJ/mol 程度に収束、 smoke で見えた 1.8× の
  差は消失

ただし z 開始位置は L vs D で大きく異なる (L: 1.48 vs D: 0.64 nm)、
これは 100 ns aggregation での cluster cohesion 差 (L=surface 接触型、
D=核埋め込み型) を反映。 release barrier 自体は両系で comparable に。

artifact + plot (`octreotide_release_pmf_l_vs_d_5ns.png`) も OneDrive
保管。

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
