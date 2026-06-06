# `octreotide_l_aggregation_10mM` — 真の論文準拠 (Hossain 2023 Table 1 主要系)

L-octreotide (FCFWKTCT + Cys2-Cys7 disulfide + ACE/NMe cap) **× 6** + Na-caprate
× 32 (16 neutral + 16 charged) + taurocholate × 2 を **完全 random** に配置し、
**100 ns NPT で bias 一切なし** で natural aggregation を観察する系。
Hossain 2023 の **主要系 (10 mM peptide、Fig 1-2 等の出典)** に **mol 数として準拠**。

## 既存 `octreotide_l_aggregation_100ns/` との違い

| 項目 | 既存 (100ns/) | **本系 (10mM/)** |
|---|---|---|
| octreotide n_copies | 2 (= 3.3 mM) | **6 (= 10 mM)** ← 修正点 |
| caprate (16+16=32) | 同 | 同 |
| taurocholate × 2 | 同 | 同 |
| box 10 nm cubic | 同 | 同 |
| atom 数 (推定) | ~87k | **~99k** |
| 論文との対応 | hybrid (caprate/tauro は主要系準拠だが peptide 量が baseline 系) | **主要系 (Fig 1-2, 4-7 出典) と一致** |
| Fig 1-2 (trimer/hexamer transition) との比較 | **原理的に不可** (peptide ≤ 2) | **可能** (peptide × 6 で trimer/hexamer 観測可) |

## なぜ peptide × 6 か

論文 p3:
> systems containing 10 mM peptides, PEs, and taurocholate

10 nm cubic box = 10⁻²¹ L、 10 mM = 6.022×10²¹ /L
→ 10 mM × 10⁻²¹ L × NA = **6.022 ≈ 6 molecule per box**

論文 p9 「3.3 mM octreotide alone」は **peptide-only baseline** (caprate/tauro なし)
の比較系で、 mixed colloid 系の出典ではない。

## Build + run

```bash
# Build (~30-60 秒、 packmol で peptide × 6 が入る box 探索に時間)
cd sample/formulation/octreotide_l_aggregation_10mM
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/oct_l_10mM

# Run (em + nvt + npt + 100 ns prod) — GPU 利用
cd /tmp/oct_l_10mM
NT=8 MDRUN_OPTS="-nb gpu -pin on" /usr/bin/time -v bash run.sh > run_100ns.log 2>&1
```

## 期待 wall time (RTX 4070 Ti + GROMACS 2026)

| 構成 | atom 数 | performance | 100 ns wall |
|---|---|---|---|
| 既存 100ns 系 | 87,647 | 242 ns/day | 9h56m |
| **本系 10mM** | ~99,000 (推定) | ~210 ns/day (atom 数増分で減速) | **~11.5 時間** |

## Post-process (基本セット)

```bash
cd /tmp/oct_l_10mM
python -m abmptools.trajectory thin_nojump \
    --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10
# → prod/prod_nojump_skip10.xtc (100 frame、 1 ns stride)
```

## 期待される観察項目 (論文 Fig 対応)

| 解析 | 論文 Fig | 我々の data から取得 |
|---|---|---|
| Aggregate size transition network | Fig 2b | `abmptools.formulation.analysis.aggregate_transition` |
| % aggregated vs time | Fig 1b | 同 (peptide × 6 で 0/16/33/50/66/83/100% の 7 状態) |
| Maximum aggregate size vs time | Fig 1c | 同 |
| Peptide residue-residue contacts / peptide | Fig 4 | `analysis.contact_map` (final 100 ns 平均) |
| Secondary structure (β-turn 維持率) | Fig 5 | `analysis.secondary_structure` (`gmx dssp`) |
| Peptide SASA | p11 | `analysis.sasa` (`gmx sasa`) |
| Caprate-peptide H-bond | Fig 7 | `analysis.hbond` (`gmx hbond`) |

論文 Fig 2b (octreotide + 50mM C10 + 3mM TC random init) の transition network
における dominant aggregate size = **trimer**、 hexamer まで観察される予測。

## Release PMF (任意)

`release_us` section は既存 100ns sample と同じ設定 (6 win × 5 ns × k=4000)。
peptide × 6 から 1 peptide を引きはがす形なので、 PMF 値域は既存と異なる可能性
(より大きい cluster からの release ≈ より高い barrier、 論文 octreotide で
~30-50 kJ/mol の range)。

## 保管先

- 入力: 本 dir (`config.json` + `README.md`)
- artifact: OneDrive `abmptools-dump/formulation-octreotide-pmf-smoke/aggregation-10mM/octreotide_l/`
  (build 完了後に backup)
