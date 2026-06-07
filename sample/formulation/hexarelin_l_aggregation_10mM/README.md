# `hexarelin_l_aggregation_10mM` — 論文 hexarelin の L 体近似 (Phase 2-B sequence build)

論文 Hossain 2023 の **hexarelin** (His-D-2-Me-Trp-Ala-Trp-D-Phe-Lys-NH₂、 6 残基、
D-amino + methylated Trp + NH₂ cap) を **natural L-AA で近似** した系。

近似:
- D-2-Me-Trp → Trp (W) — methylation を無視、 L 体
- D-Phe → Phe (F) — L 体
- 結果 sequence: `HWAWFK` (6 残基、 すべて natural L-AA)
- NH₂ cap は OXT のみ (Phase 2-B sequence build は cap 未サポート、 PeptideBuilder
  の bare peptide PDB を使用)

10 mM 濃度の論文準拠系 (peptide × 6 in 10 nm cubic box) + caprate × 32 + TC × 2。

## 論文との制約

| 項目 | 論文 hexarelin | 我々 (L 体近似) |
|---|---|---|
| 残基 | His-D-2-Me-Trp-Ala-Trp-D-Phe-Lys | His-Trp-Ala-Trp-Phe-Lys |
| chirality | D-amino 含む | 全 L (近似) |
| methylation | D-2-Me-Trp 有 | 無 (W) |
| cap | NH₂ | OXT (PeptideBuilder default) |
| force field | CHARMM36 + CGenFF | ff14SB + GAFF2 (商用 OK) |

論文との PMF / aggregation の **絶対値は乖離**する可能性があるが、 **質的傾向 (mixed
colloid 形成、 Trp/Lys 主要 contact 等) は再現できる予定**。 真の hexarelin (D 体含む) は
whole-peptide GAFF route で別途実装可。

## Build + run

```bash
cd sample/formulation/hexarelin_l_aggregation_10mM
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/hexarelin_l_10mM

cd /tmp/hexarelin_l_10mM
PATH=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256:$PATH \
NT=8 MDRUN_OPTS="-nb gpu -pin on" \
/usr/bin/time -v bash run.sh > run_100ns.log 2>&1
```

## 期待値 (octreotide L 体 10mM の実測から外挿)

| 項目 | hexarelin L 近似 (予測) |
|---|---|
| atom 数 (推定) | ~78,000 |
| wall (RTX 4070 Ti) | ~9-10 時間 |
| 期待挙動 | mixed colloid 形成、 Trp/Lys が主要 contact |

## Post-process + 解析

```bash
cd /tmp/hexarelin_l_10mM
GMX=/home/okuwaki/.local/share/mamba/envs/gmxcudaenv/bin.AVX2_256/gmx
python -m abmptools.trajectory thin_nojump --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10 --gmx $GMX

python -m abmptools.formulation analyze \
    --traj prod/prod.xtc --top prod/prod.gro --out analysis/ \
    --n-peptides 6 --enhancer-resnames "CPN,CPC" --bile-salt-resnames "TCH"
```

## 保管先

artifact: OneDrive `abmptools-dump/formulation-octreotide-pmf-smoke/aggregation-10mM/hexarelin_l/`
