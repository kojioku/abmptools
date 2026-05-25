# `octreotide_l_cluster_smoke` — Cluster-templated initial state

L-octreotide × 2 を caprate shell の **中心に最初から置く** 構成で
packmol packing。Production 200 ps でも peptide-cluster aggregate
が initial で確保され、release-PMF が物理的に意味のある形で取れる
想定の smoke。

## ポイント

- `system.packmol_layout = "cluster_core"`
- `cluster_core_radius_nm = 1.2` (peptide 2 copies が入る中心球、 半径 1.2 nm)
- `cluster_shell_radius_nm = 3.0` (caprate × 32 が入る shell の外側 半径 3.0 nm)
- bile salt は shell の外、 box 全体に random
- production 200 ps → cluster-bound state を維持しつつ side chain relax
- release_us pull-group1 = Enhancer (caprate cluster COM) vs Peptide

## Run

```bash
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config sample/formulation/octreotide_l_cluster_smoke/config.json \
    --output-dir /tmp/oct_l_cluster

micromamba run -n abmptoolsenv python -m abmptools.formulation release_us \
    --config sample/formulation/octreotide_l_cluster_smoke/config.json \
    --output-dir /tmp/oct_l_cluster

cd /tmp/oct_l_cluster
NT=4 bash run.sh        # em + nvt + npt + production 200 ps
cd us
NT=4 bash run_us.sh     # pull + 6 windows × 200 ps + WHAM
```
