# `octreotide_d_cluster_smoke` — D-octreotide cluster-templated smoke

D-octreotide (RCSB 1SOC、 D-Phe1 + D-Trp4 + Thr-ol) を caprate shell の
**中心に最初から置く** 構成で packmol packing。peptide 全体は GAFF2 +
Gasteiger charges で parametrize (AM1-BCC SCF が 137-atom peptide で
発散するため `-c gas` にフォールバック、 `octreotide_d_gaff_smoke/`
と同じ手法)。

## ポイント

- `system.packmol_layout = "cluster_core"`
- `peptides[0].parameterize_method = "gaff"` + `pdb_path` で 1SOC PDB を流用
- `box_size_nm = 12.0` (D 体 GAFF mode は production 中 NPT で box が
  collapse しやすい — L 体は 10 nm で十分だが D 体は 12 nm でも安全側)
- `release_us.n_windows = 4` (D 体 GAFF mode の cluster cohesion が弱く
  pull range が広い、 box-half rule を確実に避けるため少なめ)

詳細は `octreotide_l_cluster_smoke/README.md` および
`octreotide_d_gaff_smoke/README.md` を参照。

## Run

```bash
cd sample/formulation/octreotide_d_cluster_smoke
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/oct_d_cluster
micromamba run -n abmptoolsenv python -m abmptools.formulation release_us \
    --config config.json --output-dir /tmp/oct_d_cluster
cd /tmp/oct_d_cluster
NT=4 bash run.sh
cd us && NT=4 bash run_us.sh
```

## Smoke 限界

D 体 GAFF + Gasteiger mode は water-peptide interaction が弱く、
production 200 ps 中に box が 12 → 7-8 nm 程度に collapse することが
あります。 そうなると pull-coord 3D COM 距離が grompp の box-half
rule (= 0.49 × min box edge) に触れて release_us が fail することが
あります。 production scale (50+ ns) で実用的な PMF を取るには:

1. production を NPT → NVT に切り替え (box 固定、 ただし density 反映無し)
2. または box を更に大きく (15+ nm) して NPT 維持
3. または cluster radii を tight に縮める (core 1.5 / shell 2.5)

を検討してください。
