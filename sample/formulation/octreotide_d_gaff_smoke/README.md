# `octreotide_d_gaff_smoke` — D-octreotide via whole-peptide GAFF

Hossain 2023 system 4 で扱われている **本物の D 体 octreotide**
(D-Phe1-Cys-Phe-D-Trp-Lys-Thr-Cys-Thr-ol) を扱う smoke build。

D-amino acid + non-standard residues (Thr-ol C-term alcohol) は
tleap (ff14SB) の standard amino library で組めないため、**peptide
全体を 1 つの ligand として GAFF2 + AM1-BCC で parameterise**する。

これは `parameterize_method = "gaff"` に設定することで自動的に
small-molecule pipeline (acpype) に流される機能を使う。

## Source PDB

| File | 出典 |
|---|---|
| `octreotide_1SOC.pdb` | RCSB PDB ID **1SOC** (Sandostatin NMR structure, Pohl et al. 1995). license-free per RCSB terms. |

## Run

```bash
cd sample/formulation/octreotide_d_gaff_smoke
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/octreotide_d_gaff_smoke
```

## Force-field rationale

D-Phe1 + D-Trp4 + L-Cys2 + L-Phe3 + L-Lys5 + L-Thr6 + L-Cys7 +
L-Thr-ol8 with a disulfide bridge between Cys2 and Cys7.

By parameterising the entire peptide as one GAFF2 small molecule
via acpype:
- D stereochemistry preserved (acpype uses 3D coords as given)
- Disulfide bridge detected via bond perception in antechamber
- Thr-ol C-term alcohol is just another aliphatic OH

**Trade-off**: bonded parameters per atom type from a generic GAFF2
database; sidechain dihedrals are less accurate than ff14SB's
residue-specific terms. For peptide aggregation (Hossain 2023 question)
this is a fair approximation; for receptor binding affinity it would
not suffice.
