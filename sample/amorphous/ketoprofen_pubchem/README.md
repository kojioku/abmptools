# ketoprofen (PubChem 3D SDF input)

Amorphous ketoprofen build using the **PubChem 3D SDF** (CID 3825) as
the initial conformer, instead of `--smiles` which asks OpenFF to
generate a single conformer internally. Useful when the conformer
matters (reactive sites, hindered rotation, stereo-sensitive packings).

## Quick start

```bash
cd sample/amorphous/ketoprofen_pubchem
bash run_sample.sh

# After the build completes:
cd md
bash run_all.sh     # 5-stage annealing MD (~10 min on CPU 8 cores)
bash wrap_pbc.sh    # produces *_pbc.xtc for VMD
```

The SDF at `input/ketoprofen_pubchem_cid3825.sdf` is bundled so the
sample runs offline. If you delete it, `run_sample.sh` will
re-download it from the PubChem REST API on the next run.

## System

| 項目 | 値 |
|---|---|
| 分子 | Ketoprofen (C16H14O3, MW 254.28) |
| 分子数 | 50 |
| 原子数 | 1650 (33 atoms/mol, 19 heavy + 14 explicit H) |
| 目標密度 | 0.8 g/cm³ |
| ボックス | 2.977 nm (自動算出) |
| 力場 | OpenFF `openff_unconstrained-2.1.0.offxml` |
| 電荷 | AM1-BCC (AmberTools) |
| 乱数シード | 42 |
| 立体 | PubChem 3D は片側エナンチオマー ((S) or (R)) を返す |

## vs. SMILES-based sibling

[`../ketoprofen/`](../ketoprofen/) runs the same protocol starting from
a SMILES string. In a 1.3 ns annealing test both approaches converge
to nearly identical box / density:

| 量 | SMILES 版 | PubChem 3D 版 |
|---|---|---|
| Final Box | 2.639 nm | 2.647 nm |
| Final Density | 1138.8 kg/m³ | 1140.8 kg/m³ |
| Temperature | 299.7 K | 299.85 K |

See the private archive entry `md-archive/ketoprofen_pubchem_amorphous_20260419_n50_d0.8/README.md` for the full comparison.

## Input SDF source

Retrieved from NCBI PubChem:

```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/3825/SDF?record_type=3d
```

PubChem data is produced by NIH/NCBI and is in the public domain
(freely redistributable). If this sample is used in a derivative
publication, please cite PubChem:

> Kim S, Chen J, Cheng T, Gindulyte A, He J, He S, Li Q, Shoemaker BA,
> Thiessen PA, Yu B, Zaslavsky L, Zhang J, Bolton EE.
> **PubChem 2025 update.** Nucleic Acids Res. 2025;53(D1):D1516-D1525.
> doi:[10.1093/nar/gkae1059](https://doi.org/10.1093/nar/gkae1059)

Specific compound: PubChem CID 3825 (ketoprofen), 3D conformer computed
by PubChem with MMFF94.

## Note

The 3D SDF encodes one stereoisomer (chirality is fixed in the file).
For an explicitly racemic amorphous, mix (R)- and (S)-ketoprofen SDFs
as two components via `--mol ketoprofen_R.sdf ketoprofen_S.sdf` with
matching `--n_mol` counts.
