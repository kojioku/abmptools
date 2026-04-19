# Data Flow

## Overview

ABMPTools processes molecular data through three main pipelines:

1. **FMO Setup** — Structural data → AJF input files for ABINIT-MP.
2. **Analysis** — ABINIT-MP outputs (LOG/CPF) → CSV results with IFIE/PIEDA data.
3. **DIFIE Averaging** — Multiple CPF snapshots → Single averaged CPF.

## Main Data Flow Diagram

```
                        ┌─────────────┐
                        │  PDB file   │
                        └──────┬──────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
        ┌──────────┐    ┌──────────┐    ┌──────────────┐
        │ pdb2fmo  │    │generateajf│   │  pdbmodify   │
        │ (config) │    │          │    │  (editing)   │
        └────┬─────┘    └────┬─────┘    └──────┬───────┘
             │               │                 │
             ▼               ▼                 ▼
        ┌──────────────────────┐          ┌─────────┐
        │    AJF input file    │          │ PDB out │
        └──────────┬───────────┘          └─────────┘
                   │
                   ▼
          ┌────────────────┐
          │   ABINIT-MP    │  (external calculation)
          │   execution    │
          └───────┬────────┘
                  │
         ┌────────┴────────┐
         ▼                 ▼
   ┌──────────┐      ┌──────────┐
   │ LOG file │      │ CPF file │
   └────┬─────┘      └────┬─────┘
        │                  │
   ┌────┴────┐        ┌────┴──────────┐
   ▼         ▼        ▼               ▼
┌───────┐ ┌───────┐ ┌───────────┐ ┌────────────┐
│log2cpf│ │getifie│ │convertcpf │ │generate    │
│       │ │pieda  │ │           │ │  _difie    │
└───┬───┘ └───┬───┘ └─────┬─────┘ └─────┬──────┘
    ▼         ▼           ▼              ▼
  CPF      CSV files    CPF (filtered)  DIFIE CPF
           (IFIE/PIEDA   (version       (averaged)
            matrices)     converted)
```

## FMO Setup Pipeline (PDB → AJF)

### Protein Systems

```
PDB ──→ pdb_io.readpdb() ──→ internal atom/residue data
                                      │
                                      ▼
                              generateajf.py
                              (method, basis, solvation, etc.)
                                      │
                                      ▼
                                  AJF output
```

Modules involved: `pdb_io.py` → `generateajf.py`

### Non-Protein Systems (Molecular Aggregates)

```
                ┌───────────────┐
                │ 1-molecule AJF│  (manual fragment definition)
                └───────┬───────┘
                        │
                        ▼
                   ajf2config.py ──→ segment_data.dat (config dict)
                                           │
PDB (full system) ──→ pdb2fmo.py ◀─────────┘
                          │
                          ▼
                   setfmo.setrfmoparam()
                   (cut mode: sphere/cube/around/neutral/none)
                          │
                          ▼
                    AJF + PDB output
```

Modules involved: `ajf2config.py` → `pdb2fmo.py` → `setfmo.py` → `generateajf.py`

### UDF (OCTA COGNAC MD) Systems

```
UDF file ──→ udfrm_io.convert_udf_pdb() ──→ PDB (per-frame)
                                                   │
                                                   ▼
                                          udf2fmo.py / pdb2fmo.py
                                                   │
                                                   ▼
                                              AJF output
```

Modules involved: `udf_io.py` → `udfrm_io.py` → `udf2fmo.py`

## Analysis Pipeline (LOG/CPF → CSV)

### Single-Point IFIE Analysis

```
LOG file ──→ anlfmo.readlog()
                 │
                 ├──→ version detection (getversion)
                 ├──→ fragment info (getfraginfo)
                 ├──→ IFIE/PIEDA data extraction
                 │         │
                 │    [optional: Fortran SO for speed]
                 │         │
                 ▼         ▼
           getifiepieda.py (mode selection)
                 │
         ┌───────┼────────┬──────────┐
         ▼       ▼        ▼          ▼
      frag     mol      ffmatrix   fragids
      mode     mode     mode       mode
         │       │        │          │
         ▼       ▼        ▼          ▼
       CSV     CSV    CSV files    CSV
    (per-frag) (per-mol) (per-PIEDA  (per-frag)
                         component)
```

### Time-Series / Multi-Sample Analysis

```
LOG files (numbered) ──→ getifiepieda --multi / --tfmatrix
                              │
                         multiprocessing.Pool
                         (parallel log reading)
                              │
                              ▼
                         Aggregated CSV
                    (IFIE over time, sum, individual)
```

Key output files per mode:
- `--frag`: `log-fragi-fragj-ifie.csv`
- `--ffmatrix`: `{ES,EX,CT,DI,HF,MP2,...}-ffmatrix.csv`, `Distance-ffmatrix.csv`
- `--tfmatrix`: Time × fragment matrix CSV (one file per reference fragment)
- `--multi`: `fragi-distr-ifiedt.csv` (per-pair), `fragi-distr-ifiesum.csv` (totals)

### CPF-Based Analysis

```
CPF file ──→ CPFManager.parse()
                  │
                  ├──→ atominfo DataFrame
                  ├──→ fraginfo DataFrame
                  ├──→ diminfo DataFrame (IFIE/PIEDA components)
                  ├──→ mominfo DataFrame (monomer energies)
                  └──→ static_data dict
                          │
                          ▼
                    cpf2ifielist.py
                    (filter by fragment range)
                          │
                          ▼
                   Formatted CSV (kcal/mol, Å)
```

## DIFIE Pipeline (Multiple CPFs → Averaged CPF)

```
CPF files (trajectory snapshots)
    file-001.cpf, file-002.cpf, ..., file-N.cpf
                        │
                        ▼
              generate_difie.py
              -t start end interval
              -np (parallel reading)
                        │
              ┌─────────┴─────────┐
              ▼                   ▼
    getcpfobj() per timestep   setoutcpfstatic()
    (parallel via Pool)        (prepare output structure)
              │                   │
              ▼                   ▼
         getavestddf()      representative
         (mean, std)        structure coords
              │                   │
              └─────────┬─────────┘
                        ▼
                  DIFIE CPF output
            (M-/S- prefixed columns for
             mean and std of each property)
```

## File Conversion Paths

```
LOG ──→ log2cpf.py    ──→ CPF
LOG ──→ log2config.py ──→ Python config dict
AJF ──→ ajf2config.py ──→ segment_data.dat
CIF ──→ readcif.py    ──→ Cartesian coords (XYZ/PDB)
CPF ──→ convertcpf.py ──→ CPF (different version / fragment subset)
UDF ──→ udfrm_io.py   ──→ PDB / XYZ
PDB ──→ pdbmodify.py  ──→ PDB (edited: renamed, moved, re-numbered)
AJF ──→ ajfserial.py  ──→ Numbered AJF series (for trajectory FMO)
```

## Amorphous Build Pipeline (SMILES/SDF → GROMACS MD)

```
SMILES / SDF  ──→ molecule_prep.py     ──→ OpenFF Molecule + single-mol PDB
                      │
                      ↓
                 density.py             ──→ count & cubic-box size from target density
                      │
                      ↓
                 packing.py (packmol)   ──→ mixture.pdb (multi-component packed)
                      │
                      ↓
                 parameterizer.py       ──→ system.gro / system.top (OpenFF + AM1-BCC)
                      │                    via openff-interchange
                      ↓
                 ndx_writer.py          ──→ system.ndx  (per-component index groups)
                      │
                      ↓
                 mdp_protocol.py        ──→ 01_em.mdp / 02_nvt / 03_npt / 04_anneal /
                      │                    05_npt_final.mdp  + run_all.sh + wrap_pbc.sh
                      ↓
                 (run_all.sh)           ──→ 5-stage GROMACS MD
                      ↓
                 (wrap_pbc.sh)          ──→ *_pbc.xtc / 05_npt_final_pbc.gro (VMD ready)
```

Orchestrated end-to-end by `abmptools.amorphous.builder.AmorphousBuilder.build()`.

## Internal Data Structures

All modules converge on **pandas DataFrames** for structured data:

| DataFrame | Source Module | Key Columns |
|-----------|--------------|-------------|
| `atominfo` | `cpfmanager` | alabels, elems, elemtypes, resnames, resids, fragids, xcoords, ycoords, zcoords, MUL-HF |
| `fraginfo` | `cpfmanager` | Fragment definitions, atom ranges |
| `diminfo` | `cpfmanager` | fragi, fragj, min-dist, NR, HF, ES, EX, CT, DQ |
| `mominfo` | `cpfmanager` | fragi, DPM-HF-{X,Y,Z}, NR, HF |
| `static_data` | `cpfmanager`, `logmanager` | nuclear_repulsion_energy, total_energy, natom, nfrag, ndimer |

Unit conventions:
- Energies in IFIE/PIEDA tables: **kcal/mol** (conversion factor: 627.5095 from Hartree).
- Distances: **Angstrom** (conversion factor: 0.529177 from Bohr).
