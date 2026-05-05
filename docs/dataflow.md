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
PubChem CID / name ──→ amorphous.pubchem ──→ 3D SDF (local cache)
                                  │
                                  ↓   (or direct SMILES / SDF input)
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

## CG Builder Pipeline (Martini 3, 1.18.0+)

Two CG sub-packages live under `abmptools.cg/`. Both treat external
tools (`vermouth`, `insane`, `gmx`, `tleap`) as subprocess-only black
boxes — no source bundled or modified.

### `abmptools.cg.peptide` — Martini 3 peptide CG box (1.18.0)

```
PeptideBuildConfig (JSON)
    │
    ↓
peptide_atomistic.py     ──→ atomistic peptide PDB
   (tleap or extended-backbone fallback)
    │
    ↓
martinize_runner.py      ──→ <name>_cg.pdb + <name>.itp
   (martinize2 -ff martini3001, vermouth Apache-2.0 subprocess)
    │
    ↓
system_packer.py         ──→ packed.gro
   (gmx insert-molecules)
    │
    ↓
top_writer.py            ──→ topol.top  (M3 ITPs + peptide ITP)
    │
    ↓
water_box.py             ──→ martini_v3.0.0_water.gro (auto-generated
   if not in ff/, gmx insert-molecules with W density 8.36/nm³)
    │
    ↓
system_packer.solvate    ──→ system_solv.gro  (gmx solvate -cs)
    │
    ↓
system_packer.add_ions   ──→ system_ions.gro  (gmx grompp + genion,
                              NaCl 0.15 M, neutralize)
    │
    ↓
mdp_templates.py + run script writer
                         ──→ mdp/{em,nvt,npt,md}.mdp + index.ndx + run.sh
```

Orchestrated end-to-end by `abmptools.cg.peptide.builder.PeptideCGBuilder.build()`.

### `abmptools.cg.membrane` — Martini 3 peptide-membrane PMF (1.19.0)

```
MembraneCGBuildConfig (JSON)
    │
    ↓
_copy_ff_files            ──→ output_dir/martini_v3.0.0*.itp (4 files)
    │
    ↓
PeptideCGBuilder (sub-call,    ──→ molecules/<name>/{name}_cg.pdb + .itp
solvent_enabled=False,
mdp_*=False)
    │
    ↓
insane_runner.run_insane  ──→ bilayer.gro + insane_topol.top
   (insane GPL-2.0 subprocess; peptide -dm <z_offset> in POPC + W + NaCl)
    │
    ↓
topology_composer.compose_topology
                          ──→ topol.top  (4 M3 ITPs + peptide ITP +
                              "Protein" → moleculetype name +
                              NA+/CL- → NA/CL normalisation)
    │
    ↓
topology_composer.normalize_ion_atom_names_gro
                          ──→ system_ions.gro  (NA+/CL- → NA/CL in .gro)
    │
    ↓
system_packer.write_ndx_from_gro_cg
                          ──→ index.ndx  (Bilayer / Peptide / W / NA / CL /
                                         Non_Bilayer; no gmx make_ndx)
    │
    ↓
pulling.find_pbc_center_atom (imported from abmptools.membrane)
                          ──→ bilayer pbc-atom index
    │
    ↓
mdp_templates + pulling.write_pulling_mdp_cg + umbrella.write_window_mdps
                          ──→ mdp/{em,nvt,npt}.mdp +
                              pull/pull.mdp (NVT-chassis + direction-periodic) +
                              windows/win_NNN/window.mdp × N
                              (NPT-semiisotropic + direction + pbcatom)
    │
    ↓
umbrella.write_run_script ──→ run.sh
    │
    ↓
(bash run.sh)   em → nvt → npt → pull →
                python -m abmptools.cg.membrane make-windows →
                per-window grompp + mdrun (× N) →
                python -m abmptools.cg.membrane wham
                          ──→ analysis/pmf.xvg + histo.xvg
```

Orchestrated end-to-end by `abmptools.cg.membrane.builder.MembraneCGBuilder.build()` (build phase only) + `run.sh` (MD execution + post-analysis).

`pulling.parse_pullx_xvg / extract_window_frames / find_pbc_center_atom`,
`mdp_us_protocol.render_pull_block`, and `pmf.run_wham` are imported from
`abmptools.membrane.*` (duck-typed via `UmbrellaCGProtocol` field-name
match with `USProtocol`) — no helper duplication between CG and AA.

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
