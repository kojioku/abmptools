# ABMPTools Overview

## Project Name & Version

**ABMPTools** (ABINIT-MP Tools)
Author: Koji Okuwaki

> Latest released version on PyPI: v1.15.4 (refactor/all line).
> v1.16.0 (core / amorphous updates), v1.17.0 (membrane subpackage),
> v1.17.1 (mixed-lipid support), v1.17.2 (extended lipid table + discovery
> API), v1.17.3 (CHARMM36 backend 実機検証 + 7 件 bug fix),
> **v1.18.0 (`cg.peptide` subpackage — Martini 3 peptide CG builder)**,
> **v1.19.0 (`cg.membrane` subpackage — Martini 3 peptide-membrane
> PMF via insane + umbrella sampling)**,
> **v1.20.0 (`genesis.grest` subpackage — GENESIS gREST_SSCR replica
> exchange with solute tempering)**,
> **v1.21.0 (`fragmenter` subpackage — FMO automatic fragment splitter
> for small molecules / lipids / polymers)**,
> **v1.22.0 (`genesis.mmgbsa` subpackage — GENESIS MM/GBSA
> single-point ΔG_bind via atdyn implicit GBSA)**,
> and **v1.23.0 (`crystal` subpackage — organic-crystal FMO pipeline:
> CIF → supercell → fragment cut → ABINIT-MP AJF + HPC jobscripts;
> Phase A–D = skeleton + regression fixture + 8-subcommand CLI +
> `--run-local` + ASE PBC unwrap, includes 4 public-molecule MP2/6-31G(d)
> reference set: urea / glycine / benzene / naphthalene)** are staged on
> `develop`, not yet released to `main` / PyPI. See
> [`CHANGELOG.md`](../CHANGELOG.md) for the per-version detail.

## What is ABMPTools?

ABMPTools is a Python toolkit for pre-processing, post-processing, and analysis of Fragment Molecular Orbital (FMO) calculations performed with [ABINIT-MP](https://fmodd.jp/member_contents/manual_ABINIT-MP/). It provides:

- **IFIE/PIEDA analysis** — Extract and visualize inter-fragment interaction energies from calculation logs and CPF files.
- **CPF parsing & conversion** — Read, write, filter, and version-convert CPF (Coordinate Property File) binary/text files.
- **FMO input generation** — Automatically generate AJF input files from PDB structures with fragment assignment.
- **File format conversion** — Convert between PDB, LOG, CPF, UDF (OCTA COGNAC), CIF, and XYZ formats.
- **Dynamic IFIE (DIFIE)** — Average IFIE data across MD trajectory snapshots into a single CPF with mean/standard-deviation statistics.

## Key Capabilities

| Category | Modules | Description |
|----------|---------|-------------|
| Analysis | `getifiepieda`, `anlfmo`, `cpf2ifielist`, `getcharge` | IFIE/PIEDA extraction, fragment interaction matrices, charge analysis |
| CPF Management | `cpfmanager`, `convertcpf`, `generate_difie`, `log2cpf` | Parse/write/convert CPF files, create DIFIE averages |
| FMO Setup | `generateajf`, `pdb2fmo`, `udf2fmo`, `setfmo`, `addsolvfrag` | Generate AJF inputs, assign fragments, handle solvation |
| File Conversion | `log2config`, `ajf2config`, `readcif`, `pdbmodify`, `ajfserial` | Format conversion and PDB editing utilities |
| MD Integration | `udf_io`, `udfrm_io`, `udfcreate`, `gro2udf`, `udf2gro` | OCTA COGNAC UDF file handling and GROMACS ↔ OCTA conversion |
| Structure Optimization | `geomopt.{MacePdbOptimizer, OpenFFOpenMMMinimizer, QMOptimizerPySCF}` | MACE / OpenFF / PySCF-DFT driven geometry optimization for PDB inputs |
| Amorphous Builder | `amorphous` (`build_amorphous.py`) | SMILES / SDF / PubChem CID (via `amorphous.pubchem`) multi-component amorphous builder (Packmol + OpenFF + AM1-BCC), auto-generates 5-stage GROMACS annealing protocol and VMD-friendly trajectory post-processing |
| Membrane Builder | `membrane` (`MembraneUSBuilder`) | Lipid-bilayer + peptide umbrella-sampling PMF builder (packmol-memgen + AMBER `ff19SB`/`Lipid21`/TIP3P or CHARMM36 backend → semiisotropic NPT → z-pulling → per-window US → `gmx wham`). GPU-aware `run.sh`. Designed commercial-license-clean (no CGenFF / CHARMM-GUI). |
| Crystal-FMO Pipeline | `crystal` (`CrystalOrchestrator`, `abmp-crystal` CLI) | Organic-crystal FMO workflow: CIF → supercell PDB (legacy parser or ASE backend with PBC unwrap) → fragment cut around solute → ABINIT-MP AJF (full-precision `&XYZ` block) → PJM/SLURM/PBS/local jobscripts → optional `--run-local` execution → `getifiepieda` postprocessing. 8-subcommand CLI driven by single YAML/JSON config. |
| CG Peptide Builder | `cg.peptide` (`PeptideBuildConfig`, vermouth-martinize2) | Martini 3 peptide CG system builder. tleap (optional) で AA PDB → martinize2 で CG 化 → GROMACS solvate / genion / MDP / run.sh まで end-to-end (v1.18.0+). |
| CG Membrane PMF | `cg.membrane` (`MembraneCGBuildConfig`, insane + cg.peptide sub-call) | Martini 3 peptide-bilayer umbrella-sampling PMF builder. AA membrane より 30-100× 速い (v1.19.0+). |
| GENESIS gREST | `genesis.grest` (`GrestBuilder`, `python -m abmptools.genesis.grest` CLI) | GENESIS gREST_SSCR (Replica-Exchange Solute Tempering with Solute Side Chain Repartitioning) 5-stage builder + replica transition / acceptance / PMF 解析 (AMBER ff19SB + TIP3P; v1.20.0+). |
| GENESIS MM/GBSA | `genesis.mmgbsa` (`MmgbsaBuilder`, `python -m abmptools.genesis.mmgbsa` CLI) | GENESIS atdyn `[ENERGY] implicit_solvent=GBSA` による protein-ligand 単フレーム MM/GBSA ΔG_bind builder (AMBER ff14SB + GAFF/GAFF2 via acpype; v1.22.0+). |
| Auto Fragmentation | `fragmenter` (`AutoFragmenter`, `python -m abmptools.fragmenter` CLI) | FMO 自動フラグメント分割: PDB → 連結成分 → canonical SMILES グループ化 → C-C 単結合切断 → log2config 互換 segment_data.dat 出力。Jupyter UI (ipywidgets) と headless CLI (SVG + JSON 編集) の 2 系統 (v1.21.0+). |

## Supported ABINIT-MP Versions

- ABINIT-MP v1: Rev.10–23
- ABINIT-MP v2: Rev.4–8
- CPF format versions: 4.201, 7.0 (MIZUHO), 10, 23

## Installation

Editable install is recommended for day-to-day use and development:

```bash
pip install -e .
```

Non-editable install (e.g. for production deployment):

```bash
pip install .
```

`--user` is usually unnecessary; pip handles both virtual environments and system Python appropriately. Installation runs `make` to compile the optional Fortran shared library (`readifiepiedalib.so`). If `gfortran` is not available, the install still succeeds — the Fortran acceleration is optional.

### Requirements

- **Required**: numpy, pandas
- **Optional**: UDFManager (OCTA COGNAC), gfortran (Fortran acceleration), OpenBabel (`obabel` CLI)

## Quick Example

```bash
# Generate an AJF input from a PDB file
python -m abmptools.generateajf -i protein.pdb -basis 6-31G* --method MP2

# Extract IFIE for fragment 10, within 8 Å distance
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i calculation.log

# Convert a log file to CPF format
python -m abmptools.log2cpf -i calculation.log -o output.cpf

# Create a DIFIE-averaged CPF from trajectory snapshots
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 10 1 -f 1-100 -np 4
```

## Where to Start Reading

For new developers approaching this codebase:

1. **`abmptools/__init__.py`** — Package entry point; see what is exported.
2. **`abmptools/mol_io.py`** — Base I/O class; all coordinate handling builds on this.
3. **`abmptools/cpfmanager.py`** — Central data manager; understand CPF parsing and pandas DataFrame usage.
4. **`abmptools/getifiepieda.py`** — Primary analysis CLI; shows how modules compose together.
5. **`docs/ABMPTools-user-manual.md`** — Comprehensive user manual with option descriptions and output examples.
6. **`sample/`** — Working examples with `run.sh` scripts for each workflow.
7. **`docs/tutorial_membrane_us.md`** — Step-by-step ops tutorial for the
   `membrane` subpackage (peptide-bilayer Umbrella Sampling, AA backend).
8. **`docs/cg_membrane.md`** + **`docs/tutorial_cg_membrane_us.md`** —
   Martini 3 CG version of the membrane PMF builder
   (`cg.membrane`); 30-100× faster wall time than AA, smoke run completes
   in 5-6 min on 8-thread CPU.
9. **`docs/cg_peptide.md`** — Martini 3 peptide CG builder
   (`cg.peptide`); peptide-only system in water box, used standalone or
   sub-called by `cg.membrane`.
10. **`docs/peptide_builders.md`** — Cross-cutting selection guide for the
    three peptide-from-sequence builders (AA membrane / CG peptide /
    CG membrane); compares resolution, cost, license, and use cases.
11. **`docs/fragmenter.md`** — `abmptools.fragmenter` (1.21.0+); FMO
    automatic fragment splitter for small molecules / lipids / polymers
    (canonical SMILES grouping + C-C MW walk + Jupyter UI / headless CLI).
12. **`docs/grest.md`** + **`docs/tutorial_grest.md`** —
    `abmptools.genesis.grest` (1.20.0+); GENESIS gREST_SSCR replica
    exchange with solute tempering (AMBER ff19SB + TIP3P, 4-12 replicas
    + REST residue selection + temperature ladder).
13. **`docs/mmgbsa.md`** + **`docs/tutorial_mmgbsa.md`** —
    `abmptools.genesis.mmgbsa` (1.22.0+); GENESIS atdyn-based MM/GBSA
    single-point ΔG_bind for protein-ligand complexes (AMBER ff14SB +
    GAFF/GAFF2 via acpype, 4-stage pipeline: split PDB → parameterize →
    GBSA single-point → ΔG_bind aggregation).
    Protein / DNA stays on the existing `log2config` route.
14. **`docs/crystal.md`** + **`docs/tutorial_crystal_fmo.md`** —
    `abmptools.crystal` (1.23.0+); organic-crystal FMO pipeline
    (CIF → supercell → fragment cut → ABINIT-MP AJF + HPC jobscripts;
    8-subcommand CLI: `abmp-crystal {expand,fragment,jobs,pipeline,
    postproc,nearest,validate,example}`). Companion docs:
    **`docs/crystal_verification.md`** (verification matrix) +
    **`docs/crystal_public_molecule_references.md`** (4-molecule
    MP2/6-31G(d) reference summary: urea / glycine / benzene /
    naphthalene with E' in crystal vs MP2 total isolated, sum-IFIE,
    wall time, and CIF source attribution).
