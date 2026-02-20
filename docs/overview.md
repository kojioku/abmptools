# ABMPTools Overview

## Project Name & Version

**ABMPTools** (ABINIT-MP Tools) — v1.14.6
Author: Koji Okuwaki

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
| MD Integration | `udf_io`, `udfrm_io`, `udfcreate` | OCTA COGNAC UDF file handling and conversion |

## Supported ABINIT-MP Versions

- ABINIT-MP v1: Rev.10–23
- ABINIT-MP v2: Rev.4–8
- CPF format versions: 4.201, 7.0 (MIZUHO), 10, 23

## Installation

```bash
pip install --user .
```

This runs `make` during install to compile the optional Fortran shared library (`readifiepiedalib.so`). If `gfortran` is not available, the install still succeeds — the Fortran acceleration is optional.

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
