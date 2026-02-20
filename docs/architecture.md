# Architecture

## Design Overview

ABMPTools follows an **inheritance-based layered architecture**:

- A **core I/O layer** provides coordinate reading, molecular math, and fragment management.
- **Data managers** (`CPFManager`, `LOGManager`) handle ABINIT-MP-specific file formats independently.
- **Analysis and setup modules** compose the above layers into CLI-driven workflows.
- A **Fortran extension** accelerates the performance-critical IFIE/PIEDA reading path.

Key design patterns:
- **Class inheritance chain** for progressive capability building.
- **Pandas DataFrames** as the universal internal data representation.
- **argparse-based CLIs** invoked via `python -m abmptools.<module>`.
- **`multiprocessing.Pool`** for embarrassingly parallel workloads.

## Class Inheritance Hierarchy

```
molcalc                          (coordinate math, distance, rotation, PBC)
  └── udf_io                    (OCTA UDF file reading)
        └── udfrm_io            (UDF → PDB/XYZ conversion)

mol_io                           (base molecular I/O — XYZ, coordinate conversion)
  └── abinit_io                 (ABINIT-MP I/O — fragment definitions)
        └── pdb_io              (PDB parsing — residues, atoms, read modes)
              └── anlfmo        (advanced FMO analysis — log parsing, IFIE)

setfmo(pdb_io, udfcreate, udfrm_io)   (FMO setup — multiple inheritance)

CPFManager                       (standalone — CPF parse/write)
LOGManager                       (standalone — ABINIT-MP log parse)
```

**Assumption**: The inheritance hierarchy is inferred from class definitions in each module. `setfmo` uses multiple inheritance combining PDB, UDF creation, and UDF conversion capabilities.

## Core I/O Layer

### `mol_io.py` — Molecular I/O Base

The foundation class. Provides:
- `read_xyz()` — Read XYZ coordinate files.
- `getatoms()` — Extract atom lists from structures.
- `convert_xyzs_pdb()` — XYZ to PDB format conversion.

### `molcalc.py` — Coordinate Mathematics

Standalone utility class for:
- Distance calculation (`getdist()`), center of mass (`getCenter()`).
- Molecular translation (`moveMolTrans()`) and rotation (`moveMolEuler()`).
- Periodic boundary condition handling.
- LAMMPS data/trajectory parsing (`parse_lammps_data()`, `parse_lammps_trajectory()`).
- Uses `multiprocessing.Pool` for distance matrix calculations.

### `abinit_io.py` — ABINIT-MP I/O Base

Extends `mol_io` with:
- Fragment definition management.
- Fragment geometry extraction.
- Support for various ABINIT-MP calculation modes and options.
- Parallel processing for fragment setup.

### `pdb_io.py` — PDB Parser

Extends `abinit_io` with:
- `readpdb()` — Parse PDB files with modes: `TER`, `resnum`, `rfile`.
- Residue and atom property extraction.
- Extended atom numbering support (5+ digit atoms).
- Optional `UDFManager` integration for OCTA MD data.

## Data Managers

### `CPFManager` (`cpfmanager.py`)

Central data hub for CPF files. Key methods:
- `parse(filename)` — Read CPF files (auto-detects version).
- `read_header()`, `read_atominfo()`, `read_fraginfo()`, `read_dimer()` — Section parsers.
- `write(title, filename)` — Generate CPF output.

Stores data in pandas DataFrames:
- `atominfo` — Per-atom properties (coords, charges, residue info).
- `fraginfo` — Fragment-level data.
- `diminfo` — Dimer interaction data (distance, energy components).
- `mominfo` — Monomer/dipole moment information.

Supported CPF versions: 4.201, 7.0, 10, 23 (including `.gz` compressed files).

### `LOGManager` (`logmanager.py`)

Parses ABINIT-MP calculation log files. Key methods:
- `parse(filename)` — Full log parsing.
- `getversion()` — Detect ABINIT-MP version from log header.
- `getcondition()` — Extract calculation conditions.
- `getfraginfo()` — Extract fragment definitions.

Returns: `condition` dict, `fraginfo` dict, `static_data` dict.
Supports nucleic acids and protein-nucleic acid complexes (added v1.14.6).

## Analysis Modules

### `getifiepieda.py` — IFIE/PIEDA Extraction

The primary analysis CLI. Supports multiple modes:

| Mode | Flag | Description |
|------|------|-------------|
| Fragment pair | `--frag i j` | IFIE between specified fragments |
| Distance filter | `--frag i -d r` | Fragments within distance `r` of fragment `i` |
| Molecule-based | `--mol i -d r` | Interactions from molecule `i` within distance `r` |
| Fragment matrix | `--ffmatrix i1-i2 j1-j2` | Full i×j IFIE matrix |
| Time-fragment matrix | `--tfmatrix i-j k-l` | Time-series fragment interaction matrix |
| Multi-sample | `--multi i -t start end interval` | Time-series IFIE with parallel reading |
| In-molecule | `--fraginmol i j MOL k` | Intra-molecular fragment interactions |

Uses the Fortran library (`readifiepiedalib.so`) when available; falls back to Python with `-nof90`.

### `anlfmo.py` — Advanced FMO Analysis

Inherits from `pdb_io`. Provides:
- `readlog()` — Comprehensive log parsing with fragment interaction extraction.
- Time-dependent IFIE calculation.
- Multi-sample analysis with parallel processing.
- Column definitions for CSV output formatting.

### `generate_difie.py` — Dynamic IFIE

Creates time-averaged CPF files from MD trajectory snapshots:
- `getcpfobj()` — Load CPF at specific timestep.
- `getavestddf()` — Compute mean/std across trajectory.
- Outputs a single CPF with `M-` (mean) and `S-` (std) prefixed columns.
- Parallel processing via `-np` flag.

## FMO Setup Pipeline

### `generateajf.py` — AJF Template Generation

Full-featured CLI for ABINIT-MP input file creation:
- Input: PDB file or config dictionary.
- Options: method (HF/MP2), basis set, solvation (PB), PIEDA, CPF output, RESP, DGEMM, dispersion, BSSE.
- Version-aware output for different ABINIT-MP revisions.

### `setfmo.py` — FMO Setup Orchestrator

Multiple-inheritance class combining PDB, UDF, and configuration capabilities:
- `setrfmoparam()` — Configure FMO parameters.
- Cut modes: `sphere`, `cube`, `around`, `neutral`, `none`.
- Solute/solvent specification, ion handling (`remain`/`remove`).
- Fragment assignment for molecular aggregates.

### `pdb2fmo.py` — PDB to FMO Converter

CLI for non-protein systems (molecular aggregates, polymers):
- Reads config file with fragment definitions for one molecule.
- Assigns fragments to all molecules in the system.
- Produces AJF + PDB output set.

## File Conversion Modules

```
LOG  ──→ log2cpf   ──→ CPF
LOG  ──→ log2config ──→ Config dict
AJF  ──→ ajf2config ──→ segment_data.dat
CIF  ──→ readcif    ──→ Cartesian coordinates
CPF  ──→ convertcpf ──→ CPF (different version / filtered)
UDF  ──→ udfrm_io   ──→ PDB / XYZ
AJF  ──→ ajfserial  ──→ Numbered AJF files
PDB  ──→ pdbmodify  ──→ PDB (edited)
```

## Fortran Extension Integration

`abmptools/f90/src/readifiepiedalib.f90` (219 lines) provides a fast log file parser for IFIE/PIEDA data. It is:

1. Compiled to `readifiepiedalib.so` via `Makefile` (using gfortran with `-shared -fPIC`).
2. Loaded at runtime by `getifiepieda.py` using `ctypes`.
3. Optional — the `-nof90` flag forces pure-Python fallback.

**Note**: MP3/MP4 extraction requires the Fortran module. PB-IFIE, BSSE-IFIE, and monomer/dimer energies require the `-nof90` (pure Python) path.

## CLI Entry Points

All CLIs are invoked as `python -m abmptools.<module>`:

| Module | Purpose | Key Arguments |
|--------|---------|---------------|
| `generateajf` | Generate AJF templates | `-i` (PDB), `-basis`, `--method`, `-pb`, `-np` |
| `getifiepieda` | Extract IFIE/PIEDA | `--frag`, `--mol`, `--ffmatrix`, `--tfmatrix`, `-t`, `-np` |
| `log2cpf` | Convert log to CPF | `-i` (input), `-o` (output) |
| `log2config` | Convert log to config | `-i` (input), `-o` (output), `-np` |
| `pdb2fmo` | PDB to FMO setup | `-i` (PDB), `-p` (parameter file) |
| `generate_difie` | Create DIFIE CPF | `-i` (template), `-t` (range), `-z` (padding), `-np` |
| `convertcpf` | Convert CPF versions | `-i` (input), `-v` (version), `-f` (fragments) |
| `udf2fmo` | UDF to FMO setup | `-i` (UDF), `-p` (param), `-s` (solutes), `-r` (record) |
| `pdbmodify` | Modify PDB files | `-i`, `-move`, `-mode`, `-str` |
| `getcharge` | Extract charges | `-i` (log), `-t` (type), `-f` (fragments) |
| `addsolvfrag` | Add solvation | `-i` (PDB), `-temp` (template), `-solv` |
| `ajfserial` | Numbered AJF files | `-i` (template), `-t` (range), `-str` |
| `cpf2ifielist` | CPF to IFIE list | `-i` (CPF), `-f` (fragment range) |

## Where to Start Reading

1. **`mol_io.py`** → `abinit_io.py` → `pdb_io.py` — Follow the inheritance chain to understand the I/O foundation.
2. **`cpfmanager.py`** — Understand how CPF data is parsed into DataFrames; this is the central data structure.
3. **`getifiepieda.py`** — See how analysis CLIs compose the base classes and data managers.
4. **`generateajf.py`** — Understand the FMO input generation pipeline.
5. **`sample/*/run.sh`** — Run working examples to see inputs and outputs.
