# Architecture

## Design Overview

ABMPTools follows an **inheritance-based layered architecture**:

- A **core I/O layer** provides coordinate reading, molecular math, and fragment management.
- **Data managers** (`CPFManager`, `LOGManager`) handle ABINIT-MP-specific file formats independently.
- **Analysis and setup modules** compose the above layers into CLI-driven workflows.

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

### `anlfmo.py` — Advanced FMO Analysis

Inherits from `pdb_io`. Provides:
- `readlog()` — Comprehensive log parsing with fragment interaction extraction.
- Time-dependent IFIE calculation.
- Multi-sample analysis with parallel processing.
- Column definitions for CSV output formatting.

### `generate_difie.py` — Dynamic IFIE

Creates time-averaged CPF files from MD trajectory snapshots:
- `getcpfobj((time, template, padding, fragments))` — load the CPF at one
  timestep. Takes a single tuple so it can be handed to `Pool.map`.
- `getavestddf(cpfs, distdf, atomdf, charge_labels)` — mean/std across the
  trajectory.
- Outputs `<input stem>-DIFIE.cpf` with `M-` (mean) and `S-` (std) prefixed
  columns per pair.
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

## Shared Layer (`abmptools/core/`)

Small pieces that more than one sub-package needs, kept here so none of them
has to import another.

| module | 提供するもの | 利用側 |
|---|---|---|
| `system_model.py` | `SystemModel` / `SimulationParams` — MD 系の記述と実行条件 | `amorphous` |
| `acpype.py` | acpype (GAFF2 / AM1-BCC) の subprocess ラッパ。`run_acpype` / `AcpypeResult` / `LigandParameterization` | 中立層 (out-of-tree consumers) |
| `_subprocess.py` | 外部コマンド実行の薄いラッパ (`CommandError` / `run_command`) | `core.acpype` |

`acpype.py` と `_subprocess.py` は v2.9.0 で `core` に移した。以前はそれぞれ別の
サブパッケージの内部にあり、そこへ手を伸ばす形の相互依存が残っていた。

## IFIE/PIEDA Reading

`anlfmo.read_ifiepieda()` parses the log in Python. The column set follows the
log's `Method` keyword — `MP2`, `MP3` and `CCPT` (MP4) each name their columns
differently, and MP4 carries one extra (`GRIMME-MP4`); see `_IFIE_COLUMNS`.

A Fortran shared library used to offer an alternative reader. It was removed
once the Python reader became the faster of the two, and it never supported
`GRIMME-MP4`, PB-IFIE, BSSE-IFIE or the monomer/dimer energies.

**Note**: MP3/CCPT logs are only handled by the `--ffmatrix` mode; the other
modes still assume the MP2 column names.

## Where to Start Reading

1. **`mol_io.py`** → `abinit_io.py` → `pdb_io.py` — Follow the inheritance chain to understand the I/O foundation.
2. **`cpfmanager.py`** — Understand how CPF data is parsed into DataFrames; this is the central data structure.
3. **`getifiepieda.py`** — See how analysis CLIs compose the base classes and data managers.
4. **`generateajf.py`** — Understand the FMO input generation pipeline.
5. **`sample/*/run.sh`** — Run working examples to see inputs and outputs.
