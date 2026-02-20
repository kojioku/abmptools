# Developer Quickstart

## Prerequisites

- Python 3.8+ (recommended)
- numpy
- pandas
- (Optional) gfortran — for Fortran-accelerated IFIE/PIEDA reading
- (Optional) OpenBabel — for XYZ→PDB conversion
- (Optional) UDFManager — for OCTA COGNAC MD workflows

## Installation from Source

```bash
git clone <repository-url>
cd abmptools

# Install in development/editable mode
pip install --user -e .

# Or standard install
pip install --user .
```

The install triggers `make` to compile the Fortran shared library. If gfortran is absent, this step fails silently and the package installs without Fortran acceleration.

## Building the Fortran Extension

```bash
# From repository root
make

# Verify the library was built
ls -la abmptools/f90/bin/readifiepiedalib.so
```

The `Makefile` compiles `abmptools/f90/src/readifiepiedalib.f90` into `abmptools/f90/bin/readifiepiedalib.so` using gfortran with `-shared -fPIC` flags.

## Running Sample Workflows

Each `sample/` subdirectory has a `run.sh` script:

```bash
# Generate an AJF input file
cd sample/generateajf && bash run.sh

# Convert a log to CPF
cd sample/log2cpf && bash run.sh

# Generate a DIFIE-averaged CPF
cd sample/generate_difie/TrpCage && bash run.sh

# Convert CPF versions
cd sample/convertcpf && bash run.sh
```

## Project Structure at a Glance

```
abmptools/
├── abmptools/          ← All source code here (27 modules)
│   ├── __init__.py     ← Package exports
│   ├── mol_io.py       ← START HERE: base I/O class
│   ├── cpfmanager.py   ← Central data manager (CPF files)
│   ├── getifiepieda.py ← Primary analysis CLI
│   └── ...
├── sample/             ← Working examples with run.sh
├── docs/               ← Documentation
└── setup.py            ← Build configuration
```

## Where to Start Reading

Recommended reading order for understanding the codebase:

1. **`abmptools/__init__.py`** — See what is exported and how classes are organized.
2. **`abmptools/mol_io.py`** — Base class for all molecular I/O. Understand coordinate handling.
3. **`abmptools/abinit_io.py`** → **`abmptools/pdb_io.py`** — Follow the inheritance chain upward.
4. **`abmptools/cpfmanager.py`** — Central data structure. Understand how CPF data maps to DataFrames.
5. **`abmptools/logmanager.py`** — How ABINIT-MP logs are parsed.
6. **`abmptools/getifiepieda.py`** — Main analysis entry point. See how all pieces compose.
7. **`abmptools/generateajf.py`** — FMO input generation pipeline.

## Code Conventions

Observations from the codebase (not formal style guide):

- **Class inheritance** is the primary composition mechanism (not mixins or protocols).
- **pandas DataFrames** are the universal internal data format for tabular data.
- **argparse** is used for all CLI entry points, invoked via `python -m abmptools.<module>`.
- **`multiprocessing.Pool`** is used for parallelism; worker functions are module-level for pickling.
- **Version detection** is common — modules auto-detect ABINIT-MP and CPF versions from file headers.
- **Unit conversions** use constants: `627.5095` (Hartree→kcal/mol), `0.529177` (Bohr→Å).
- Comments and variable names are a mix of English and Japanese.
- No formal test suite exists (noted as a future plan in `CHANGELOG.md`).

## Adding a New Module

To add a new analysis or conversion module:

1. Create `abmptools/new_module.py`.
2. Add CLI via `argparse` with an `if __name__ == '__main__':` block.
3. Import shared data structures from `cpfmanager` or `logmanager` as needed.
4. Inherit from `pdb_io` or `mol_io` if coordinate handling is needed.
5. Add a sample in `sample/new_module/` with input data and `run.sh`.
6. Update `__init__.py` if the module should be importable directly from `abmptools`.

Invoke as:
```bash
python -m abmptools.new_module [options]
```

## Common Development Tasks

### Parse a CPF file interactively

```python
import abmptools
cpf = abmptools.CPFManager()
cpf.parse("path/to/file.cpf")
print(cpf.atominfo.head())
print(cpf.diminfo.head())
print(cpf.static_data)
```

### Parse a LOG file interactively

```python
import abmptools
log = abmptools.LOGManager()
condition, fraginfo, static_data = log.parse("path/to/file.log")
print(condition)
print(fraginfo)
```

### Debug a CLI module

```bash
# Run with -h to see all options
python -m abmptools.getifiepieda -h

# Run with verbose output (most modules print progress to stdout)
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i test.log -nof90
```

## Known Limitations & Future Plans

From `CHANGELOG.md`:
- No formal test suite yet (planned).
- Refactoring is planned.
- Some local ABINIT-MP output variants may not be fully supported.
- BSSE reading functionality is under development.
