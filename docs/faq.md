# FAQ

## How do I install ABMPTools?

```bash
pip install --user .
```

This installs the `abmptools` package and attempts to compile the Fortran shared library. If `gfortran` is not available, the install still succeeds — the Fortran acceleration is optional.

Ensure `numpy` and `pandas` are installed in your environment (they are not auto-installed; see [dependencies.md](dependencies.md)).

## What Python version is required?

Python 3.6+ is assumed based on code features. Python 3.8+ is recommended. No explicit version constraint is declared in `setup.py`.

## How do I run a specific tool?

All tools are invoked as Python modules:

```bash
python -m abmptools.<module_name> [options]
```

Examples:
```bash
python -m abmptools.generateajf -i protein.pdb --method MP2 -basis 6-31G*
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i calc.log
python -m abmptools.log2cpf -i calc.log -o output.cpf
python -m abmptools.convertcpf -i large.cpf -f 1-100 -v 23
```

Use `-h` with any module to see available options:
```bash
python -m abmptools.getifiepieda -h
```

## What is CPF and which versions are supported?

CPF (Coordinate Property File) is ABINIT-MP's output format containing atomic coordinates, fragment definitions, monomer/dimer energies, and PIEDA components. Supported versions:

| Version | Notes |
|---------|-------|
| 4.201 | Legacy format |
| 7.0 | MIZUHO variant |
| 10 | Default since v1.8.0 |
| 23 | Latest (Open 1.0 Rev.23); used for DIFIE output |

`CPFManager.parse()` auto-detects the version. Use `convertcpf` to convert between versions:
```bash
python -m abmptools.convertcpf -i input.cpf -v 23
```

## What is DIFIE?

DIFIE (Dynamic IFIE) is a time-averaged CPF file generated from multiple MD trajectory snapshots. It computes the mean and standard deviation of IFIE/PIEDA values across snapshots and stores them in a single CPF file with `M-` (mean) and `S-` (std) prefixed column headers.

```bash
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 10 1 -f 1-100 -np 4
```

The output CPF can be visualized in Biostation Viewer.

## Do I need gfortran?

No. gfortran is optional and only needed to compile the Fortran shared library (`readifiepiedalib.so`) that accelerates IFIE/PIEDA reading from log files.

- **With gfortran**: Run `make` or install via `pip install --user .` to compile the library.
- **Without gfortran**: All functionality works via pure Python. Use the `-nof90` flag with `getifiepieda` to explicitly use the Python path.

Note that some features are only available via specific paths:
- MP3/MP4 extraction requires the Fortran module.
- PB-IFIE, BSSE-IFIE, monomer/dimer energies require the pure Python mode (`-nof90`).

## How do I use the Fortran-accelerated IFIE reader?

If the Fortran library is compiled (present at `abmptools/f90/bin/readifiepiedalib.so`), it is loaded automatically by `getifiepieda.py`. To disable it:

```bash
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i calc.log -nof90
```

To compile manually:
```bash
make    # in repository root
```

## What are the supported ABINIT-MP versions?

| ABINIT-MP Version | Revisions | AJF Version Flag |
|--------------------|-----------|------------------|
| v1 | Rev.10–23 | `rev22`, `rev23`, etc. |
| v2 | Rev.4–8 | `v2rev4`, `v2rev8` |

Specify the version when generating AJF files:
```bash
python -m abmptools.generateajf -i protein.pdb -ajfv v2rev8
```

## How do I handle nucleic acid systems?

Nucleic acid support was added in v1.14.6. The `LOGManager` (`logmanager.py`) and `log2config` module handle:
- Pure nucleic acid systems.
- Protein-nucleic acid complexes.
- CYS disulfide bridge handling in mixed systems.

Use the standard `log2config` or `log2cpf` workflows — nucleic acid detection is automatic.

## Where are sample workflows?

The `sample/` directory contains working examples:

| Directory | Workflow | Run |
|-----------|----------|-----|
| `sample/convertcpf/` | CPF conversion/filtering | `bash run.sh` |
| `sample/generate_difie/` | DIFIE averaging (TrpCage, CS4) | `bash run.sh` |
| `sample/gbsa/` | GBSA solvation setup | `bash run.sh` |
| `sample/generateajf/` | AJF generation | `bash run.sh` |
| `sample/log2config/` | Log to config extraction | `bash run.sh` |
| `sample/log2cpf/` | Log to CPF conversion | `bash run.sh` |
| `sample/rmsd/` | RMSD analysis | `bash run.sh` |

## How do I use CPFManager programmatically?

```python
import abmptools

cpf = abmptools.CPFManager()
cpf.parse("input.cpf")

# Access data as pandas DataFrames
print(cpf.atominfo)     # Per-atom data
print(cpf.fraginfo)     # Fragment definitions
print(cpf.diminfo)      # Dimer interaction data (IFIE/PIEDA)
print(cpf.mominfo)      # Monomer energies
print(cpf.static_data)  # Summary statistics

# Write to new CPF
cpf.write("Title", "output.cpf")
```

## How do I speed up multi-sample analysis?

Use the `-np` flag to enable parallel processing:

```bash
python -m abmptools.getifiepieda --multi 10 -d 8.0 -t 1 100 1 -i template.log -np 8
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 50 1 -np 8
```

This uses Python's `multiprocessing.Pool` to read and process log/CPF files in parallel. See [parallelization.md](parallelization.md) for details.

## What is the `tips/` directory?

The `tips/` directory contains auxiliary scripts for MD post-processing workflows with AMBER, GROMACS, and NAMD. These are standalone helper scripts, not part of the core `abmptools` package.
