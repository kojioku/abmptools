# ABMPTools (ABINIT-MP Tools)

A Python toolkit for pre-processing, post-processing, and analysis of Fragment Molecular Orbital (FMO) calculations with [ABINIT-MP](https://fmodd.jp/member_contents/manual_ABINIT-MP/).

## Features

### IFIE/PIEDA Analysis (`getifiepieda`, `anlfmo`, `cpf2ifielist`, `getcharge`)

- Distance-filtered IFIE tables for target fragments or molecules
- Fragment–fragment interaction matrices (1:1, 1:N, N:1, N:N)
- Time-series IFIE from MD-FMO trajectory snapshots
- SVD-based interaction decomposition
- Charge extraction from ABINIT-MP logs

### CPF Management (`cpfmanager`, `convertcpf`, `generate_difie`, `log2cpf`)

- Parse and write CPF files (versions 4.201, 7.0 MIZUHO, 10, 23)
- Version conversion between CPF formats
- Residue-based CPF extraction
- Dynamic IFIE (DIFIE) averaging across MD snapshots with mean/σ statistics
- Generate CPF from ABINIT-MP log files

### FMO Input Generation (`generateajf`, `pdb2fmo`, `udf2fmo`, `setfmo`, `addsolvfrag`)

- Auto-generate AJF input files from PDB structures
- Fragment assignment for proteins and molecular assemblies
- Solvation fragment addition
- Support for sp2 fragmentation and various basis sets

### File Format Conversion

- CIF → PDB/XYZ (`readcif`) with symmetry operations
- ABINIT-MP log → fragment config (`log2config`, `ajf2config`)
- PDB editing and serial AJF generation (`pdbmodify`, `ajfserial`)

### GROMACS ↔ OCTA COGNAC Conversion

- **udf2gro**: Convert OCTA UDF files to GROMACS format (`.gro`, `.top`, `.mdp`, `.itp`)
- **gro2udf**: Convert GROMACS files to OCTA UDF format (supports `--from-top` mode)
- **udfcharge**: Transfer per-atom partial charges from a single-molecule UDF to bulk molecules (`transfer`), or restore a neutralized UDF's charges to a target integer formal charge (`restore`)

### Geometry Optimization (`geomopt`)

- **MacePdbOptimizer**: MACE/ASE-based PDB structure optimization
- **OpenFFOpenMMMinimizer**: OpenFF force-field minimization via OpenMM
- **QMOptimizerPySCF**: Quantum chemistry optimization with PySCF

### Amorphous Structure Building (`amorphous`, `build_amorphous.py`)

- Multi-component amorphous system construction (API + polymer / API + API / binary mixture)
- Initial structures from either **SMILES** (OpenFF conformer generation), external **3D SDF/MOL files** (`--mol`), or **PubChem CID / name** (`--pubchem_cid` / `--pubchem_name`, auto-downloads MMFF94 3D SDF; raises `PubChemNo3DError` when no 3D conformer exists)
- Packmol-based packing + OpenFF force field parameterization + AM1-BCC charges
- Auto-generates GROMACS inputs and a 5-stage annealing protocol
  (EM → high-T NVT → high-T NPT → simulated annealing → low-T NPT equilibration)
- Bundled `md/run_all.sh` drives the MD run; `md/wrap_pbc.sh` post-processes
  trajectories with `gmx trjconv -pbc mol -ur compact` for VMD-friendly
  `*_pbc.xtc` / `_pbc.gro` outputs

## Supported ABINIT-MP Versions

- ABINIT-MP v1: Rev.10–23
- ABINIT-MP v2: Rev.4–8

## Installation

### From PyPI (normal use)

```bash
pip install abmptools
```

Updating:

```bash
pip install --upgrade abmptools
```

Checking what you actually got — run this **outside the repository**, since a
source tree in the current directory shadows the installed package:

```bash
cd /tmp
pip show abmptools | head -2
python -c "import abmptools; print(abmptools.__file__)"
```

If `Location` and `__file__` disagree, another copy is winning.

Optional extras pull in the plotting and chemistry dependencies a given
subpackage needs, e.g. `pip install 'abmptools[amorphous]'`.

### From source (development)

Editable install is recommended for day-to-day use and development:

```bash
pip install -e .
```

Non-editable install (e.g. for production deployment):

```bash
pip install .
```

`--user` is usually unnecessary; pip handles both virtual environments and system Python appropriately.

Installation runs `make` to compile the optional Fortran shared library for accelerated IFIE/PIEDA reading. If `gfortran` is not available, the install still succeeds without Fortran acceleration.

### Requirements

- **Required**: Python 3.8+, numpy, pandas
- **Optional**: UDFManager (OCTA COGNAC), gfortran, OpenBabel, PySCF, ASE, OpenMM, Packmol

## Testing

```bash
pytest tests/ -v                     # 927 tests collected (2.12.0 時点)
pytest tests/ -v -k molcalc          # specific module
pytest tests/test_regression.py -v   # regression tests (60 bundled + 16 gated)
```

See [tests/TEST_COVERAGE.md](tests/TEST_COVERAGE.md) for details.

### Regression Tests

`tests/test_regression.py` compares current CLI output against reference
fixtures stored in `tests/regression/reference/` (generated from the
pre-refactor state). This guards against behavior drift during refactoring.

Covered tools: `generateajf`, `log2cpf`, `convertcpf`, `udf2gro`, `gro2udf`,
and `getifiepieda`.

**Developer-only tests**: the 16 `getifiepieda` regression cases require
external sample data (the internal `abmptools-sample` repository) at:

```
../abmptools-sample/sample/getifiepieda/
├── 6lu7-multi-fmolog/    (extracted from abmptools-fmolog-sample.tar.bz2)
├── cd7-fmolog/
├── 6m0j-pb-fmolog/
└── xyzfile/
```

These tests are automatically skipped when the data is not available, so
public CI runs are unaffected.

## Samples

Each `sample/` subdirectory contains input data and a `run.sh` / `run_sample.sh` script:

```bash
# FMO / IFIE / CPF samples
cd sample/generateajf            && bash run.sh
cd sample/log2cpf                && bash run.sh
cd sample/generate_difie/TrpCage && bash run.sh
cd sample/convertcpf             && bash run.sh

# Amorphous structure builder samples — see sample/amorphous/README.md for the full index
cd sample/amorphous/pentane_benzene     && bash run_sample.sh   # pentane / benzene mixture (SMILES)
cd sample/amorphous/ketoprofen          && bash run_sample.sh   # ketoprofen (SMILES)
cd sample/amorphous/ketoprofen_pubchem  && bash run_sample.sh   # ketoprofen via PubChem 3D SDF (CID 3825)
cd sample/amorphous/mixture_json        && bash run_sample.sh   # multi-component via JSON config
```

See [`docs/amorphous_tutorial.md`](docs/amorphous_tutorial.md) for the hands-on walk-through and [`sample/amorphous/ketoprofen/README.md`](sample/amorphous/ketoprofen/README.md) for an annotated run log of the ketoprofen build.

## Beyond this package

`abmptools` covers the generic layer: file formats, established-method
utilities, and the analysis and conversion tools around ABINIT-MP. Workflows
for methods still under active development are maintained separately, in a
private package (`moldeck`) that builds on this one. It currently covers:

| Area | What it automates |
|---|---|
| DPD / coarse-grained input | Building Cognac UDF and OCTA viewer inputs from a segmented structure, assigning interaction parameters, composing systems |
| FMO fragmentation | Proposing fragment splits for arbitrary molecules and exporting `segment_data` |
| Martini 3 builders | Peptide systems, and peptide-membrane PMF via umbrella sampling |
| Organic-crystal FMO | CIF to FMO inputs, lattice energy and deformation |
| Enhanced sampling / binding free energy | GENESIS gREST_SSCR and MM/GBSA end to end |
| Hydrogen-bond analysis | Detection, per-functional-group roles, lifetimes, and trajectory colouring |
| Peptide formulation | Solvated peptide + excipient systems end to end: build, equilibrate, then aggregation, contacts, secondary structure, SASA and release PMF |
| FMO interaction maps | Ligand-pocket contribution maps from IFIE/PIEDA |

The dependency runs one way — those workflows import `abmptools`, never the
reverse — so nothing here depends on having them.

It is **not publicly distributed**: it is proprietary and handed out
individually, on request and subject to approval. Contact the author below if
you have a use for it.

## How to cite

If you use ABMPTools in academic or scientific work, please cite the
project. GitHub's "Cite this repository" button on the repo home page
generates BibTeX / APA / etc. from [`CITATION.cff`](CITATION.cff).

A peer-reviewed publication and Zenodo DOI will be added on the first
release tag; until then, use:

> Okuwaki, K. (2026). *ABMPTools: a Python toolkit for ABINIT-MP
> Fragment Molecular Orbital pre/post-processing.*
> https://github.com/kojioku/abmptools

## Author

[Koji Okuwaki](mailto:koujioku81@gmail.com)
