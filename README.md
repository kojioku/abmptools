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

### Martini 3 Peptide CG System (`cg.peptide`)

- End-to-end Martini 3 peptide CG system builder: residue sequence (1-letter) + counts + box → `martinize2 -ff martini3001` で CG mapping → `gmx insert-molecules` で peptide 配置 → `gmx solvate` で Martini W → `gmx genion` で NaCl 中和 + 0.15 M → em/nvt/npt/md.mdp + `run.sh` を生成
- atomistic PDB は **tleap** (推奨; full sidechain) または extended-backbone fallback (tleap 不在時; 芳香族残基では NaN bead artifact の可能性あり、警告ログ)
- データクラス: `PeptideSpec` / `PeptideBuildConfig` (`@dataclass`、JSON 往復、YAML optional)
- CLI: `python -m abmptools.cg.peptide {build,validate,example}` (argparse)
- **Apache-2.0 互換のみ**: `vermouth-martinize` を `subprocess` で呼ぶのみ、改変・同梱なし。Martini 3 force field `.itp` は本パッケージ未同梱で `validate` サブコマンドが取得手順を表示
- `abmptools/cg/` namespace は MO-AAMD-CGMD マルチスケール基盤の CG 系統として新設。後続で `cg/polymer/` (polyply 経由) や `cg/smallmol/` (Auto-Martini 経由) を計画

### Peptide-Bilayer Umbrella Sampling (`membrane`)

- End-to-end PMF builder for peptide membrane permeation: bilayer + peptide + water + ions → AMBER (`ff19SB` + `Lipid21` + TIP3P / Joung-Cheatham) or CHARMM36 backend → semiisotropic NPT equilibration → z-pulling → per-window umbrella MDPs → `gmx wham` PMF
- packmol-memgen lipid placement (no CHARMM-GUI dependency); peptide built from one-letter sequence via `tleap`, capped with ACE/NME by default
- Two parameterisation routes:
  - `backend="amber"` — fully commercial-OK (`tleap` + `parmed` → GROMACS top/gro)
  - `backend="charmm36"` — MacKerell-free CHARMM36 parameter values via Klauda lab GROMACS port (`pdb2gmx`); CGenFF / CHARMM-GUI **forbidden by design** to keep the route commercial-clean
- GPU acceleration hook in the generated `run.sh` (`MDRUN_OPTS` env var)
- Bundled tutorial walks through poly-Ala 5-mer + POPC bilayer end-to-end

## Supported ABINIT-MP Versions

- ABINIT-MP v1: Rev.10–23
- ABINIT-MP v2: Rev.4–8

## Installation

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

## Quick Start

```bash
# Extract IFIE for fragment 10, within 8 Å
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i calculation.log

# Generate AJF input from PDB
python -m abmptools.generateajf -i protein.pdb -basis 6-31G* --method MP2

# Convert log to CPF
python -m abmptools.log2cpf -i calculation.log -o output.cpf

# Create DIFIE-averaged CPF from trajectory
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 10 1 -f 1-100 -np 4

# Convert UDF to GROMACS
python -m abmptools.udf2gro.cli -i system.udf -o output

# Convert GROMACS to UDF
python -m abmptools.gro2udf.cli -i system.gro -t system.top -o output.udf

# Build an amorphous mixture from SMILES (50 ketoprofen molecules, density 0.8 g/cm^3)
python -m abmptools.amorphous --smiles "OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1" \
    --name ketoprofen --n_mol 50 --density 0.8 --output_dir ./ketoprofen

# Or use an external 3D SDF (e.g. from PubChem) as the initial conformer
python -m abmptools.amorphous --mol ketoprofen_pubchem_cid3825.sdf \
    --name ketoprofen --n_mol 50 --density 0.8 --output_dir ./ketoprofen_pubchem

# Or let abmptools fetch the 3D SDF straight from PubChem (1.15.3+)
python -m abmptools.amorphous --pubchem_cid 3825 \
    --name ketoprofen --n_mol 50 --density 0.8 --output_dir ./ketoprofen_pubchem

# Build a Martini 3 peptide CG system (KGG x5 + RGG x5 in 10 nm box)
python -m abmptools.cg.peptide example > kgg.json   # 最小 example
python -m abmptools.cg.peptide validate --config kgg.json --ff-dir ./ff
python -m abmptools.cg.peptide build    --config kgg.json --ff-dir ./ff -o ./out
# 上記で out/run/run.sh が生成される。Martini 3 .itp は cgmartini.nl から
# 別途取得が必要 (本パッケージ未同梱)。詳細は abmptools/cg/peptide/README.md。

# Build a peptide-bilayer Umbrella Sampling system (AMBER backend)
python - <<'PY'
from abmptools.membrane import (MembraneConfig, MembraneUSBuilder,
    LipidSpec, PeptideSpec, USProtocol)
cfg = MembraneConfig(
    backend="amber",
    lipids=[LipidSpec(resname="POPC", n_per_leaflet=32)],
    peptide=PeptideSpec(name="aa5", sequence="AAAAA"),
    output_dir="./membrane_run", seed=42,
    umbrella=USProtocol(z_min_nm=-1.5, z_max_nm=+1.5, window_spacing_nm=0.25),
)
print(MembraneUSBuilder(cfg).build()["run_script"])
PY
```

Use `-h` with any module for full option details.

## Documentation

- **[User Manual](docs/ABMPTools-user-manual.md)** — CLI options, output formats, and workflow examples
- **[Architecture](docs/architecture.md)** — Class hierarchy and design overview
- **[Developer Quickstart](docs/dev_quickstart.md)** — Setup and code conventions
- **[I/O Spec](docs/io_spec.md)** — File format specifications
- **[gro2udf](docs/gro2udf.md)** / **[udf2gro](docs/udf2gro.md)** — GROMACS ↔ OCTA conversion
- **[geomopt](docs/geomopt.md)** / **[amorphous](docs/amorphous.md)** — Optimization and structure building
- **[membrane](docs/membrane.md)** / **[tutorial_membrane_us](docs/tutorial_membrane_us.md)** — Peptide-bilayer umbrella-sampling PMF (reference + step-by-step ops tutorial)

## Testing

```bash
pytest tests/ -v                     # 671 tests across 30 files
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

# Amorphous structure builder samples
cd sample/amorphous                    && bash run_sample.sh   # pentane / benzene mixture (SMILES)
cd sample/amorphous/ketoprofen_pubchem && bash run_sample.sh   # ketoprofen via PubChem 3D SDF (CID 3825)
```

See [`sample/amorphous/ketoprofen/README.md`](sample/amorphous/ketoprofen/README.md) for a step-by-step walk-through of the ketoprofen amorphous workflow (SMILES input + 5-stage MD + VMD post-processing).

## Author

[Koji Okuwaki](mailto:koujioku81@gmail.com)
