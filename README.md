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

### H-bond Analyzer for COGNAC Trajectories (`hbond`)

- **OCTA COGNAC `.udf` / `.bdf` 専用** の H-bond 解析サブパッケージ。非晶質 MD で
  カルボキシル基同士の dual H-bond (環状二量体) と COOH→アミド C=O の single H-bond
  を区別して数え、gourmet で 3 色可視化できる UDF を出力する
- 検出基準は **Luzar-Chandler** (`d(D-A) ≤ 3.5 Å`, `∠(D-H-A) ≥ 120°`) を default、
  **strict** (`d(H-A) ≤ 2.5 Å`, `∠ ≥ 150°`)、**custom** (任意閾値) も選択可能。
  直交 cubic box の minimum image PBC 対応
- 官能基自動検出: GAFF2 atomtype (`c`/`oh`/`ho`/`o`/`n`) + bond graph で carboxyl /
  amide / hydroxyl を機械的に同定 (SMARTS 不要)。Tertiary amide 判定付き
- 3 経路: CLI (`python -m abmptools.hbond <bdf> -o prefix`) / Python API
  (`Analyzer`, `AnalyzerConfig`) / Jupyter ipywidgets UI (`open_panel(bdf_path)`、
  RDKit 2D 構造図上で carboxyl/amide ハイライト + matplotlib count plot)
- 出力: per-record summary CSV (官能基単位の dual/single/free 数 + 比率) +
  per-functional-group classification CSV + H-bond pair CSV + colored BDF
  + Mol_Name 維持 plain BDF + count vs record PNG
- 色付け: `<prefix>_colored.bdf` は `Set_of_Molecules.molecule[i].Mol_Name` を
  3 グループ (`IMC_DUAL` / `IMC_SINGLE` / `IMC_FREE`) にリネームし、
  `Draw_Attributes.Molecule[]` に named color (Red/Blue/Gray) を書き込む
  (**GOURMET Draw_Attributes の color は select 型 9 色名のみ、RGBA tuple 不可**)。
  `<prefix>.bdf` は Mol_Name 維持コピーで **J-OCTA プリ描画でも分子が空表示にならない**
- バンドル sample (IMC): `sample/hbond/imc_amorphous/` (非晶質インドメタシン
  T=450 K、125 分子。**4-species per-COOH**: dual=10 / chain=41 / single=38 /
  free=36; per-amide: accept=49 / free=76。4-species の枠組みは Yuan et al. 2015
  Mol. Pharm. 12, 4518 の固体 NMR 帰属に倣う。掲載は MD 値のみ・論文の図表は非再現)
- バンドル sample (PVA, v1.28+): `sample/amorphous/pva_amorphous/` — PVA 10-mer
  × 30、OpenFF Sage + AM1-BCC + 5-stage MD + xtc→UDF + hbond generic mode の
  end-to-end 例 (平均 198.8 H-bonds/record、ratio_donor_busy=65.2%)
- 依存: `pip install abmptools[hbond]` (matplotlib for plot)、Jupyter UI を使うなら
  `[jupyter]` + `[fragmenter]` (rdkit) を併用。UDFManager は OCTA に同梱
- **v1.26.0+ 拡張**: FF 抽象化 (GAFF2/OPLS-AA/CHARMM36/OpenFF)、任意官能基対選択
  (donor: carboxyl/amide_donor/amine_donor/hydroxyl × acceptor: carboxyl_O/amide_O/
  hydroxyl_O/ether_O)、secondary amide N-H donor 対応、multi-record lifetime +
  Luzar-Chandler 自己相関 `C(t)` + τ_HB 算出
- **v1.27.0 候補**: per-functional-group classification (1 分子内に複数 COOH/amide
  がある場合に役割が混在するケースに正しく対応)、`<prefix>.bdf` 併出
  (J-OCTA プリ描画用)、`<prefix>_classification.csv` 新規追加、
  **4-species 分類** (dual/chain/single/free、枠組みは Yuan 2015 IMC NMR 帰属に
  倣う)。MD species-fraction plot script (`plot_nmr_comparison.py`、MD 値のみ・
  論文図表は非再現) 同梱
- **v1.28.0 候補**: **generic mode** (`--classify-mode generic`) で COOH を持た
  ない任意系 (PVA / peptide / アルコール / 混合系) の donor-type × acceptor-type
  pair 統計に対応 (新規 `<prefix>_pair_stats.csv` + atom-role `Donor/Acceptor/Both`
  色付け)。**element + bond-graph fallback** (default ON) で OpenFF SMIRNOFF UDF
  (per-atom unique `MOL0_X`) を **antechamber GAFF patch 不要**で直接解析可能。
  Jupyter UI に `Mode:` dropdown 追加。J-OCTA Viewer 用 plain `.py` script 併出
  (autorun crash 回避用)、`<prefix>.bdf` Attributes に `hbond=Dual/Chain/Single/
  Free/Accept` (imc mode) または `Donor/Acceptor/Both` (generic mode) を append
  (J-OCTA Attribute フィルタでカテゴリ可視化)

### Peptide-Bilayer Umbrella Sampling (`membrane`)

- End-to-end PMF builder for peptide membrane permeation: bilayer + peptide + water + ions → AMBER (`ff19SB` + `Lipid21` + TIP3P / Joung-Cheatham) or CHARMM36 backend → semiisotropic NPT equilibration → z-pulling → per-window umbrella MDPs → `gmx wham` PMF
- packmol-memgen lipid placement (no CHARMM-GUI dependency); peptide built from one-letter sequence via `tleap`, capped with ACE/NME by default
- Two parameterisation routes:
  - `backend="amber"` — fully commercial-OK (`tleap` + `parmed` → GROMACS top/gro)
  - `backend="charmm36"` — MacKerell-free CHARMM36 parameter values via Klauda lab GROMACS port (`pdb2gmx`); CGenFF / CHARMM-GUI **forbidden by design** to keep the route commercial-clean
- GPU acceleration hook in the generated `run.sh` (`MDRUN_OPTS` env var)
- Bundled tutorial walks through poly-Ala 5-mer + POPC bilayer end-to-end
- Sample driver/config (Phase D = L9 verification, both backends in parallel):
  `sample/membrane/amber_phaseD/` (AMBER ff19SB + Lipid21 + TIP3P, PMF +86.7 kJ/mol) and `sample/membrane/charmm_phaseD/` (CHARMM36 Klauda port, PMF +97.9 kJ/mol, Δ-11.3 kJ/mol vs AMBER — typical FF gap)

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

## Testing

```bash
pytest tests/ -v                     # 1300 tests collected (2.9.0 時点)
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
