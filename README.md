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

### Martini 3 Peptide-Membrane PMF (`cg.membrane`)

- End-to-end Martini 3 peptide-bilayer **umbrella sampling** builder: cg.peptide で M3 CG ペプチド生成 → `insane` (GPL-2.0、subprocess only) で POPC bilayer に埋め込み → topology composer で 4 ITP includes / `Protein → molecule_0` / `NA+/CL- → NA/CL` 正規化 → `index.ndx` (Bilayer / Peptide / W / NA / CL / Non_Bilayer) → em / nvt / npt-semiisotropic + pull-direction-periodic + 13 window static umbrella + `run.sh` を生成
- AA 系の `membrane` (CHARMM36 / Lipid21) と並走する CG 版。`abmptools.membrane.{pulling,pmf,mdp_us_protocol}` の generic helpers を **import 経由で再利用** (コード重複ゼロ)
- データクラス: `MembraneCGBuildConfig` / `LipidMix` / `PeptideMembraneSpec` / `EquilibrationCGProtocol` / `PullingCGProtocol` / `UmbrellaCGProtocol` (5 段 nested JSON 往復)
- CLI: `python -m abmptools.cg.membrane {build,validate,example,make-windows,wham}` (argparse)
- Default umbrella: 13 windows (z = -1.5 to +1.5 nm), k = 1000 kJ/mol/nm², 1 ns/window (50,000 steps × dt=20 fs); pulling 5 ns × 1 nm/ns
- **GPL-2.0 / Apache-2.0 / MIT 互換のみ**: `insane` (GPL-2.0) と `vermouth-martinize` (Apache-2.0) は subprocess only -- abmptools 自体は MIT のまま (mere aggregation)

### GENESIS gREST_SSCR (`genesis.grest`)

- End-to-end **GENESIS gREST_SSCR** builder + analysis: protein PDB → `tleap` で AMBER ff19SB + TIP3P → REST 残基確定 (explicit / around モード両対応) → 温度ラダー生成 (auto geometric / manual 両対応) → 4 つの GENESIS `.inp` (`step1_minimize` / `step2_equilibrate` / `step3_grest` / `step5_remd_convert`) + `mpirun -np N spdyn` 用 `run.sh` を生成
- 解析: `analyze` サブコマンドで replica transition plot / acceptance ratio plot / `remd_convert` で param sort / 1D 距離 PMF (`-kT log P(r)`) を実装
- データクラス: `GrestBuildConfig` / `RESTSelectionSpec` / `ReplicaTemperatureSpec` / `MinimizationStage` / `EquilibrationStage` / `GrestStage` (5 段 nested JSON 往復)
- CLI: `python -m abmptools.genesis.grest {build,validate,example,analyze}` (argparse)
- **LGPL-3.0-or-later 互換**: GENESIS (`spdyn` / `atdyn` / `remd_convert`) は subprocess only -- abmptools 自体は MIT (1.20.0) → 将来 Apache-2.0 化後も互換 (mere aggregation per LGPL §5/§6)
- `abmptools/genesis/` は GENESIS 系統 namespace の最初の occupant。後続で `genesis/reus/` / `genesis/fep/` を計画

### GENESIS MM/GBSA (`genesis.mmgbsa`)

- End-to-end **GENESIS MM/GBSA** ΔG_bind builder + analysis: protein-ligand complex PDB + ligand 残基番号 → Biopython で receptor/ligand 分割 → acpype (GAFF/GAFF2 + AM1-BCC) + tleap で 3 系 (complex/ligand/receptor) の AMBER prmtop+inpcrd → `mpirun -np 1 atdyn` で `[ENERGY] implicit_solvent=GBSA` 単フレーム評価 → `[STEP4]` log パース → `ΔG_bind = (E+S)_complex - (E+S)_ligand - (E+S)_receptor` を CSV + 棒グラフ出力
- データクラス: `MMGBSABuildConfig` / `TargetSpec` / `ForceFieldSet` / `LigandParameterization` / `EnergyProtocol` / `MinimizationProtocol` (5 段 nested JSON 往復)
- CLI: `python -m abmptools.genesis.mmgbsa {build,validate,example,divide,parameterize,run,analyze,pipeline}` (7 sub-command + folder-mode shortcut `-i / -r / -c` で POC 互換)
- 力場 default: AMBER **ff14SB** + DNA.OL15 + RNA.OL3 + TIP3P + GAFF/GAFF2 (POC 通り、grest の ff19SB とは意図的に別)
- **ΔG_bind 計算**: GENESIS `ENERGY` 列は `U = U_FF + ΔG_solv` の合計 (doc 05_Energy.rst:564) なので `ΔG_bind = E_c - E_l - E_r` で全 MM/GBSA 寄与込み。POC `4_analyse.py` は等価な分解形 `(egas + S)_c - (egas + S)_l - (egas + S)_r` を使用 (`egas = E - S`)、両者は代数的に同一。`compute_dg_bind` (合計形) + `compute_dg_components` (分解報告 `{dg_mm, dg_solv, dg_bind}`) を提供
- **LGPL-3.0+ / GPL-3.0 互換**: GENESIS (LGPL-3.0+) と acpype (GPL-3.0) は subprocess only -- abmptools 自体は MIT (1.22.0) → 将来 Apache-2.0 化後も互換 (mere aggregation)

### FMO Fragment Auto-splitter (`fragmenter`)

- End-to-end **FMO 自動フラグメント分割** ツール: PDB → RDKit で 4-strategy bond perception (proximity → CONECT-only → `sanitize=False` → `obabel`) → heavy-atom-only canonical SMILES でグループ化 → graph diameter (heavy-atom-only 2-pass BFS) で主鎖検出 → 累積 MW + 側鎖 MW を加算しながら `target_mw` (default 200 g/mol) ごとに C-C 切断候補を提案 → `RWMol` で破壊的に bond 削除 → `log2config` 互換 `segment_data.dat` 出力 (pdb2fmo がそのまま読む)
- 対象は **小分子 / 脂質 / ポリマー** (タンパク質・DNA は対象外、既存 `log2config` 経路へ流す方針)
- フィルタ: 環内 / 多重結合 / ヘテロ隣接 (N/O/S/P/F/Cl/Br/I) を除外、いずれも config で off 可能
- ポリマー γ 経路 (`declare_same_pattern`): 異なる SMILES (PE N=10 / N=11 等) を明示的に同一視、master (最も cut の多い group) のパターンを atom-path-index 対応で短い chain にも転送
- 3 つの UI 経路:
  - **A: Jupyter UI** (`AutoFragmenter.from_pdb` + `open_panel`、ipywidgets dropdown / SVG / checkbox)
  - **C: ヘッドレス CLI** (`python -m abmptools.fragmenter {suggest,apply,example}` + SVG+JSON review bundle)
  - **API**: 関数 chain (`load_pdb_molecules` → `group_by_smiles` → `suggest_cuts_for_groups` → `export_to_system`)
- データクラス: `FragmenterConfig` / `CutSite` / `MoleculeGroup` / `FragmentResult` (JSON 往復可能)
- 14 unit tests (10 basic + 4 polymer)、5 実 PDB 検証ケース (ketoprofen / PE N=20 / PP N=10 / propane×5+acetone×3 / antibody+ligand)
- 依存: `pip install abmptools[fragmenter]` (rdkit-pypi >= 2022.09、BSD-3-Clause)、Jupyter UI を使うなら `[jupyter]` extras も

### Crystal-FMO Pipeline (`crystal`)

- End-to-end **organic-crystal FMO** workflow: CIF → supercell PDB → fragment cut around a target solute → ABINIT-MP AJF (with full-precision `&XYZ` block) → optional `--run-local` invocation → `getifiepieda` postprocessing (IFIE/PIEDA + MonomerEnergy)
- Two CIF backends: `engine='legacy'` (the historical hand-rolled parser, byte-equivalent with v1.22.0 csp7 outputs) and `engine='ase'` (ASE space-group expansion, arbitrary `layer`)
- 8-subcommand CLI: `abmp-crystal {expand,fragment,jobs,pipeline,postproc,nearest,validate,example}` driven by a single YAML/JSON config (`CrystalBuildConfig` + 7 leaf dataclasses)
- HPC scheduler templates: PJM / SLURM / PBS / `local`. The `local` scheduler combined with `--run-local` invokes `abinitmp` directly per AJF (smoke / single-shot reference runs)
- Numeric reference frozen for csp7 R00001 layer3 HF/6-31G (9h on 1 core, abinitmp v2r8) as `frag1-dimer-es-false-{ifiesum,ifiedt}.csv` from the in-tree `getifiepieda` post-processor — no parallel parser duplication
- `abmptools.anlfmo` gained HF-log support along the way (5 defensive edits, MP2 production path unchanged)
- Bundled tutorial: `docs/tutorial_crystal_fmo.md` (9 sections, including reference-establishment recipe); design notes: `docs/crystal.md`; verification matrix: `docs/crystal_verification.md`
- Sample driver/config: `sample/crystal/csp7_smoke/` (cif and `UNK.ajf` template are private and live in `abmptools-sample` — staged automatically when `ABMPTOOLS_SAMPLE_DIR` is set)
- Dependencies: `pip install abmptools[crystal]` (`ase >= 3.22` / `pyyaml >= 6.0`); ABINIT-MP for `--run-local` only

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

# Martini 3 peptide-membrane PMF (umbrella sampling) — abmptools.cg.membrane
python -m abmptools.cg.membrane example > kgg_popc.json
python -m abmptools.cg.membrane validate --config kgg_popc.json --ff-dir ./ff
python -m abmptools.cg.membrane build    --config kgg_popc.json --ff-dir ./ff -o ./out
bash ./out/run.sh                                    # em → nvt → npt → pull → 13 windows → wham
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
- **[membrane](docs/membrane.md)** / **[tutorial_membrane_us](docs/tutorial_membrane_us.md)** — Peptide-bilayer umbrella-sampling PMF (AA, CHARMM36 / Lipid21)
- **[cg_membrane](docs/cg_membrane.md)** / **[tutorial_cg_membrane_us](docs/tutorial_cg_membrane_us.md)** — Martini 3 peptide-bilayer PMF (CG, 30-100× faster than AA, KGG-POPC smoke 5 min / production 45 min)
- **[cg_peptide](docs/cg_peptide.md)** — Martini 3 peptide CG builder (peptide-only in water box, sub-called by `cg.membrane` or standalone)
- **[peptide_builders](docs/peptide_builders.md)** — Selection guide across 3 peptide-from-sequence builders (AA membrane / CG peptide / CG membrane)
- **[fragmenter](docs/fragmenter.md)** — FMO automatic fragment splitter for small molecules / lipids / polymers (canonical SMILES grouping + C-C MW walk + Jupyter UI / headless CLI; v1.21.0+)

## Testing

```bash
pytest tests/ -v                     # 1530 tests collected (1.22.0 時点)
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
