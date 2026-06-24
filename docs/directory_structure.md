# Directory Structure

## Top-Level Layout

```
abmptools/                  ← Repository root
├── abmptools/              ← Python package (27+ modules)
│   ├── __init__.py
│   ├── *.py
│   ├── core/               ← 共通データモデル (system_model)
│   ├── gro2udf/            ← GROMACS → OCTA UDF 変換
│   ├── udf2gro/            ← OCTA UDF → GROMACS 変換
│   ├── udfcharge/          ← 単分子 UDF → バルク UDF 電荷転写
│   ├── amorphous/          ← アモルファス系構築 (packmol, OpenMM)
│   ├── membrane/           ← AA 膜透過 PMF 用 US ビルダー (packmol-memgen, AMBER/CHARMM36)
│   ├── cg/                 ← Martini 3 CG 系統 (1.18.0+)
│   │   ├── peptide/        ← Martini 3 ペプチド CG ビルダー (vermouth-martinize)
│   │   └── membrane/       ← Martini 3 ペプチド-膜 PMF (insane + umbrella)
│   ├── genesis/            ← GENESIS 系統 (1.20.0+)
│   │   ├── grest/          ← gREST_SSCR replica exchange (atdyn/spdyn + remd_convert)
│   │   └── mmgbsa/         ← GENESIS MM/GBSA single-point ΔG_bind (atdyn + acpype + tleap)
│   ├── fragmenter/         ← FMO 自動フラグメント分割 (rdkit + ipywidgets) (1.21.0+)
│   ├── geomopt/            ← 構造最適化 (PySCF, MACE)
│   └── f90/                ← Fortran extension
│       ├── src/
│       │   └── readifiepiedalib.f90
│       └── bin/
│           └── readifiepiedalib.so
├── tests/                  ← pytest テストスイート (1530 tests collected, 80+ files)
│   ├── conftest.py
│   ├── test_cli_scripts.py ← 全14 CLIスクリプトのargparseテスト
│   ├── test_grest_*.py     ← genesis.grest (1.20.0、159 tests + 1 slow)
│   ├── test_mmgbsa_*.py    ← genesis.mmgbsa (1.22.0、117 tests + 1 slow)
│   ├── fragmenter/         ← fragmenter (1.21.0、14 tests)
│   ├── fixtures/           ← 共通テスト fixture (gbsa_logs/ 等)
│   └── test_*.py           ← 各モジュールの単体テスト
├── docs/                   ← Documentation
│   ├── ABMPTools-user-manual.md
│   ├── overview.md
│   ├── architecture.md
│   ├── directory_structure.md
│   ├── dev_quickstart.md
│   ├── io_spec.md
│   ├── dataflow.md
│   ├── dependencies.md
│   ├── parallelization.md
│   ├── faq.md
│   ├── gro2udf.md
│   ├── udf2gro.md
│   ├── udfcharge.md
│   ├── tutorial_udfcharge.md
│   ├── amorphous.md
│   ├── membrane.md
│   ├── tutorial_membrane_us.md
│   ├── cg_membrane.md
│   ├── tutorial_cg_membrane_us.md
│   ├── geomopt.md
│   ├── qmopt.md
│   └── licenses_third_party.md
├── sample/                 ← Sample data & run scripts
│   ├── convertcpf/
│   ├── generate_difie/
│   ├── gbsa/
│   ├── generateajf/
│   ├── log2config/
│   ├── log2cpf/
│   └── rmsd/
├── tips/                   ← Auxiliary MD post-processing scripts
├── build/                  ← Build artifacts
├── setup.py                ← Package setup (setuptools)
├── Makefile                ← Fortran compilation
├── CHANGELOG.md            ← Version history
├── LICENSE                 ← License
├── README.md               ← Project overview (Japanese)
├── input_param             ← Parameter template for pdb2fmo
└── segment_data.dat        ← Fragment definition format example
```

## `abmptools/` — Python Package

### Core I/O Modules (Base Classes)

| File | Description |
|------|-------------|
| `__init__.py` | Package initialization; exports `CPFManager`, `LOGManager` and other public classes |
| `mol_io.py` | Base molecular I/O class — XYZ reading, coordinate conversions |
| `molcalc.py` | Coordinate math — distances, rotations, Euler angles, PBC handling, LAMMPS parsing |
| `abinit_io.py` | ABINIT-MP I/O base — fragment definitions, geometry management (inherits `mol_io`) |
| `pdb_io.py` | PDB file parser — residue/atom parsing, multiple read modes (inherits `abinit_io`) |

### Data Managers

| File | Description |
|------|-------------|
| `cpfmanager.py` | **CPFManager** class — parse/write CPF files (v4.201/7.0/10/23), pandas DataFrame storage |
| `logmanager.py` | **LOGManager** class — parse ABINIT-MP log files, extract conditions/fragments/energies |

### Analysis Modules

| File | Description |
|------|-------------|
| `getifiepieda.py` | Primary IFIE/PIEDA CLI — single, matrix, time-series, multi-sample modes |
| `anlfmo.py` | Advanced FMO analysis — log parsing, fragment interaction, multi-sample (inherits `pdb_io`) |
| `generate_difie.py` | Dynamic IFIE generation — average CPFs across trajectory snapshots |
| `cpf2ifielist.py` | CPF to formatted IFIE list — filter by fragment range, CSV output |
| `getcharge.py` | Charge extraction from logs — NBO and other charge models |

### FMO Setup Modules

| File | Description |
|------|-------------|
| `generateajf.py` | AJF template generation CLI — method, basis set, solvation, PIEDA options |
| `setfmo.py` | FMO setup orchestrator — cut modes, solute/solvent, fragment assignment |
| `pdb2fmo.py` | PDB to FMO converter — config-based fragment assignment |
| `addsolvfrag.py` | Add solvation molecules to AJF — template-based solvent addition |

### File Conversion Modules

| File | Description |
|------|-------------|
| `log2cpf.py` | Log to CPF conversion CLI |
| `log2config.py` | Log to fragment config dictionary |
| `convertcpf.py` | CPF version converter with fragment filtering |
| `readcif.py` | CIF (crystallographic) to Cartesian coordinate conversion |
| `ajf2config.py` | AJF to `segment_data.dat` config converter |
| `ajfserial.py` | Generate numbered AJF files from a template |
| `pdbmodify.py` | PDB editing — move, rename, chain assignment, residue refresh |

### UDF/MD File Handling

| File | Description |
|------|-------------|
| `udf_io.py` | UDF (OCTA COGNAC) file reader — molecular coordinates from MD trajectories |
| `udfcreate.py` | UDF parameter setup — MD simulation configuration |
| `udfrm_io.py` | UDF to PDB/XYZ conversion with PBC handling (inherits `udf_io`) |
| `udf2fmo.py` | UDF to FMO converter CLI |

## `abmptools/f90/` — Fortran Extension

- **`src/readifiepiedalib.f90`** (219 lines) — High-performance IFIE/PIEDA reader for log files.
- **`bin/readifiepiedalib.so`** (21 KB) — Pre-compiled shared library loaded at runtime by `getifiepieda.py`.
- Compiled via `Makefile` using `gfortran`. Falls back to pure Python if unavailable.

## `tests/` — テストスイート

pytest ベースのテストスイート。**1530 tests collected** (1.22.0 時点、`pytest --collect-only -q` 結果)。

| カテゴリ | ファイル数 | テスト数 | 概要 |
|---------|-----------|---------|------|
| コアクラス | 7 | 300 | molcalc, mol_io, abinit_io, pdb_io, udf_io, udfrm_io, udfcreate |
| 複合クラス | 2 | 72 | setfmo, anlfmo |
| マネージャ | 2 | 64 | cpfmanager, logmanager |
| スタンドアロン | 2 | 61 | readcif, getifiepieda |
| サブパッケージ (AA) | 14+ | 127+ | core, gro2udf, udf2gro, **udfcharge** (10), amorphous, geomopt, **membrane** |
| **サブパッケージ (CG, 1.18.0+)** | 23 | **287** | `cg.peptide` 117 + `cg.membrane` 170 (+ 4 slow integration) |
| **サブパッケージ (GENESIS, 1.20.0+)** | 17 | **276** | `genesis.grest` 159 + `genesis.mmgbsa` 117 (+ 2 slow integration、tools 不在時 auto-skip) |
| **`fragmenter` (1.21.0+)** | 2 | **14** | basic + polymer の 2 ファイル、RDKit 必須 |
| CLIスクリプト | 1 | 50 | 全 CLI の argparse テスト |

詳細は `tests/TEST_COVERAGE.md` を参照。

## `docs/` — Documentation

| ファイル | 概要 |
|---------|------|
| `ABMPTools-user-manual.md` | ユーザーマニュアル (CLI オプション、出力例、ワークフロー) |
| `overview.md` | プロジェクト概要 |
| `architecture.md` | アーキテクチャ設計 (継承チェーン、データフロー) |
| `directory_structure.md` | ディレクトリ構成 (本ファイル) |
| `dev_quickstart.md` | 開発者向けクイックスタート |
| `io_spec.md` | I/O仕様 (ファイルフォーマット) |
| `dataflow.md` | データフロー図 |
| `dependencies.md` | 依存関係 |
| `parallelization.md` | 並列処理 |
| `faq.md` | FAQ |
| `gro2udf.md` / `udf2gro.md` | GROMACS ↔ OCTA 変換ドキュメント |
| `udfcharge.md` / `tutorial_udfcharge.md` | 単分子 UDF → バルク UDF 電荷転写 |
| `amorphous.md` | アモルファス系構築 |
| `membrane.md` | 膜透過 PMF 用 US ビルダー (AA reference; AMBER/CHARMM36) |
| `tutorial_membrane_us.md` | 同上の step-by-step 操作チュートリアル |
| `cg_membrane.md` | Martini 3 ペプチド-膜 PMF ビルダー (CG reference) |
| `tutorial_cg_membrane_us.md` | 同上の step-by-step (smoke 5-6 分 + production 45 分) |
| `cg_peptide.md` | Martini 3 ペプチド CG ビルダー (CG reference、1.18.0+) |
| `peptide_builders.md` | 3 種ペプチドビルダー横断比較 (AA membrane / CG peptide / CG membrane) |
| `grest.md` / `tutorial_grest.md` | GENESIS gREST_SSCR ビルダー (1.20.0+、subsystem reference + step-by-step) |
| `mmgbsa.md` / `tutorial_mmgbsa.md` | GENESIS MM/GBSA single-point ΔG_bind (1.22.0+、subsystem reference + step-by-step) |
| `fragmenter.md` | FMO automatic fragment splitter (1.21.0+、graph-diameter MW walk + canonical SMILES + γ polymer) |
| `geomopt.md` / `qmopt.md` | 構造最適化 |
| `licenses_third_party.md` | サードパーティライセンス |

## `sample/` — Sample Data & Workflows

Each subdirectory contains input data and a `run.sh` script:

| Directory | Workflow |
|-----------|----------|
| `convertcpf/` | CPF version conversion and fragment filtering |
| `generate_difie/` | DIFIE averaged CPF from trajectory (TrpCage, CS4 examples) |
| `gbsa/` | GBSA solvation setup |
| `generateajf/` | AJF template generation from PDB |
| `log2config/` | Log to fragment configuration extraction |
| `log2cpf/` | Log to CPF conversion |
| `rmsd/` | RMSD analysis post-processing |

## Build & Configuration Files

| File | Purpose |
|------|---------|
| `setup.py` | setuptools configuration; triggers `make` for Fortran compilation |
| `Makefile` | Compiles `readifiepiedalib.f90` → `readifiepiedalib.so` |
| `input_param` | Parameter template for `pdb2fmo` (cut mode, solute, criteria, etc.) |
| `segment_data.dat` | Fragment definition format (Python dict / binary) |
