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
│   ├── amorphous/          ← アモルファス系構築 (packmol, OpenMM)
│   ├── membrane/           ← 膜透過 PMF 用 US ビルダー (packmol-memgen, AMBER/CHARMM36)
│   ├── geomopt/            ← 構造最適化 (PySCF, MACE)
│   └── f90/                ← Fortran extension
│       ├── src/
│       │   └── readifiepiedalib.f90
│       └── bin/
│           └── readifiepiedalib.so
├── tests/                  ← pytest テストスイート (658 tests, 28 files)
│   ├── conftest.py
│   ├── test_cli_scripts.py ← 全14 CLIスクリプトのargparseテスト
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
│   ├── amorphous.md
│   ├── membrane.md
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

pytest ベースのテストスイート。658テスト、28ファイル。

| カテゴリ | ファイル数 | テスト数 | 概要 |
|---------|-----------|---------|------|
| コアクラス | 7 | 300 | molcalc, mol_io, abinit_io, pdb_io, udf_io, udfrm_io, udfcreate |
| 複合クラス | 2 | 72 | setfmo, anlfmo |
| マネージャ | 2 | 64 | cpfmanager, logmanager |
| スタンドアロン | 2 | 61 | readcif, getifiepieda |
| サブパッケージ | 14 | 127 | core, gro2udf, udf2gro, amorphous, geomopt |
| CLIスクリプト | 1 | 50 | 全14 CLIの argparse テスト |

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
| `amorphous.md` | アモルファス系構築 |
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
