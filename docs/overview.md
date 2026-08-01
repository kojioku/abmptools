# ABMPTools Overview

## What is ABMPTools?

ABMPTools is a Python toolkit for pre-processing, post-processing, and analysis of Fragment Molecular Orbital (FMO) calculations performed with [ABINIT-MP](https://fmodd.jp/member_contents/manual_ABINIT-MP/). It provides:

- **IFIE/PIEDA analysis** — Extract and visualize inter-fragment interaction energies from calculation logs and CPF files.
- **CPF parsing & conversion** — Read, write, filter, and version-convert CPF (Coordinate Property File) binary/text files.
- **FMO input generation** — Automatically generate AJF input files from PDB structures with fragment assignment.
- **File format conversion** — Convert between PDB, LOG, CPF, UDF (OCTA COGNAC), CIF, and XYZ formats.
- **Dynamic IFIE (DIFIE)** — Average IFIE data across MD trajectory snapshots into a single CPF with mean/standard-deviation statistics.

## Supported ABINIT-MP Versions

- ABINIT-MP v1: Rev.10–23
- ABINIT-MP v2: Rev.4–8
- CPF format versions: 4.201, 7.0 (MIZUHO), 10, 23

## Installation

Editable install is recommended for day-to-day use and development:

```bash
pip install -e .
```

Non-editable install (e.g. for production deployment):

```bash
pip install .
```

`--user` is usually unnecessary; pip handles both virtual environments and system Python appropriately. Installation runs `make` to compile the optional Fortran shared library (`readifiepiedalib.so`). If `gfortran` is not available, the install still succeeds — the Fortran acceleration is optional.

### Requirements

- **Required**: numpy, pandas
- **Optional**: UDFManager (OCTA COGNAC), gfortran (Fortran acceleration), OpenBabel (`obabel` CLI)

## Quick Example

```bash
# Generate an AJF input from a PDB file
python -m abmptools.generateajf -i protein.pdb -basis 6-31G* --method MP2

# Extract IFIE for fragment 10, within 8 Å distance
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i calculation.log

# Convert a log file to CPF format
python -m abmptools.log2cpf -i calculation.log -o output.cpf

# Create a DIFIE-averaged CPF from trajectory snapshots
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 10 1 -f 1-100 -np 4
```

## Key Capabilities

| 分類 | モジュール | 概要 |
|---|---|---|
| FMO 入力生成 | `generateajf` / `pdb2fmo` / `addsolvfrag` | PDB から ajf を作る。溶媒フラグメント付加も |
| FMO 結果解析 | `anlfmo` / `getifiepieda` / `generate_difie` | IFIE / PIEDA の抽出・集計、trajectory 平均 |
| ファイル変換 | `cpfmanager` / `logmanager` / `abinit_io` | CPF / log / ajf の読み書き |
| 構造最適化 | `geomopt` / `qmopt` | 座標最適化のラッパ |
| アモルファス構築 | `amorphous` | packmol + OpenFF による非晶質系の構築 |
| 水素結合解析 | `hbond` | 4 力場対応、寿命・4-species 分類・2D diagram |
| 製剤系 MD | `formulation` | ペプチド製剤の構築・解析 |
| 膜透過 PMF | `membrane` | umbrella sampling による PMF (AMBER / CHARMM36) |
| OCTA 連携 | `gro2udf` / `udf2gro` / `udfcharge` | GROMACS ⇄ COGNAC UDF の変換、電荷転写 |
| trajectory 操作 | `trajectory` | gmx trjconv / energy のラッパ |

## Where to Start Reading

1. **`docs/overview.md`** (本ページ) — 全体像とクイックスタート
2. **`docs/architecture.md`** — モジュール構成と設計方針
3. **`docs/directory_structure.md`** — ディレクトリの地図
4. **`docs/dataflow.md`** — 入出力の流れ
5. **`docs/io_spec.md`** — ファイル形式の仕様
6. サブシステム別リファレンス — `amorphous.md` / `hbond.md` / `formulation.md` /
   `membrane.md` / `gro2udf.md` / `udf2gro.md` / `udfcharge.md` / `trajectory.md` /
   `geomopt.md` / `qmopt.md`
7. チュートリアル — `amorphous_tutorial.md` / `tutorial_hbond_imc.md` /
   `tutorial_formulation_smoke.md` / `tutorial_membrane_us.md` / `tutorial_udfcharge.md`
8. **`docs/dependencies.md`** / **`docs/platform_support.md`** — 依存と動作環境
9. **`docs/faq.md`** — よくある詰まりどころ
