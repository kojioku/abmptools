# Dependencies

## Python Version

ABMPTools does not declare a minimum Python version in `setup.py`. Based on language features used (f-strings, `pathlib` usage patterns, `multiprocessing.Pool`), **Python 3.6+** is assumed. Python 3.8+ is recommended for best compatibility.

> **Assumption**: No explicit Python version constraint is declared; the above is inferred from code inspection.

## Required Dependencies

| Package | Purpose | Used In |
|---------|---------|---------|
| **numpy** | Numerical arrays, coordinate transformations, linear algebra, rotations | `molcalc.py`, `abinit_io.py`, `mol_io.py`, `readcif.py`, `cpfmanager.py`, and others |
| **pandas** | DataFrame storage for atom/fragment/dimer data, CSV I/O | `cpfmanager.py`, `getifiepieda.py`, `anlfmo.py`, `generate_difie.py`, `cpf2ifielist.py` |

These are listed in `README.md` but **commented out** in `setup.py` (`install_requires` is not active):

```python
# install_requires=['numpy', 'pandas'],   ← commented out in setup.py
```

Users must ensure numpy and pandas are installed manually or via their environment.

## Optional Dependencies

| Package | Purpose | Used In | Fallback |
|---------|---------|---------|----------|
| **UDFManager** | OCTA COGNAC UDF file I/O | `udf_io.py`, `pdb_io.py` | `try/except` import; UDF features unavailable if missing |
| **gfortran** | Compile Fortran shared library | `Makefile`, `setup.py` (build-time) | Install proceeds; Fortran acceleration unavailable |
| **OpenBabel** (`obabel` CLI) | XYZ → PDB format conversion | `mol_io.py` (via `subprocess`) | Conversion fails if not installed; not required for most workflows |

## Standard Library Usage

ABMPTools makes extensive use of Python standard library modules:

| Module | Purpose |
|--------|---------|
| `argparse` | CLI argument parsing (all entry-point modules) |
| `multiprocessing` | Parallel processing (`Pool`) |
| `os`, `sys` | Path handling, system interaction |
| `re` | Regular expression parsing (log files, PDB) |
| `math` | Trigonometric functions, constants |
| `copy` | Deep copy of data structures |
| `csv` | CSV file writing |
| `gzip` | Compressed file reading (.gz) |
| `glob` | File pattern matching |
| `datetime` | Timestamp generation (CPF headers) |
| `collections` | `OrderedDict`, `defaultdict` |
| `itertools` | Combinatorial iteration |
| `subprocess` | External tool execution (obabel) |
| `statistics` | Mean/stdev calculations |
| `shutil`, `tarfile` | File operations, archive handling |
| `time` | Performance timing |
| `ctypes` | Fortran shared library loading |

## Build-Time Dependencies

### Fortran Compiler (gfortran)

The `Makefile` compiles `readifiepiedalib.f90` into a shared library:

```makefile
# Compilation command (inferred from Makefile):
gfortran -shared -fPIC -o abmptools/f90/bin/readifiepiedalib.so abmptools/f90/src/readifiepiedalib.f90
```

- Triggered automatically during `pip install --user .` via `subprocess.call('make')` in `setup.py`.
- **Not required**: If gfortran is absent, the build step silently fails (`try/except` in `setup.py`), and the package installs without the Fortran library.
- At runtime, `getifiepieda.py` attempts to load the `.so`; if unavailable, it falls back to pure Python. Use `-nof90` flag to explicitly force pure Python mode.

### setuptools

The build system uses `setuptools` with `find_packages()`:

```python
setup(
    name='ABMPTools',
    version='1.14.6',
    packages=find_packages(exclude=('tests', 'docs', 'sample')),
    package_data={'abmptools': ['f90/bin/*', '../*', '../tips/*']}
)
```

> **Note**: There is no `pyproject.toml` or `requirements.txt`. The package relies solely on `setup.py` for build configuration.

## Fortran Shared Library

| File | Size | Purpose |
|------|------|---------|
| `abmptools/f90/src/readifiepiedalib.f90` | 219 lines | Source: fast IFIE/PIEDA log parser |
| `abmptools/f90/bin/readifiepiedalib.so` | ~21 KB | Compiled shared library |

The library is loaded via `ctypes` at runtime and provides significant speedup for parsing large IFIE/PIEDA data from log files.

**Important**: Some features require specific paths:
- MP3/MP4 extraction → requires Fortran module (`readifiepiedalib.so`).
- PB-IFIE, BSSE-IFIE, monomer/dimer energies → requires pure Python mode (`-nof90`).

## Optional Dependencies — abmptools.geomopt

`abmptools.geomopt` はすべての重い依存を実行時まで遅延インポートするため、
未インストールでも `import abmptools` 自体は成功します。
各バックエンドに必要なパッケージは以下の通りです。

### pdbopt バックエンド（`MacePdbOptimizer` / `OpenFFOpenMMMinimizer`）

| パッケージ | バックエンド | 用途 | インストール |
|---|---|---|---|
| `ase` | MACE | 原子シミュレーション環境 | `pip install ase` |
| `mace-torch` | MACE | MACE ML ポテンシャル | `pip install mace-torch` |
| `torch` | MACE | PyTorch（mace-torch の依存） | `pip install torch` |
| `openmm` | OpenFF | 分子動力学エンジン | `conda install -c conda-forge openmm` |
| `openff-toolkit` | OpenFF | OpenFF 力場ツールキット | `conda install -c conda-forge openff-toolkit` |
| `rdkit` | OpenFF | 分子グラフ（fallback） | `conda install -c conda-forge rdkit` |

### qmopt バックエンド（`QMOptimizerPySCF`）

| パッケージ | 必須/任意 | 用途 | インストール |
|---|---|---|---|
| `pyscf` | **必須** | DFT 計算エンジン | `pip install pyscf` |
| `geometric` | 必須（推奨） | geomeTRIC 構造最適化ドライバ | `pip install geometric` |
| `pyberny` | 任意 | berny 構造最適化ドライバ（geometric の代替） | `pip install pyberny`\* |
| `simple-dftd3` | 任意 | D3(BJ) 分散補正（推奨） | `pip install simple-dftd3` |
| `dftd3` | 任意 | D3(BJ) 分散補正（代替） | `pip install dftd3` |

\* pyberny 0.6.3 は setuptools 82+ 環境で `pkg_resources` インポートエラーが発生する既知の問題があります（2026-02 確認）。
`geometric` を推奨します。

分散補正ライブラリがいずれも未インストールの場合は警告を出して
dispersion なしで実行されます（エラーにはなりません）。

### 動作確認済みバージョン（2026-02）

| パッケージ | バージョン |
|---|---|
| pyscf | 2.12.1 |
| geometric | 1.1 |
| numpy | 1.26.4 |

## Dependency Summary

```
ABMPTools
├── [required] numpy
├── [required] pandas
├── [optional] UDFManager       (OCTA COGNAC MD workflows)
├── [optional] gfortran         (build-time, Fortran acceleration)
├── [optional] OpenBabel/obabel (XYZ→PDB conversion)
├── [optional/geomopt-mace]   ase, mace-torch, torch
├── [optional/geomopt-openff] openmm, openff-toolkit, rdkit
├── [optional/geomopt-qm]     pyscf + geometric (or pyberny)
│                             + simple-dftd3 (or dftd3)  [D3 dispersion]
└── [build]    setuptools
```
