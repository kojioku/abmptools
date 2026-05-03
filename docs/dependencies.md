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

## Optional Dependencies — abmptools.amorphous

`abmptools.amorphous` もすべての重い依存を実行時まで遅延インポートするため、
未インストールでも `import abmptools` は成功します。amorphous ビルダーを動かすには以下が必要です。

### 必須ランタイム

| パッケージ | 用途 | インストール |
|---|---|---|
| `openff-toolkit` | 分子準備・力場適用 (SMILES / SDF) | `conda install -c conda-forge openff-toolkit` |
| `openff-interchange` | GROMACS 形式 (gro/top) への書き出し | `conda install -c conda-forge openff-interchange` |
| `openmm` | バックエンド (Interchange が内部で使用) | `conda install -c conda-forge openmm` |
| `rdkit` | SDF 読み込みと化学的サニタイズ | `conda install -c conda-forge rdkit` |
| `packmol` | 初期配置 (外部バイナリ、CLI で呼ばれる) | `conda install -c conda-forge packmol` |

### 電荷バックエンド (どちらか一方)

| パッケージ | 特徴 |
|---|---|
| `ambertools` (AM1-BCC、**推奨**) | `conda install -c conda-forge ambertools` |
| `openff-nagl` (ML 電荷、AmberTools 不要) | `pip install openff-nagl` |

### 既知の制約: `setuptools<81`

`openff.amber_ff_ports` は `pkg_resources` を `import` するため、
setuptools 82+ (2025 以降デフォルト) では `ModuleNotFoundError: No module named 'pkg_resources'` が発生します。
対処:

```bash
pip install "setuptools<81"
```

(この注意は geomopt の `pyberny` と同じく、将来的に上流が対応したら不要になります)

### 外部ツール (任意、後処理)

| ツール | 用途 |
|---|---|
| `gromacs` (2020+ 推奨、2025.4 で動作確認) | `md/run_all.sh` の実行 + `wrap_pbc.sh` の `gmx trjconv` |
| `vmd` | 生成された `*_pbc.xtc` / `05_npt_final_pbc.gro` の可視化 |

### 追加依存 (`abmptools.amorphous.pubchem`)

PubChem から 3D SDF や SMILES を取得する補助モジュール。
Python 追加依存は **なし** (`urllib` 標準ライブラリのみ) ですが、使用時には
`https://pubchem.ncbi.nlm.nih.gov` への HTTPS アクセスが必要です。オフライン
環境では事前に SDF をダウンロードしておくか `--mol` で指定してください。

## Optional Dependencies — abmptools.membrane

`abmptools.membrane` (peptide-bilayer Umbrella Sampling builder) も重い依存を
遅延インポートするため、未インストールでも `import abmptools` は成功します。
ビルダーを動かすには以下が必要です。

### 必須ランタイム

| パッケージ | 用途 | インストール |
|---|---|---|
| `ambertools` (≥ 22) | `tleap` (ペプチド構築 + 力場割り当て)、`packmol-memgen` (bilayer 構築)、`antechamber` (新規小分子の GAFF2 パラメータ化) | `conda install -c conda-forge 'ambertools>=22'` |
| `parmed` | AMBER `prmtop`/`inpcrd` → GROMACS `top`/`gro` 変換 | `conda install -c conda-forge parmed` |
| `gromacs` (2021+) | `gmx grompp` / `gmx mdrun` / `gmx wham` / `gmx trjconv` | `conda install -c conda-forge 'gromacs[build=nompi*]'` |

### 推奨環境

```bash
# 最小構成 (CPU のみ、AMBER backend)
micromamba create -n abmptoolsenv -c conda-forge \
    python=3.10 gromacs=2021 'ambertools>=22' parmed openff-toolkit numpy
```

### 任意: GPU 加速 (NVIDIA)

`bioconda::gromacs=2021.3` は OpenCL ビルドで、WSL2 / Linux で NVIDIA GPU を
直接認識できないことがあります。CUDA ビルドの GROMACS だけを別 env に
インストールすることで、`run.sh` 経由の MD だけを GPU offload できます:

```bash
micromamba create -n gmxcudaenv -c conda-forge 'gromacs[build=nompi_cuda*]' -y
# run.sh 起動時に GMX env var を CUDA env の gmx に向け、
# MDRUN_OPTS="-nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu -pin on" を渡す
```

詳細は [`tutorial_membrane_us.md`](tutorial_membrane_us.md) §1.3 / §3.2、
[`membrane.md`](membrane.md) の "GPU 加速 (NVIDIA / CUDA)" セクション参照。

### 必須一回限り: packmol-memgen NumPy ≥ 1.24 互換性パッチ

`packmol-memgen` 2023.2.24 は内部の `pdbremix/v3numpy.py` で削除済の
`np.float` を参照するため、NumPy 1.24+ 環境では下記の 2 行 sed が必要:

```bash
ENV=~/.local/share/mamba/envs/abmptoolsenv
LIB=$ENV/lib/python3.10/site-packages/packmol_memgen/lib/pdbremix
sed -i 's/np\.zeros(3, dtype=np\.float)/np.zeros(3, dtype=float)/' $LIB/v3numpy.py
sed -i 's/np\.array(args, dtype=np\.float, copy=True)/np.array(args, dtype=float, copy=True)/' $LIB/v3numpy.py
```

冪等。再実行しても問題なし。

### CHARMM36 backend を使う場合の追加依存

`MembraneConfig.backend="charmm36"` で動かすには、Klauda lab の **GROMACS port**
(`charmm36-jul2022.ff/` ディレクトリ) を MacKerell 研の公式ページから別途
ダウンロードして配置し、`MembraneConfig.charmm_ff_dir` で絶対パスを指定する
必要があります。本パッケージは ff ディレクトリを同梱しません (license 配布元の
差し替えを許容するため)。

取得手順は [`membrane.md`](membrane.md#charmm36-gromacs-port-の取得-phase-c-用) 参照。
**CGenFF Web server / CHARMM-GUI の利用は禁止されている** (本パッケージの
商用利用 OK ルール) ため、新規小分子のパラメータ化は AMBER 経路 (`antechamber`
+ GAFF2) を使う。

### 任意: PMF 解析の代替

| パッケージ | 用途 |
|---|---|
| `pymbar` | WHAM の代替 / 検証用 MBAR 解析 (現在 stub、将来本実装予定) |

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
├── [optional/amorphous]      openff-toolkit, openff-interchange, openmm,
│                             rdkit, packmol (binary),
│                             ambertools (AM1-BCC) or openff-nagl (ML)
│                             + gromacs (for run_all.sh / wrap_pbc.sh)
├── [optional/membrane]       ambertools (tleap + packmol-memgen + antechamber),
│                             parmed, gromacs (2021+),
│                             + CHARMM36 .ff dir (Klauda port, manual DL)
│                             + optional: gromacs[cuda] for GPU offload
├── [optional/geomopt-mace]   ase, mace-torch, torch
├── [optional/geomopt-openff] openmm, openff-toolkit, rdkit
├── [optional/geomopt-qm]     pyscf + geometric (or pyberny)
│                             + simple-dftd3 (or dftd3)  [D3 dispersion]
└── [build]    setuptools (<81 if using openff.amber_ff_ports)
```
