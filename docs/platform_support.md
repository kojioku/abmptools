# Platform support matrix (Linux / macOS / Windows / WSL2)

`abmptools` を OS 別にどこまで動かせるかをまとめた早見表。 配布時の前提共有 +
Windows ネイティブ運用したい組織への参照資料を兼ねる。

> **要約**: Linux / macOS は全機能 OK。 Windows ネイティブは **AmberTools が install
> できない** ため `formulation` / `genesis.*` の build pipeline は不可。 WSL2 が
> 使える組織なら Linux と同等。 WSL2 不可の組織には **OpenFF route** (Phase 1
> 開発中) で formulation を Windows native 化する予定 (詳細下記)。

## OS 別 sub-package 対応

| Sub-package | Linux | macOS | **Windows native** | WSL2 |
|---|:---:|:---:|:---:|:---:|
| `abmptools.trajectory` (gmx trjconv/energy wrapper) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.amorphous` (OpenFF route) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.crystal` (OpenFF route) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.cg.peptide` (Martini 3 + vermouth) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.cg.membrane` (insane + cg.peptide) | ✅ | ✅ | ⚠ | ✅ |
| `abmptools.cg.dpd` (Cognac UDF) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.fragmenter` (FMO 自動分割) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.hbond` (J-OCTA UDF colorize) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.gro2udf` (GROMACS → cognac UDF) | ✅ | ✅ | **✅** | ✅ |
| `abmptools.formulation` (Amber route) | ✅ | ✅ | ❌ | ✅ |
| `abmptools.formulation` (**OpenFF route**, Phase 1 開発中) | ✅ | ✅ | **✅ (予定)** | ✅ |
| `abmptools.membrane` (CHARMM/AMBER backend) | ✅ | ✅ | ❌ | ✅ |
| `abmptools.genesis.mmgbsa` (acpype + tleap + atdyn) | ✅ | ✅ | ❌ | ✅ |
| `abmptools.genesis.grest` (tleap + atdyn) | ✅ | ✅ | ❌ | ✅ |
| `abmptools.geomopt` (DFT 経由) | ✅ | ✅ | △ | ✅ |
| `abmptools.fragmenter.cg_segmenter` (CG segment + DPDgen) | ✅ | ✅ | **✅** | ✅ |
| FMO 解析 CLI (`generateajf` / `getifiepieda` / `convertcpf` 等) | ✅ | ✅ | **✅** | ✅ |
| `abinitmp` ビルド (Fortran 90) | ✅ | ✅ | ❌ (gfortran/MSYS2 経由なら可) | ✅ |

凡例: ✅ 動作確認済 / **太字** Windows native で動く重要項目 / ⚠ 一部制約あり /
❌ ネイティブ不可 / △ 設定次第

## Module 別 install 可否

| Module | Linux | macOS | Windows native | WSL2 | 備考 |
|---|:---:|:---:|:---:|:---:|---|
| `numpy` / `pandas` / `scipy` / `matplotlib` | ✅ | ✅ | ✅ | ✅ | pip wheel あり |
| `rdkit` | ✅ | ✅ | ✅ | ✅ | `pip install rdkit` (Windows wheel あり) |
| `parmed` | ✅ | ✅ | ✅ | ✅ | pip OK |
| `MDAnalysis` | ✅ | ✅ | ✅ | ✅ | pip OK |
| `openff-toolkit` / `openff-interchange` | ✅ | ✅ | ✅ | ✅ | conda-forge / pip 両対応 |
| `openff-amber-ff-ports` (ff14SB SMIRNOFF) | ✅ | ✅ | ✅ | ✅ | pip OK、 Windows route の鍵 |
| `vermouth` (Martini martinize2) | ✅ | ✅ | ✅ | ✅ | pip OK |
| `pypdf` / `pdfminer` | ✅ | ✅ | ✅ | ✅ | pip OK |
| **`packmol`** | ✅ | ✅ | ✅ | ✅ | 公式 binary (Windows .exe あり) |
| **`gmx` (GROMACS)** | ✅ | ✅ | ✅ (CPU) / △ (GPU) | ✅ | Windows native installer あり、 ただし CUDA GPU は WSL2 推奨 |
| **`tleap` / `antechamber` / `sqm` (AmberTools)** | ✅ | ✅ | **❌** | ✅ | 公式は Linux/macOS のみ、 conda-forge にも Windows パッケージ無し |
| **`acpype`** | ✅ | ✅ | **❌** | ✅ | Python パッケージだが `antechamber` + `sqm` 依存 |
| `martinize2` | ✅ | ✅ | ✅ | ✅ | vermouth 経由で OK |
| `insane` (CG membrane builder) | ✅ | ✅ | ⚠ | ✅ | Python 2/3 互換 script、 Windows パス区切り注意 |
| `xtb` (semiempirical QM) | ✅ | ✅ | ✅ | ✅ | Windows binary あり、 OpenFF route の AM1-BCC 代替候補 |
| GROMACS の CUDA GPU offload | ✅ | △ | △ | ✅ | Windows native でも build できるが install 手間大、 WSL2 推奨 |

## 推奨 setup 別シナリオ

### A. Linux / macOS / WSL2 ユーザー — 全機能 OK

```bash
micromamba create -n abmptoolsenv -c conda-forge ambertools rdkit parmed openff-toolkit \
                   openff-interchange openff-amber-ff-ports vermouth
micromamba activate abmptoolsenv
pip install -e <abmptools repo>
```

### B. Windows native ユーザー (WSL2 不可組織) — Windows route で formulation 可

```powershell
# Anaconda / Miniforge を install
conda create -n abmptoolsenv -c conda-forge rdkit parmed openff-toolkit openff-interchange ^
                                            openff-amber-ff-ports vermouth packmol
conda activate abmptoolsenv
pip install abmptools

# GROMACS は公式 Windows installer (CPU run; GPU は限定)
# https://manual.gromacs.org/current/install-guide/index.html
```

`abmptools.formulation` の `config.json` で:

```json
"force_field_route": "openff"
```

を指定すると、 AmberTools 経由を回避して全 stage が Windows native で動作する
(**Phase 1 開発中**、 進捗は CHANGELOG 参照)。

### C. Windows native ユーザー (build は別マシン) — post-process / 可視化のみ

build + mdrun は WSL2 or Linux マシンで実行し、 `prod.xtc / prod.tpr / *.xvg` を
Windows にコピーして:

```powershell
pip install abmptools
python -m abmptools.trajectory thin_nojump --traj prod\prod.xtc --tpr prod\prod.tpr --skip 10
python -m abmptools.trajectory wrap_pbc --traj prod\prod.xtc --tpr prod\prod.tpr
```

→ VMD on Windows / OCTA Viewer for Windows で可視化。

## `abmptools.trajectory` の Windows 互換設計

今回 (v2.0.0+) 新規追加した `abmptools.trajectory` は全 sub-package で **唯一**
Windows native 動作を明示的に保証している module。 設計選択:

| 選択 | 理由 |
|---|---|
| `subprocess.run(shell=False, input=stdin.encode())` | `echo \| gmx ...` パイプを排除、 Windows cmd でも動作 |
| `pathlib.Path` で全 path 操作 | `/` `\` 差を自動吸収 |
| `shutil.which(gmx)` で実行 path 解決 | PATH 解決を OS-agnostic に |
| bash 依存ゼロ | Windows native cmd / PowerShell でも実行可 |

旧 sample `wrap_pbc.sh` / `gen_for_udf.sh` は **deprecated**、 生成物は **`wrap_pbc.py` / `gen_for_udf.py`** に統一済 (`amorphous.mdp_protocol.write_wrap_script` / `write_udf_export_script` が Python script を出力)。

## Phase 1 計画: `formulation` の Windows route

### 現状 (Amber route)

| Stage | tool | OS 制約 |
|---|---|---|
| 1 peptide_atomistic | `tleap` (AmberTools) | Linux/macOS のみ |
| 2 small_mol_parameterize | `acpype` + `sqm` (AmberTools) | 同上 |
| 3 packmol | `packmol` | OK |
| 4 solvate_ions | `tleap solvatebox + addions` + `parmed.to_gromacs()` | Linux/macOS のみ |
| 5-7 index / mdp / run_script | Python 内 | OK |

### 新 OpenFF route (Phase 1 実装予定)

| Stage | tool | Windows OK? |
|---|---|---|
| 1 peptide_atomistic | **OpenFF Toolkit + `openff-amber-ff-ports` (ff14SB SMIRNOFF)** | ✅ |
| 2 small_mol_parameterize | **OpenFF Sage 2.x (SMIRNOFF) + Interchange** | ✅ |
| 2.5 charges | OpenFF `am1bcc` (内部 OE/`xtb`/事前計算) | ✅ (xtb 経路) |
| 3 packmol | 共通 | ✅ |
| 4 solvate_ions | `Interchange.to_gromacs()` + Joung-Cheatham ions | ✅ |
| 5-7 | 共通 | ✅ |

### 設定切替

`config.json` に新フィールド:

```json
"force_field_route": "amber"    // 現行 (default、 Linux/macOS)
"force_field_route": "openff"   // 新規 (全 OS)
```

### 既知の課題と対応案

| 課題 | 対応 |
|---|---|
| OpenFF の `am1bcc` charges 計算で `sqm` が要るケース | (a) `xtb` で代替 (Windows native binary)、 または (b) 事前計算 charges を `.mol2` で渡す |
| Cys2-Cys7 disulfide bond declare | OpenFF Topology に `top.add_bond(...)` で手動追加 (amorphous でも同パターン) |
| 非標準残基 (D-Phe1 / D-Trp4 / Thr-ol) | `parameterize_method = "gaff"` (whole-peptide SMIRNOFF) で対応 (現行 D 体経路と同じ) |
| `openff-amber-ff-ports` の peptide サポート範囲 | natural L-amino + ACE/NME cap + disulfide が基本サポート。 ff14SB 全機能ではない |

### 工数

| 作業 | 規模 |
|---|---|
| `formulation/peptide_atomistic_openff.py` 新規 | ~150 行 |
| `formulation/small_molecule_openff.py` 新規 | ~100 行 |
| `formulation/topology.py` に Interchange merge stage 追加 | ~80 行 |
| builder に route 分岐 + tests | ~50 行 |
| Windows smoke test (CI 含む) | ~100 行 |
| **合計** | **~480 行 + tests** |

amorphous で確立した OpenFF パターンを再利用するので、 reusable code 多い。

## 関連 docs

- [`docs/dependencies.md`](dependencies.md) — sub-package 別の依存リスト
- [`docs/amorphous_tutorial_windows.md`](amorphous_tutorial_windows.md) — amorphous の Windows native tutorial (既存)
- [`docs/formulation.md`](formulation.md) — formulation pipeline 詳細
- [`docs/licenses_third_party.md`](licenses_third_party.md) — 各依存ツールの license 一覧

## 変更履歴

- **2026-06-06**: 初版作成。 v2.0.0 時点の対応状況 + Phase 1 (`formulation` Windows route) 計画
