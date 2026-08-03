# Platform support matrix (Linux / macOS / Windows / WSL2)

`abmptools` を OS 別にどこまで動かせるかをまとめた早見表。 配布時の前提共有 +
Windows ネイティブ運用したい組織への参照資料を兼ねる。

> **要約**: Linux / macOS は全機能 OK。 Windows ネイティブは **AmberTools が install
> Linux と同等。 OpenFF ベースの経路
> で Windows native 化完了** — multi-chain protein + disulfide (insulin) を含め
> tleap/acpype 不要で全 OS build 可能 (詳細下記 "OpenFF route" 節)。

## Module 別 install 可否

| Module | Linux | macOS | Windows native | WSL2 | 備考 |
|---|:---:|:---:|:---:|:---:|---|
| `numpy` / `pandas` / `scipy` / `matplotlib` | ✅ | ✅ | ✅ | ✅ | pip wheel あり |
| `rdkit` | ✅ | ✅ | ✅ | ✅ | `pip install rdkit` (Windows wheel あり) |
| `parmed` | ✅ | ✅ | ✅ | ✅ | pip OK |
| `MDAnalysis` | ✅ | ✅ | ✅ | ✅ | pip OK |
| `openff-toolkit` / `openff-interchange` | ✅ | ✅ | ✅ | ✅ | conda-forge / pip 両対応。 **Windows conda は `openff-toolkit-base` を使う** (メタパッケージ `openff-toolkit` は AmberTools に hard-depend し Windows で solve 不能。 base は RDKit backend で ambertools 非依存) |
| `openff-amber-ff-ports` (ff14SB SMIRNOFF) | ✅ | ✅ | ✅ | ✅ | pip OK、 Windows route の鍵 |
| `openff-nagl` / `openff-nagl-models` | ✅ | ✅ | ✅ | ✅ | 小分子の ML-AM1-BCC 電荷 (graph NN、 pure-Python)。 **Windows で sqm/AmberTools の代替**。 torch を引く |
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

### B. Windows native ユーザー (build は別マシン) — post-process / 可視化のみ

build + mdrun は WSL2 or Linux マシンで実行し、 `prod.xtc / prod.tpr / *.xvg` を
Windows にコピーして:

```powershell
pip install abmptools
python -m abmptools.trajectory thin_nojump --traj prod\prod.xtc --tpr prod\prod.tpr --skip 10
python -m abmptools.trajectory wrap_pbc --traj prod\prod.xtc --tpr prod\prod.tpr
```

→ VMD on Windows / OCTA Viewer for Windows で可視化。

## `abmptools.trajectory` の Windows 互換設計

`abmptools.trajectory` (v2.0.0+) は全 sub-package で **唯一**
Windows native 動作を明示的に保証している module。 設計選択:

| 選択 | 理由 |
|---|---|
| `subprocess.run(shell=False, input=stdin.encode())` | `echo \| gmx ...` パイプを排除、 Windows cmd でも動作 |
| `pathlib.Path` で全 path 操作 | `/` `\` 差を自動吸収 |
| `shutil.which(gmx)` で実行 path 解決 | PATH 解決を OS-agnostic に |
| bash 依存ゼロ | Windows native cmd / PowerShell でも実行可 |

旧 sample `wrap_pbc.sh` / `gen_for_udf.sh` は **deprecated**、 生成物は **`wrap_pbc.py` / `gen_for_udf.py`** に統一済 (`amorphous.mdp_protocol.write_wrap_script` / `write_udf_export_script` が Python script を出力)。

## 関連 docs

- [`docs/dependencies.md`](dependencies.md) — sub-package 別の依存リスト
- [`docs/amorphous_tutorial_windows.md`](amorphous_tutorial_windows.md) — amorphous の Windows native tutorial (既存)
- [`docs/licenses_third_party.md`](licenses_third_party.md) — 各依存ツールの license 一覧

## 変更履歴

- **2026-06-06**: 初版作成。 v2.0.0 時点の対応状況
- **2026-08-02**: v2.9.0 で開発中手法のサブパッケージを分離。 対応表は公開分のみを載せる
- **2026-08-03**: v2.10.0 / v2.11.0 でさらに分離。 Windows route の詳細は移設先へ
