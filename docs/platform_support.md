# Platform support matrix (Linux / macOS / Windows / WSL2)

`abmptools` を OS 別にどこまで動かせるかをまとめた早見表。 配布時の前提共有 +
Windows ネイティブ運用したい組織への参照資料を兼ねる。

> **要約**: Linux / macOS は全機能 OK。 Windows ネイティブは **AmberTools が install
> できない** ため `genesis.*` の build pipeline は不可。 WSL2 が使える組織なら
> Linux と同等。 `formulation` は **OpenFF route (`force_field_route="openff"`)
> で Windows native 化完了** — multi-chain protein + disulfide (insulin) を含め
> tleap/acpype 不要で全 OS build 可能 (詳細下記 "OpenFF route" 節)。

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
| `abmptools.formulation` (Amber route, tleap) | ✅ | ✅ | ❌ | ✅ |
| `abmptools.formulation` (**OpenFF route**, `force_field_route="openff"`) | ✅ | ✅ | **✅** | ✅ |
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

### B. Windows native ユーザー (WSL2 不可組織) — Windows route で formulation 可

```powershell
# Anaconda / Miniforge を install
# 注: Windows では `openff-toolkit` (メタパッケージ) は AmberTools に hard-depend し
#     conda solve が "ambertools does not exist" で失敗する。 ambertools 非依存
#     (RDKit backend) の `openff-toolkit-base` を使う。 ff14SB library charges のみ
#     使う formulation OpenFF route はこれで十分 (AM1-BCC=sqm は呼ばない)。
conda create -n abmptoolsenv -c conda-forge python=3.11 rdkit parmed ^
             openff-toolkit-base openff-interchange openff-amber-ff-ports ^
             openff-nagl openff-nagl-models pdbfixer openmm vermouth packmol
conda activate abmptoolsenv
pip install abmptools PeptideBuilder biopython

# GROMACS は公式 Windows installer (CPU run; GPU は限定)
# https://manual.gromacs.org/current/install-guide/index.html
```

> この Windows recipe (`openff-toolkit-base` + `openff-nagl` 経路) は
> `.github/workflows/windows-native.yml` の `windows-openff-smoke` ジョブが
> `windows-latest` runner で実機検証している (2026-06-16 green):
> - protein: sequence → PeptideBuilder → PDBFixer → `Topology.from_pdb` → ff14SB
>   library charges
> - small mol: SMILES → **NAGL ML-AM1-BCC 電荷** → Sage Interchange (sqm 不要)
>
> **電荷の扱い (Windows で AmberTools が無いことの帰結)**: protein は ff14SB の
> library charges なので電荷計算が要らない。 小分子 (caprate / taurocholate 等) は
> AM1-BCC 相当が要るが、 `sqm` (AmberTools) は Windows に無いため **`openff-nagl`
> (graph neural net の ML 電荷、 total charge を formal charge に保存) で代替**する。
> `config.json` の small-molecule charge_method は `"nagl"` を指定 (Linux/macOS で
> 厳密な AM1-BCC が要るなら `"am1bcc"`、 軽量 fallback は `"gasteiger"`)。

`abmptools.formulation` の `config.json` で:

```json
"force_field_route": "openff"
```

を指定すると、 AmberTools 経由を回避して全 stage が Windows native で動作する
(**実装完了**、 insulin 含め実証済。 詳細は下記 "OpenFF route (実装済)" 節)。

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

## OpenFF route (実装済): `formulation` の Windows native build

`force_field_route="openff"` で AmberTools (tleap/acpype/sqm) を一切使わず
全 OS で build できる。 multi-chain protein + disulfide (insulin) を含め実証済。

### route 別の stage (Amber=Linux/macOS、 OpenFF=全 OS)

| Stage | Amber route (tleap) | **OpenFF route (全 OS)** |
|---|---|---|
| 1 peptide | `tleap` (ff14SB) | **PDBFixer (water除去+H付加) → `Topology.from_pdb`** (multi-chain 1 分子認識 + disulfide 自動検出) → ff14SB SMIRNOFF (`ff14sb_off_impropers_0.0.4.offxml`) library charges。 sequence は PeptideBuilder で 3D 生成 |
| 2 small mol | `acpype` + `sqm` (AM1-BCC) | **OpenFF Sage 2.x SMIRNOFF** + RDKit ETKDGv3、 charge は gasteiger/nagl precomputed |
| 3 packmol | 共通 | 共通 (water 1 個を typing template に含める) |
| 4 topology | `tleap solvatebox+addions` + parmed | **`build_protein_route_topology`**: 単一コピー Interchange → `to_top` → `[molecules]` count 複製 → `gmx editconf` で .gro。 wet 化は `gmx solvate`+`genion` (Joung-Cheatham 0.15 M) |
| 5-7 | 共通 | 共通 |

### 設定切替

```json
"force_field_route": "amber"    // default (Linux/macOS)
"force_field_route": "openff"   // 全 OS (Windows native)
```

### 実装上の重要ポイント (insulin で確立)

| 課題 | 対応 |
|---|---|
| `Molecule.from_polymer_pdb` が multi-chain+disulfide 非対応 | **`Topology.from_pdb`** (OpenFF 0.14+) — 1 分子認識 + S-S 自動検出 (手動宣言不要) |
| OpenFF は explicit H 必須、 crystal water が邪魔 | **PDBFixer** `removeHeterogens(keepWater=False)` + `addMissingHydrogens(7.0)` |
| **full mixture の `Interchange.from_smirnoff` が O(N²) で 2 時間+** | **単一コピー parametrize (~4 秒) → `[molecules]` count 複製** (GROMACS の moleculetype count 参照の本来の使い方) |
| `am1bcc` は `sqm` (Linux 専用) | protein=ff14SB **library charges**、 small mol=gasteiger/nagl precomputed (`protein_flags` で振り分け) |
| `gmx genion` が ion moleculetype を要求 | Na+/Cl- を単一 topology に含めて moleculetype 定義のみ生成、 `[molecules]` 除外 |
| 非標準残基 (D-Phe1/D-Trp4/Thr-ol) | Amber route の whole-peptide GAFF (`parameterize_method="gaff"`) で対応 (D-octreotide)。 OpenFF route は標準 protein のみ |
| >99999 atom 系の gro serial wrap | `ndx.parse_gro_residues` を 1-based 行番号に修正済 (insulin 156k で発覚) |

### 実証済 (実機 build)

| 系 | atoms | route | 結果 |
|---|---|---|---|
| insulin 2G4M × 6 | 168,899 (wet) | openff | S-S 3 本/copy 自動検出、 grompp em/nvt PASS |
| kggggg × 2 | 20,916 (wet) | openff | sequence build (PeptideBuilder) |
| octreotide D × 6 | 90,494 | amber GAFF | 実薬構造 (DPN1/DTR4)、 100 ns 完走 |

### 残課題 (Phase 2-D)

- GitHub Actions `windows-latest` runner での実機 smoke (現状は Linux/WSL2 で
  OS-agnostic 実装を確認、 Windows 実機 end-to-end は未検証)
- GROMACS GPU は Windows native だと制約あり (CPU run は可、 GPU は WSL2 推奨)

## 関連 docs

- [`docs/dependencies.md`](dependencies.md) — sub-package 別の依存リスト
- [`docs/amorphous_tutorial_windows.md`](amorphous_tutorial_windows.md) — amorphous の Windows native tutorial (既存)
- [`docs/formulation.md`](formulation.md) — formulation pipeline 詳細
- [`docs/licenses_third_party.md`](licenses_third_party.md) — 各依存ツールの license 一覧

## 変更履歴

- **2026-06-06**: 初版作成。 v2.0.0 時点の対応状況 + Phase 1 (`formulation` Windows route) 計画
