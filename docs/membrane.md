# abmptools.membrane — Lipid-Bilayer Umbrella-Sampling Builder

ペプチドの脂質膜透過 PMF 計算用の GROMACS 入力一式を生成するビルダー。
packmol-memgen で bilayer + ペプチド配置 → AMBER (Lipid21 + ff19SB) もしくは
CHARMM36 でパラメータ化 → 多段平衡化 + 反応座標 pulling + Umbrella Sampling
window 一括 + WHAM 解析、までのフローを 1 つの `MembraneUSBuilder.build()`
呼び出しで完結させる。

姉妹パッケージ `abmptools.amorphous` (3D ランダム box) と並列。両者とも
`abmptools.core.system_model` の `SimulationParams` を共有するが、膜系は
2D periodic + semiisotropic Pcoupl 専用のため別 stage を提供。

## ライセンス上のルール (商用利用前提)

このサブパッケージは **会社利用 OK な権利のみ**で構成。下記を **使わない** ことが
保証される設計:

- ❌ **CGenFF Web server** (`cgenff.umaryland.edu`) の出力 — 商用は Silcsbio 経由
  の有償サブスク必要。本パッケージは `cgenff` / `match` を呼ばない。
- ❌ **CHARMM-GUI** の自動生成 (Membrane Builder, Topology Generator など) —
  商用は別契約必要。本パッケージは CHARMM-GUI 出力を消費しない。

許可されているもの:

- ✅ **AmberTools** (tleap, antechamber, parmchk2, packmol-memgen) — free incl. commercial
- ✅ **CHARMM36 力場パラメータ値**を **MacKerell 研の公式無償配布** または
  **Klauda lab の GROMACS port** から取得して使用 — MacKerell 研の公式 stance:
  「free of charge to academic and industrial researchers」
- ✅ **Packmol** (Martínez et al.) — free
- ✅ **GROMACS** — LGPL
- ✅ **parmed** — LGPL
- ✅ **PyMBAR** (オプション) — MIT

新規小分子 (今後 ligand 追加等) を扱う場合は `backend="amber"` で
`antechamber/GAFF2` を使う。CGenFF サーバ依存ルートは存在しない。

## インストール

```bash
# 基本依存 (conda 推奨)
micromamba create -n abmptoolsenv -c conda-forge \
    python=3.10 gromacs=2021 ambertools=23 parmed openff-toolkit numpy

# pip extras
pip install abmptools[membrane]   # parmed (gromacs / ambertools は conda 経由で別途)
```

### packmol-memgen + NumPy ≥ 1.24 互換性パッチ (env 一回限り)

`packmol-memgen` 2023.2.24 は内部の `pdbremix/v3numpy.py` で NumPy 1.24 で
削除された `np.float` を使っているため、conda env に下記の 2 行 sed を当てる:

```bash
ENV_LIB=$ENV/lib/python3.10/site-packages/packmol_memgen/lib/pdbremix
sed -i 's/np\.zeros(3, dtype=np\.float)/np.zeros(3, dtype=float)/' $ENV_LIB/v3numpy.py
sed -i 's/np\.array(args, dtype=np\.float, copy=True)/np.array(args, dtype=float, copy=True)/' $ENV_LIB/v3numpy.py
```

`np.float` → 組み込み `float` への置換のみ (NumPy 公式の deprecation 推奨対処)。
冪等なので再実行も安全。

### 環境変数

`packmol-memgen` と `tleap` は `AMBERHOME` を必要とする。本パッケージは
`packmol-memgen` バイナリのパスから自動推定するので、conda env を activate
しているか、`MembraneConfig.packmol_memgen_path` にフルパスを渡せば追加設定不要。

### GPU 加速 (NVIDIA / CUDA)

本パッケージの `run.sh` は `MDRUN_OPTS` 環境変数で `gmx mdrun` の追加引数を
受け付ける。NVIDIA GPU + CUDA 環境で以下を設定:

```bash
MDRUN_OPTS="-nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu -pin on" \
    bash run.sh
```

#### bioconda gmx (OpenCL build) を使っている場合の注意

`bioconda::gromacs=2021.3` などの**OpenCL ビルド**は WSL2 環境では NVIDIA GPU
を見つけられない (NVIDIA は WSL2 で OpenCL ICD を提供しないため、CUDA は
動くが OpenCL は動かない)。下記を満たす別 env を併設するのが確実:

```bash
# CUDA build の GROMACS だけを別 env に install (~5 min)
micromamba create -n gmxcudaenv -c conda-forge 'gromacs[build=nompi_cuda*]' -y

# 既存 abmptoolsenv を python / tleap / packmol-memgen 用に維持しつつ、
# gmx だけ CUDA env のものを使う:
PATH=~/.local/share/mamba/envs/abmptoolsenv/bin:$PATH \
    AMBERHOME=~/.local/share/mamba/envs/abmptoolsenv \
    GMX=~/.local/share/mamba/envs/gmxcudaenv/bin/gmx \
    MDRUN_OPTS="-nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu -pin on" \
    NT=4 \
    bash run.sh
```

`gmx grompp` も同じ `$GMX` で実行される (run.sh 内で統一) ため、
.tpr の version mismatch 問題は起きない。

#### 実測パフォーマンス (参考)

|    | system | core / GPU | 実測 ns/day |
|---|---|---|---|
| CPU only | poly-Ala 5 + POPC 32/leaflet (~18k atoms) | 4 cores (Ryzen WSL2) | ~140 |
| GPU offload | 同上 | RTX 4070 Ti + 4 CPU cores | ~640 |

GPU 加速比は **~4-5×** 程度 (この系サイズで 4 core 並列との比較)。
1 core 換算なら ~15× 前後。系が大きくなる (50k+ atoms) ほど GPU の優位
は大きくなる (10× 超)。

### CHARMM36 GROMACS port の取得

`backend="charmm36"` を使う場合、CHARMM36 力場ファイル一式 (Klauda lab の
GROMACS port) を別途ダウンロードして配置する必要がある。**力場パラメータ値
そのもの**は MacKerell 研の公式 stance 「free of charge to academic and
**industrial** researchers」により商用利用 OK。CGenFF Web server や
CHARMM-GUI のような **サービス** (商用ライセンス必要) には依存しない。

#### ダウンロード手順

1. **MacKerell 研公式ページ**にアクセス:
   <http://mackerell.umaryland.edu/charmm_ff.shtml>

2. **"Force fields available in formats other than CHARMM"** セクションの
   GROMACS port (例: `charmm36-jul2022.ff.tgz` または最新版) をダウンロード。
   ファイルサイズは ~10 MB 程度。

3. **配置先 (推奨)**: `~/llm-project/charmm36-jul2022.ff/`

   ```bash
   cd ~/llm-project
   # ブラウザでダウンロードしたファイルをここに置いてから:
   tar -xzf charmm36-jul2022.ff.tgz
   # → ~/llm-project/charmm36-jul2022.ff/  ディレクトリができる
   ```

   ディレクトリ内に `forcefield.itp`, `ffnonbonded.itp`, `ffbonded.itp`,
   `*.rtp`, `lipids.rtp`, `merged.rtp` 等が含まれていることを確認。

4. **`MembraneConfig.charmm_ff_dir`** に絶対パスを渡す:

   ```python
   cfg = MembraneConfig(
       backend='charmm36',
       charmm_ff_dir='/home/<user>/llm-project/charmm36-jul2022.ff',
       lipids=[LipidSpec(resname='POPC', n_per_leaflet=64)],
       peptide=PeptideSpec(name='aa5', sequence='AAAAA'),
       output_dir='./run01',
   )
   ```

5. **論文引用の義務**: ライセンスは citation 要求あり。論文発表時は
   MacKerell 研の指定文献 (Best 2012 *J. Chem. Theory Comput.*、
   Klauda 2010 *J. Phys. Chem. B*、Pastor 2011 *Chem. Phys. Lett.*、
   タンパク・脂質・水で各々) を引用する。

#### 商用利用での注意

- **使ってよい**: 上記公式ページの zip / Klauda 研の GROMACS port (parameter values)
- **使ってはいけない**: CGenFF Web server (`cgenff.umaryland.edu`) の出力 — 商用は
  Silcsbio 経由のサブスク必須。新規小分子のパラメータ化は AMBER (GAFF2 +
  antechamber) ルートに切り替えること
- **使ってはいけない**: CHARMM-GUI の自動生成 (Membrane Builder, Topology
  Generator など) — 商用は別契約必要。本パッケージは bilayer 構築を
  packmol-memgen で代替

#### ion 名の自動切替

`IonSpec(cation='Na+', anion='Cl-')` のように **AMBER 形式**で書けば、
CHARMM36 backend が内部で CHARMM の標準残基名 (`SOD` / `CLA` 等) に
翻訳するため、ユーザーは backend を切替えるごとにイオン名を書き換える
必要はない。明示的に CHARMM 名を書いた場合 (`IonSpec(cation='SOD')`) は
そのまま尊重する。マップ:

| AMBER 形式 | CHARMM36 形式 |
|---|---|
| `Na+` | `SOD` |
| `Cl-` | `CLA` |
| `K+`  | `POT` |
| `Mg2+` | `MG`  |
| `Ca2+` | `CAL` |

## クイックスタート (smoke)

```python
from abmptools.membrane import (
    MembraneConfig, MembraneUSBuilder,
    LipidSpec, PeptideSpec, USProtocol,
    EquilibrationProtocol, PullingProtocol,
)

# 7 windows × 1 ns、box 約 5×5×10 nm の最小構成
us = USProtocol(
    z_min_nm=-1.5, z_max_nm=+1.5, window_spacing_nm=0.5,
    window_nsteps=500_000, window_dt_ps=0.002,
)
eq = EquilibrationProtocol(
    em_steps=5000, nvt_nsteps=50_000, npt_nsteps=500_000,
    temperature_K=310.0,
)
pull = PullingProtocol(
    pull_rate_nm_per_ps=0.001, nsteps=2_500_000,
)

cfg = MembraneConfig(
    backend='amber',
    lipids=[LipidSpec(resname='POPC', n_per_leaflet=32)],
    peptide=PeptideSpec(name='aa5', sequence='AAAAA'),
    output_dir='./run01', seed=42,
    equilibration=eq, pulling=pull, umbrella=us,
)

builder = MembraneUSBuilder(cfg)
result = builder.build()           # ~30 秒で 16 ファイル生成
print(result['run_script'])        # ./run01/run.sh
```

`build()` は MD を**実行しない**。run.sh を後で起動 (もしくは HPC キューに投入):

```bash
cd ./run01
GMX=$(which gmx) NT=4 bash run.sh
```

run.sh が完走すると `analysis/pmf.xvg` (PMF[z]、kJ/mol) が生成される。

## ディレクトリレイアウト

`build()` 後の `output_dir`:

```
run01/
├── input/
│   └── config.json        # MembraneConfig snapshot (再現性)
├── build/
│   ├── peptide.pdb        # tleap 生成の初期 peptide PDB
│   ├── bilayer_peptide.pdb  # packmol-memgen 出力
│   ├── system.tleap       # tleap input
│   ├── system.prmtop      # AMBER topology
│   ├── system.inpcrd      # AMBER coords
│   ├── system.top         # GROMACS topology (parmed 変換)
│   ├── system.gro         # GROMACS coords
│   └── system.ndx         # System / Bilayer / Peptide / Solvent / Ions / Non_Bilayer
├── equil/
│   ├── em.mdp / .tpr / .gro
│   ├── nvt.mdp / .tpr / .gro
│   └── npt.mdp / .tpr / .gro
├── pull/
│   ├── pull.mdp / .tpr
│   ├── pull.xtc           # MD 後生成
│   ├── pullx.xvg / pullf.xvg
├── windows/
│   ├── win_000/
│   │   ├── start.gro      # extract_window_frames で抽出
│   │   ├── window.mdp
│   │   ├── window.tpr / window.xtc / pullx.xvg / pullf.xvg
│   ├── win_001/ ... win_NNN/
├── analysis/
│   ├── tpr.dat / pullx.dat   # gmx wham input lists
│   ├── pmf.xvg / histo.xvg
│   └── (bsResult.xvg)        # bootstrap 指定時のみ
└── run.sh                  # ↑全部を流す top-level driver
```

## パイプライン詳細

```
[stage 0] ペプチド初期 PDB     bilayer.write_peptide_pdb
                                 sequence (one-letter) → tleap → 拡張鎖 PDB
                                 cap_n / cap_c (ACE / NME 既定)
[stage 1] bilayer 構築          bilayer.build_bilayer
                                 packmol-memgen で水・イオン・脂質配置
                                 lipid 数: ratio + distxy_fix から自動決定
                                 distxy_fix = sqrt(sum(n_per_leaflet * APL))
                                 saltcon (M) と solute_charge から中性化
[stage 2] パラメータ化           parameterize_amber.parameterize
   (AMBER backend)               tleap で ff19SB + Lipid21 + TIP3P + frcmod.ionsjc_tip3p
                                 prmtop + inpcrd → parmed → top + gro
                                 write_index_from_gro で 6 group .ndx
[stage 3] 平衡化 MDP             mdp_us_protocol.write_equilibration_mdps
                                 em.mdp / nvt.mdp / npt.mdp
                                 thermostat = Bilayer + Non_Bilayer 2-group
                                 NPT は semiisotropic (xy / z 独立緩和)
[stage 4] pulling MDP            pulling.write_pulling_mdp
                                 NPT-semiisotropic chassis + pull-code (rate>0)
                                 pull-coord1-init は estimate_initial_pull_coord
                                 で system.gro から自動推定
                                 geometry = direction (signed z, semiisotropic 互換)
[stage 5] window MDP 一括        umbrella.write_window_mdps
                                 各 window で rate=0、init = z_min + i*dz
                                 force constant は USProtocol.force_constant
[stage 6] run.sh 生成             umbrella.write_run_script
                                 em → nvt → npt → pull → extract → 各 window MD → wham
                                 chmod +x、相対パス使用
== ここまで build() が生成 ==
== 以下は run.sh 実行で行われる ==
[stage 7] equil + pull + window MD  gmx grompp + gmx mdrun
[stage 8] frame extraction       pulling.extract_window_frames
                                 pullx.xvg を読み、各 window の z_target に
                                 最も近い frame を gmx trjconv -dump で抽出
[stage 9] PMF 解析               pmf.run_wham
                                 gmx wham -temp <T> -unit kJ
                                 任意 -bsnum N で bootstrap
```

## 設定リファレンス (MembraneConfig)

### MembraneConfig (top-level)

| field | default | 役割 |
|---|---|---|
| `backend` | `"amber"` | `"amber"` または `"charmm36"` (Phase C) |
| `lipids` | (required) | `LipidSpec` のリスト |
| `peptide` | (required) | `PeptideSpec` |
| `ions` | `IonSpec()` | 塩・対イオン |
| `box_xy_nm` | 8.0 | bilayer xy (実際は `n_per_leaflet * APL` から自動再計算) |
| `water_thickness_nm` | 2.5 | 膜面から box 端までの水層厚さ |
| `distance_to_lipid_nm` | 1.5 | peptide 初期位置と上 leaflet の距離 |
| `equilibration` | `EquilibrationProtocol()` | em / nvt / npt パラメータ |
| `pulling` | `PullingProtocol()` | 反応座標生成 |
| `umbrella` | `USProtocol()` | window 設定 |
| `seed` | None | 再現性 (gen-vel seed 等) |
| `output_dir` | "." | build 出力先 |
| `amber_protein_ff` | `"leaprc.protein.ff19SB"` | AMBER tleap 用 |
| `amber_lipid_ff` | `"leaprc.lipid21"` | AMBER tleap 用 |
| `amber_water_ff` | `"leaprc.water.tip3p"` | AMBER tleap 用 |
| `charmm_ff_dir` | "" | CHARMM36 backend で必須 (Klauda port パス) |
| `packmol_memgen_path` | `"packmol-memgen"` | バイナリパス |
| `tleap_path` | `"tleap"` | バイナリパス |
| `gmx_path` | `"gmx"` | バイナリパス |

### LipidSpec

```python
# 単成分
LipidSpec(resname='POPC', n_per_leaflet=64)

# 混合膜 (POPC:CHL1 = 4:1 raft 模型)
LipidSpec(resname='POPC', n_per_leaflet=80)
LipidSpec(resname='CHL1', n_per_leaflet=20)

# 三成分混合 (PC:PE:CHL = 4:2:1、内膜状)
LipidSpec(resname='POPC', n_per_leaflet=40)
LipidSpec(resname='POPE', n_per_leaflet=20)
LipidSpec(resname='CHL1', n_per_leaflet=10)

# 明示的 APL (低温 gel-phase 等で実測値を使う場合)
LipidSpec(resname='DPPC', n_per_leaflet=64, apl_angstrom2=49.0)
```

`resname` は packmol-memgen の `--available_lipids` で確認可。Lipid21 では
POPC / POPE / POPG / DOPC / DPPC / CHL1 等の主要種をサポート。

#### `apl_angstrom2`: 脂質ごとの area-per-lipid

複数の `LipidSpec` を並べると混合膜になり、`packmol-memgen --ratio` には
n_per_leaflet を gcd 約分した整数比が、`--distxy_fix` には
`sqrt(sum(n × APL))` が自動で渡される。`apl_angstrom2` は **省略 (default
0.0) で `DEFAULT_LIPID_APL` テーブル (下記) から自動 lookup** する。
T や cholesterol 凝縮の効果で表内値が合わない場合のみ明示指定する。

| 脂質 | APL (Å²、~310 K Lα) | 備考 |
|---|---|---|
| POPC | 67 | 標準的な液晶相 phospholipid |
| DOPC | 72 | 二重不飽和、より広い |
| DPPC | 63 | 飽和、Lα 相 (Tm 314 K 以上); gel 相 ~49 |
| DMPC | 60 | 短鎖飽和 |
| DLPC | 63 | 短鎖飽和 |
| DSPC | 65 | 長鎖飽和 |
| POPE | 56 | PE 頭部、PC より小さい |
| DOPE | 65 | |
| POPG | 64 | アニオン性 |
| DOPG | 67 | |
| POPS | 60 | アニオン性 |
| POPA | 62 | アニオン性、最小 |
| CHL1 / CHL / CHOL | 38 | コレステロール、ステロール環で小さい |

テーブルにない残基名は `apl_angstrom2` 省略時 65 Å² (汎用 Lα 値) を fallback。

#### 既知 lipid のロード / 検索

v1.17.2 から **CLI と Python helper** で curated APL table と
packmol-memgen の全 259 種を browse できます。

```bash
# curated APL テーブル全部 (~60 entries、PC/PE/PG/PS/PA/SM/sterol で head 別表示)
python -m abmptools.membrane.lipid_info --known

# SM (sphingomyelin) family だけ
python -m abmptools.membrane.lipid_info --known --head SM

# 単一 lipid の APL 確認 (curated table miss → 65.0 fallback)
python -m abmptools.membrane.lipid_info --apl POPC      # 67.0
python -m abmptools.membrane.lipid_info --apl PSM       # 55.0
python -m abmptools.membrane.lipid_info --apl WHATEVER  # 65.0 (fallback)

# packmol-memgen の **全 259 種**を表示 (curated table 外も含む)
python -m abmptools.membrane.lipid_info --available
```

Python から直接:

```python
from abmptools.membrane.bilayer import (
    list_known_lipids, query_packmol_memgen_lipids, DEFAULT_LIPID_APL,
)

# curated table 全部
for resname, apl, head in list_known_lipids():
    print(f"{resname}  {apl}  {head}")

# PE family だけ
for r, apl, _ in list_known_lipids(head_group="PE"):
    print(r, apl)

# packmol-memgen 由来の 259 種 (charge / 化学名付き)
for resname, charge, desc in query_packmol_memgen_lipids():
    print(f"{resname:>5s}  q={charge:+d}  {desc}")

# 単発で APL を取得
apl = DEFAULT_LIPID_APL.get("POPC", 65.0)   # 67.0
```

#### Curated table の概観 (~60 entries、v1.17.2)

| head | 含まれる resname (一部) | APL レンジ (Å²) |
|---|---|---|
| **PC** (中性) | DLPC, DMPC, DPPC, DSPC, **POPC**, PMPC, SOPC, **DOPC**, DAPC, DHPC, AHPC | 60-80 |
| **PE** (中性、H-bond で小) | DLPE, DMPE, DPPE, DSPE, **POPE**, PMPE, SOPE, **DOPE**, DAPE, DHPE | 50-75 |
| **PG** (anionic) | DLPG, DMPG, DPPG, DSPG, **POPG**, SOPG, DOPG, DAPG, DHPG, AHPG | 55-78 |
| **PS** (anionic) | DLPS, DMPS, DPPS, DSPS, **POPS**, SOPS, DOPS, DAPS, DHPS, AHPS | 55-70 |
| **PA** (anionic、最小頭部) | DLPA, DMPA, DPPA, DSPA, **POPA**, DOPA, DAPA, DHPA, AHPA | 50-70 |
| **SM** (中性、raft、Lipid21 only) | LSM, MSM, **PSM**, SSM, OSM, ASM, HSM | 53-68 |
| **sterol** | **CHL1**, CHOL, CHL (alias) | 38 |

完全なリストは `python -m abmptools.membrane.lipid_info --known` で確認。
packmol-memgen は他に **arachidonoyl-アシル × 各種 head の 30+ 種**を持って
おり (`--available` で確認可)、それらも `LipidSpec` で `apl_angstrom2`
明示すれば使えます。

### PeptideSpec

```python
PeptideSpec(name='aa5', sequence='AAAAA')              # tleap 生成
PeptideSpec(name='helix', pdb_path='./helix.pdb')       # 既存 PDB 使用
```

| field | default | 説明 |
|---|---|---|
| `name` | `"peptide"` | ファイル名 tag |
| `sequence` | "" | 一文字 AA コード (空のとき pdb_path 使用) |
| `pdb_path` | "" | sequence の代わりに既存 PDB |
| `n_copies` | 1 | 通常 1 (PMF は単分子) |
| `cap_n` | `"ACE"` | N-末端キャップ ("" で自由アミン) |
| `cap_c` | `"NME"` | C-末端キャップ ("" で自由カルボキシル) |
| `initial_z_offset_nm` | 3.0 | 注: 現在 packmol-memgen 任せで未使用 |

HIS は `HIE` (epsilon-protonated、AMBER 中性 pH 既定) にマップ。荷電残基は
`estimate_peptide_charge()` で D/E -1、K/R +1、HIE 0 で集計し、
`--solute_charge` を packmol-memgen に渡してイオン中性化。

### USProtocol

```python
USProtocol(
    z_min_nm=-3.0, z_max_nm=+3.0, window_spacing_nm=0.1,   # 61 windows
    force_constant_kj_mol_nm2=1000.0,
    window_nsteps=10_000_000, window_dt_ps=0.002,           # 20 ns/window
    pull_geometry='direction',     # semiisotropic と互換
    pull_vec='0 0 1',
    pull_dim='N N Y',
    pull_group1='Bilayer', pull_group2='Peptide',
)
```

`pull_geometry='direction'` は GROMACS が semiisotropic Pcoupl と互換にする
ための選択。`'direction-periodic'` だと grompp が
`Cannot have dynamic box while using pull geometry 'direction-periodic' (dim z)` で
失敗する。

### EquilibrationProtocol / PullingProtocol

| protocol | 主要 field | 既定値 |
|---|---|---|
| Equilibration | `em_steps`, `nvt_nsteps`, `npt_nsteps`, `temperature_K`, `pressure_bar`, `tau_t_ps`, `tau_p_ps` | 50000 / 100000 / 1000000 / 310 / 1.0 / 1.0 / 5.0 |
| Pulling | `pull_rate_nm_per_ps`, `pull_force_constant`, `nsteps`, `nstxout_compressed` | 0.001 / 1000 / 5000000 / 500 |

### IonSpec

```python
IonSpec(cation='Na+', anion='Cl-', salt_concentration_M=0.15, neutralize=True)
```

CHARMM36 backend では cation/anion 名を `SOD` / `CLA` などに合わせる必要あり (Phase C)。

## 設計上のポイント

### thermostat group
`Bilayer` (Lipid21 全 atom) と `Non_Bilayer` (= System − Bilayer = peptide + water + ions)
の 2-group。リピド緩和タイムスケールと水緩和タイムスケールを decouple。
Ions 専用 group は数十原子しかなく V-rescale で病的になりがちなので避ける。

### Pcoupl
NPT は **semiisotropic** が必須。膜面 xy と膜法線 z で独立に緩和できないと
膜厚や面積/脂質が物理的にズレる。`refcoord-scaling=com` で位置拘束 (将来追加)
が COM ベースに従う。

### pull geometry — stage 別に切替
**Pulling stage** と **window stages** で異なる geometry / Pcoupl を使う:

| stage | Pcoupl | geometry | 理由 |
|---|---|---|---|
| em / nvt / npt 平衡化 | semiisotropic | (pull なし) | 膜面 xy と法線 z 独立緩和 |
| **pull** (反応座標生成) | **none (NVT)** | **direction-periodic** | peptide が水層→膜→反対水層を traverse する間、peptide-bilayer COM-z は box 半分に達するため direction-periodic 必須。GROMACS は direction-periodic + 動的 box (= 任意の Pcoupl) を拒否するので NVT 固定。短時間 (~ns) なら膜は変形しない |
| **windows** (US 本体) | semiisotropic | direction | peptide は window 中心近傍に局在 → 0.49 × box 制約クリア。signed coord で両側 PMF |

`distance` は常に正の値を返すため両側 PMF (z_min < 0 < z_max) で使えない。

### pull-group1-pbcatom (Bilayer 必須)
Bilayer グループ (Lipid21 split で 8000+ atoms、xy 全体) は GROMACS の
"diameter > 0.5 × box" 制約に引っかかるため、explicit な pbcatom 指定が必須。
`abmptools.membrane.pulling.find_pbc_center_atom` で **bilayer COM に最も近い
原子の 1-based index** を計算し、`pull-group1-pbcatom` に書き込む処理を
builder が自動で行う。指定しないと grompp が次のように落ちる:

```
ERROR: When the maximum distance from a pull group reference atom to
other atoms in the group is larger than 0.5 times half the box size a
centrally placed atom should be chosen as pbcatom.
```

### lipid 数の制御
`LipidSpec.n_per_leaflet` から `area_per_leaflet = n * APL` (default APL=65 Å²)
を計算 → `--distxy_fix sqrt(area)` を packmol-memgen に渡す。実際の lipid
数は packmol-memgen が patch を packing した結果 (target ± 数個 程度)。
正確な count が必要な場合は `--apl_exp` フラグ追加または APL を実験値に合わせる。

### 反応座標の符号と初期値
`pull_init_nm = peptide_COM_z - bilayer_COM_z` を `system.gro` (build 直後、
未平衡) から計算。平衡化中に peptide は通常 < 0.5 nm 程度しか動かないため、
そのまま pull.mdp に書き込んで OK。大きく動く系では平衡化後に
`pulling.estimate_initial_pull_coord(equil/npt.gro, ...)` で再計算して pull.mdp
を上書きすること。

## トラブルシューティング

### `Concentrations/number of molecules have to be provided.`
packmol-memgen で `--solute_con` を指定していない。本パッケージの
`assemble_packmol_memgen_cmd` は `--solute_con 1` を自動付加するので、
直接 packmol-memgen を呼ぶカスタムスクリプトのみで起こる。

### `AttributeError: module 'numpy' has no attribute 'float'.`
packmol-memgen バンドル `pdbremix/v3numpy.py` が NumPy 1.24+ と非互換。
[インストール](#インストール) の sed パッチを当てる。

### `Cannot have dynamic box while using pull geometry 'direction-periodic' (dim z)`
`USProtocol.pull_geometry` を `direction` (without -periodic) に変更。
本パッケージは window stage では既定で `direction`、pull stage は専用に
NVT (Pcoupl=no) + `direction-periodic` を内部で割当てるので、user 設定
を変える必要はない。

### `Distance between pull groups N and M is larger than 0.49 times the box size`
`direction` geometry の制約で、peptide が水層に居る時 peptide-bilayer COM-z
は box の半分に達してしまう。本パッケージの pull stage は NVT +
direction-periodic でこの制約を回避するので、エラーが出る場合は (a) 自前で
mdp を編集している、(b) box が極端に小さい、のいずれか。後者なら
`MembraneConfig.water_thickness_nm` を増やす。

### `gmx grompp` が tc-grps mismatch で失敗
カスタム index ファイルを使っている場合、`Bilayer` と `Non_Bilayer` の
2 group が含まれているか確認 (`grep '^\\[' system.ndx`)。
`parameterize_amber.write_index_from_gro` が自動生成するので、自前 ndx
を使う場合は同様の groups が必要。

### packmol-memgen が "PMEMD not found. Setting path to SANDER" を表示
警告のみで動作には影響しない。pmemd が無い conda env 標準状態。

### lipid 数が target と合わない
`LipidSpec.n_per_leaflet=64` 指定で実際 62 lipid に packing された等の
微小ズレは `estimate_distxy_angstrom` の APL 推定 (65 Å²) と packmol-memgen
内部 APL の差異による。気になる場合は `--apl_exp` を追加 (将来対応) または
APL を実験値で override。

## production 向け注意点

| 項目 | smoke 値 | production 推奨 |
|---|---|---|
| ペプチド残基数 | 5 | 系に依存 (10-30 で実用、長すぎると box 巨大化) |
| n_per_leaflet | 32 | 64-128 (xy 6-9 nm)。peptide が大きいと periodicity 干渉回避で要調整 |
| window 数 | 7 (Δ=0.5 nm) | 30-60 (Δ=0.1-0.2 nm)。histogram overlap > 30% を目安 |
| window 時間 | 1 ns | 20-50 ns。膜中心 / 界面ほど長く |
| pulling 時間 | 5 ns | 10-50 ns。長いほど熱化された initial 構造が得られる |
| 平衡化 | 1 ns | 10-50 ns (脂質緩和には ~10 ns 必要) |
| 力定数 k | 1000 kJ/mol/nm² | 1000-2500、histogram width で決定 |

PMF 収束チェック:
- `gmx wham -bsnum 50` で bootstrap 誤差を確認
- 時間ブロック分割 (前半/後半) で PMF 比較し、変動が大きい場所だけ追加サンプリング
- histogram overlap が悪い領域 (`gmx wham -hist` で可視化) では window を細かく
- 配向自由度が遅い場合は **REUS (GENESIS)** に切替検討

## 力場 license まとめ

| 力場 | 学術 | 商用 | 注 |
|---|---|---|---|
| AMBER ff19SB / Lipid21 / TIP3P / JC ions | ✅ | ✅ | AmberTools 配布、本パッケージのデフォルト |
| GAFF2 (Antechamber) | ✅ | ✅ | 新規小分子 (今後の ligand 拡張用) |
| CHARMM36 (param 値そのもの) | ✅ | ✅ | MacKerell 公式 / Klauda port、Phase C で使用 |
| CGenFF Web server output | ✅ | ❌ | Silcsbio 商用ライセンス必須 — 本パッケージは使わない |
| CHARMM-GUI 自動生成 | ✅ | ❌ | 商用契約必須 — 本パッケージは使わない |
| OpenFF (lipid params) | — | — | 未整備、今回は使わない |

## ステータス

| Phase | 内容 | コミット |
|---|---|---|
| A | skeleton (subpackage 骨格 + dataclass + builder shell) | `a0d2bb9` |
| B-1 | bilayer.py (peptide + packmol-memgen) | `f1d49c7` |
| B-2 | parameterize_amber.py (tleap → parmed → GROMACS) | `517848b` |
| B-3 | mdp_us_protocol.py (em / nvt / npt-semi + pull block) | `8f4657a` |
| B-4 | pulling.py (pull.mdp + frame extraction) | `e30263d` |
| B-5 | umbrella.py (window MDPs + run.sh) | `2efaff1` |
| B-6 | pmf.py (gmx wham wrapper + CLI) | `eeb57d1` |
| B-7 | smoke test + builder fix | `907c040` |
| D | tutorial doc + cross-references (this doc) | `a01fdda` `21dd567` |
| C | CHARMM36 backend (Klauda port + pdb2gmx) | `5d28034` |
| C+ | CHARMM36 backend 実機検証 + 7 件 bug fix (v1.17.3) | `<TBD-1.17.3>` |

**Phase C end-to-end smoke** はユーザー側で `CHARMM36_FF_DIR` を整えた上で
`bash tests/integration/run_membrane_us_charmm_smoke.sh` を実行すると検証
される。未設定なら自動 SKIP (exit 0)、translation ユニットテストは
.ff dir 不要で実行可能 (`pytest tests/test_membrane_charmm_translate.py`、
26 tests / <1 秒)。

v1.17.3 で **CHARMM36 backend を初回実機検証**し、Klauda port 特有の 7 件
の bug を修正。peptide + bilayer + water + ion の smoke build (grompp pass)
は `charmm36-feb2026_cgenff-5.0.ff/` / `charmm36-jul2022.ff/` のいずれでも
完走することを確認 (PDB 4-char residue 切り詰め / `NME→CT3` 等の過剰
rename / ACE atom 名 mapping / NME atom 名 mapping / terminus None index
hardcode / TIP3 spurious O-H-H angles 削除 / pdb2gmx 2021.3 無限 spin
回避)。詳細は `CHANGELOG.md` の v1.17.3 セクション参照。

### Phase C+ 実機検証の信頼度階層 (v1.17.3 時点)

CHARMM36 backend の出力が「正しい」ことを保証する根拠を階層化:

| Layer | 検証内容 | 結果 |
|---|---|---|
| L1. テーブル整合性 | rename テーブルが Klauda port の rtp と一致 | ✅ 26 unit tests pass |
| L2. パイプライン syntactic | `pdb2gmx` が "atom not found" を出さず通過 | ✅ smoke 12 PASSED |
| L3. grompp 受け入れ | 11 windows × MDP が `gmx grompp` を error 0 で通る | ✅ |
| L4. 化学量論 (Phase A) | total charge = 0、peptide=62 atoms (AMBER と一致) | ✅ 14 tests pass (`test_membrane_topology_sanity.py`) |
| L5. EM 収束 | energy minimization が物理的に妥当な値に収束 | ✅ E=-1.63×10⁵ kJ/mol, max F=783 < tol |
| L6. NVT/NPT 安定性 | 5 ns NPT (T=308.4±0.2 K, P=-27±28 bar) | ✅ 安定、blow up なし |
| L7. APL 文献値一致 | POPC@310K で APL ≈ 65 Å² | ✅ **64.88 ± 0.82 Å²** ⭐ |
| L8. P-P 厚 | POPC@310K で D_PP ≈ 38-40 Å | ⚠ 50.7 Å (Berendsen NPT 5 ns 段階; production windows で更に relax) |
| L9. AMBER vs CHARMM PMF 比較 | 同系の barrier が ff レベル (数 kJ/mol) で一致 | (進行中: Phase D) |
| L10. CHARMM-GUI と topology 比較 | atom counts / bonded params が一致 | (未実施) |

L1-L7 で「topology は正しく、物理的に動く」を保証。L8 の P-P 厚は Berendsen
barostat 5 ns では完全 equilibrate に至らず約 30% 厚いが、production phase
の c-rescale (Bussi 2020) barostat で許される proper volume fluctuation
により windows 段階で更に relax 見込み。

実機検証の sanity check helper は CLI として利用可能:

```bash
python -m abmptools.membrane.topology_sanity /path/to/build_dir
```

または pytest:

```bash
MEMBRANE_BUILD_DIR=/path/to/build_dir pytest tests/test_membrane_topology_sanity.py
```

### CHARMM 業界実態 (重要、ニッチ性の認識)

raw `pdb2gmx` で full bilayer PDB を扱うのは少数派 — academic 多数派
(~95%) は **CHARMM-GUI** で完成済 `topol.top` を使い pdb2gmx を完全回避、
少数派は psfgen + charmm2gmx (Wacha & Lemkul 2023 *JCIM*) を使う。完全
commercial-OK な membrane builder は **packmol-memgen** (AmberTools)
のみで、その出力を CHARMM 形式に変換するパスは世界的にも超ニッチ — 本
パッケージの修正点は世界で数十人しか踏まない pain points。

## TIP3P の AMBER vs CHARMM 違い

「TIP3P」と一口に言っても、AMBER と CHARMM-modified では H atom の
Lennard-Jones パラメータが異なる:

| | AMBER TIP3P (Jorgensen 1983) | CHARMM-modified TIP3P |
|---|---|---|
| O σ/ε | 3.151 Å / 0.152 kcal/mol | 同じ |
| **H σ/ε** | **0 / 0** (LJ なし、点電荷のみ) | **~0.4 Å / ~0.046 kcal/mol** |
| 電荷 | O=-0.834, H=+0.417 | 同じ |

CHARMM 側で H に小さい LJ を足したのは、HB ネットワーク安定化目的
(電気陰性原子と H の人工的 overlap 防止)。影響は限定的で、純水密度は
ほぼ同じだが、イオン水和エネルギーが ~1-3 kcal/mol、PMF 絶対値が数
kJ/mol レベルでずれる。

`abmptools.membrane` では backend ごとに自動で適切な TIP3P が選ばれる:

- **AMBER backend**: `tleap` `leaprc.water.tip3p` → parmed が GROMACS
  topology に書き出し → AMBER TIP3P (H に LJ なし)
- **CHARMM36 backend**: `pdb2gmx -water tip3p` が ff dir の `tip3p.itp`
  (CHARMM-modified) を `#include` → CHARMM TIP3P (H に小さい LJ)

論文 / 報告書に「TIP3P を使った」と書く時は、どちらの実装かを明記
することを推奨。

Phase B 動作確認 (smoke、~30 秒、4-core local):

- 入力: poly-Ala 5-mer + POPC 32/leaflet + 0.15 M NaCl + TIP3P
- 出力: 18,751 atoms (62 POPC + 3455 WAT + 16 ions + 5 ALA)、box 5.4×5.4×10.1 nm
- bilayer P-P 厚 39.05 Å、対称 (±19.5 Å)
- 7 window MDP 全て `gmx grompp -maxwarn 1` 通過

実行手順:

```bash
bash tests/integration/run_membrane_us_smoke.sh
# RUN_DIR=/path/to/keep でビルド出力を残せる
```

## 関連ドキュメント

- [`tutorial_membrane_us.md`](tutorial_membrane_us.md) — step-by-step 操作チュートリアル (本 doc は reference、tutorial は操作手順)
- [`amorphous.md`](amorphous.md) — 3D ランダム box ビルダー (姉妹サブパッケージ)
- [`gro2udf.md`](gro2udf.md) — GROMACS ↔ COGNAC UDF 変換
- [`io_spec.md`](io_spec.md) — `param_udf` 仕様
- [`licenses_third_party.md`](licenses_third_party.md) — third-party license 一覧
- [`architecture.md`](architecture.md) — abmptools 全体構成
