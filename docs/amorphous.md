# abmptools.amorphous — Amorphous Structure Builder

API + 短鎖ポリマー（オリゴマー）混合系のアモルファス構造ビルダー。
OpenFF でパラメータ化 → Packmol で初期配置 → GROMACS アニーリング MD 用ファイル一式を生成。

> **初めて使う方へ**: [amorphous_tutorial.md](amorphous_tutorial.md) の
> hands-on チュートリアル (同梱サンプル→自前系統→PubChem→JSON の 4 パターン
> + 出力の歩き方 + つまずきポイント対処) を先に読むのがおすすめ。
> 本ドキュメントは CLI/API のリファレンスです。
>
> **Sister subpackage**: [`abmptools.membrane`](membrane.md) は 2D 周期の脂質
> 二重膜系 (peptide-bilayer Umbrella Sampling / PMF) 向け。本モジュールは
> 3D ランダムボックス系 (アモルファス・混合溶媒) 向け。両者とも
> `abmptools.core.SystemModel` を共有するが、stage 構成と Pcoupl 規約が異なる。

## インストール

```bash
# 基本依存 (conda推奨)
conda install -c conda-forge openff-toolkit openff-interchange openmm numpy

# 電荷計算 (どちらか一方)
conda install -c conda-forge ambertools    # AM1-BCC (推奨)
# or
pip install openff-nagl                     # ML ベース (AmberTools 不要)

# Packmol
conda install -c conda-forge packmol
# or ソースから: http://m3g.iqm.unicamp.br/packmol/

# pip extras
pip install abmptools[amorphous]
```

## クイックスタート

### CLI

```bash
# SMILES入力
python build_amorphous.py \
    --smiles "CCCCC" "c1ccccc1" \
    --name pentane benzene \
    --n_mol 200 50 \
    --density 0.8 \
    --temperature 300 \
    --output_dir ./pentane_benzene

# SDF入力
python build_amorphous.py \
    --mol api.sdf polymer.sdf \
    --n_mol 200 40 \
    --box 6.0 \
    --temperature 300

# PubChem 3D SDF を直接取得 (abmptools >= 1.15.3)
python build_amorphous.py \
    --pubchem_cid 3825 \
    --name ketoprofen \
    --n_mol 50 --density 0.8 \
    --output_dir ./ketoprofen_pubchem
# 同様に名前指定も可能 (--pubchem_name "aspirin")
# ダウンロードした SDF は <output_dir>/input/pubchem_cid3825.sdf にキャッシュ

# 重量分率モード
python build_amorphous.py \
    --mol api.sdf polymer.sdf \
    --weight_fraction 0.7 0.3 \
    --total_molecules 300

# JSON設定
python build_amorphous.py --config mixture.json

# python -m でも実行可能
python -m abmptools.amorphous --smiles C --name methane --n_mol 10 --box 2.0
```

### Python API

```python
from abmptools.amorphous import AmorphousBuilder, BuildConfig, ComponentSpec

config = BuildConfig(
    components=[
        ComponentSpec(smiles="CCCCC", name="pentane", n_mol=200),
        ComponentSpec(smiles="c1ccccc1", name="benzene", n_mol=50),
    ],
    density_g_cm3=0.8,
    temperature=300,
    output_dir="./output",
)
builder = AmorphousBuilder(config)
result = builder.build()

print(result["gro"])        # system.gro のパス
print(result["top"])        # system.top のパス
print(result["box_nm"])     # ボックスサイズ [nm]
```

### Force field override (water 等の専用 FF)

`BuildConfig.forcefield` は `str` (単一 FF) または `list[str] | tuple[str]`
(stacked FF) を受け付ける。後者は OpenFF ForceField の SMIRKS-overlay 機構で、
**後ろの FF が前を上書き** する。water 含む系で OpenFF/GAFF organic FF の
σ/ε が repulsive 過剰になり、純 water box が 1.0 → 0.26 g/cm³ に膨張する
既知問題への根本対策に使う:

```python
from abmptools.amorphous import AmorphousBuilder, BuildConfig, ComponentSpec

config = BuildConfig(
    components=[
        ComponentSpec(smiles="CO", name="A_methanol", n_mol=100),
        ComponentSpec(smiles="O",  name="B_water",    n_mol=100),
    ],
    density_g_cm3=0.6,
    temperature=300,
    forcefield=[
        "openff_unconstrained-2.1.0.offxml",  # organic 全般
        "tip3p.offxml",                       # water 分子のみ上書き
    ],
    output_dir="./output_meoh_water",
)
result = AmorphousBuilder(config).build()
```

利用可能な OpenFF 公式 water FF (`openff-forcefields` パッケージ同梱):
- `tip3p.offxml` — 古典 TIP3P
- `tip3p_fb.offxml` — TIP3P-FB (refit)
- `spce.offxml` — SPC/E

実機検証 (2026-04-29): pure water (200 TIP3P, 300 K, 5-stage MD) の anneal
stage で **density = 0.991 g/cm³**、TIP3P literature (~0.99 g/cm³ at 298 K)
と完全一致。同条件の OpenFF organic 単独だと 0.26 g/cm³ に縮む。

### Rigid cluster center (`cluster_pdb_path` + `frozen_atom_indices`)

cluster center (e.g. water trimer の H-bond triangle) を MD 中構造維持したい
場合の 2 段構え:

1. **`cluster_pdb_path`** — pre-built rigid cluster の PDB ファイル。 packmol が
   `fixed <center> 0 0 0` constraint で box 中央に rigid block 配置する。
   first component の `n_mol` から自動で cluster 分減算 (cluster がその mol を
   提供するため、 重複防止)。
2. **`frozen_atom_indices`** — 1-based global atom indices (GROMACS ndx
   convention)。 build 時に:
   - `system.ndx` に `[ FrozenAtoms ]` index group を追加
   - first moleculetype 内に `[ position_restraints ]` を `#ifdef
     POSRES_TRIMER` ガード付きで追加 (`posres_force_constant` で `kx`/`ky`/`kz`)
   - 02_nvt_highT 以降全 mdp に `define = -DPOSRES_TRIMER` を追加 (EM は skip、
     packmol overlap relax 段階での numerical 不安定を回避)
   - homo pair で first moletype 件数が cluster mol 数より多い場合は moletype
     を split: `<name>_TRIMER` (cluster 分、 posres 入り) + `<name>` (残り、 posres
     無し)。 split 時に GROMACS の "Only one moltype with [settles] allowed"
     制約を回避するため、 TRIMER 側の `[ settles ]` を等価な `[ constraints ]`
     (LINCS-based、 同じ O-H/H-H 距離) に変換。

```python
from abmptools.amorphous import AmorphousBuilder, BuildConfig, ComponentSpec

config = BuildConfig(
    components=[
        ComponentSpec(smiles="O",  name="B_water",    n_mol=3),
        ComponentSpec(smiles="CO", name="A_methanol", n_mol=300),
    ],
    forcefield=["openff_unconstrained-2.1.0.offxml", "tip3p.offxml"],

    # 水 trimer (H-bond triangle、 O-O = 0.28 nm) を box 中央 rigid 配置
    cluster_pdb_path="./cluster_xyz/water_trimer.pdb",

    # 3 waters の O 原子 (atom 1, 4, 7) を posres で固定
    frozen_atom_indices=[1, 4, 7],
    posres_force_constant=10000.0,  # kJ/mol/nm² (default 10000)

    density_g_cm3=0.6, temperature=300, T_high=400,
    seed=42, output_dir="./output_trimer_meoh",
)
result = AmorphousBuilder(config).build()
```

**Rigid cluster PDB の format 要件** — Interchange の substructure match と
互換にするため、 packmol input PDB は:
- residue name は base component と同じ (e.g. `UNL`、 `HOH` は不可)
- chain ID 必須 (col 22 に `A` 等)
- atom name は PDB strict alignment (1-char element は col 13 空白 + col 14
  に element + col 15-16 に index、 例: `' O1 '`)
- `CONECT` records 必須 (O-H1, O-H2 bonds、 Interchange の bond perception 用)

obabel の `xyz -> pdb` 直変換だと resname=`HOH` / 各 H が resid=0 で fail する
ので手書き推奨。 fcews-manybody repo の
`cluster_xyz/water_trimer.pdb` が water trimer の互換版 (UNL+A+CONECT)。

**実機検証** (methanol-water B-A pair, 500 ps NVT, posres k=10000):
- 初期 (packmol fixed): O-O = 0.275, 0.275, 0.276 nm (H-bond triangle)
- 500 ps MD 後: O-O = 0.275, 0.292, 0.311 nm (構造維持、 wiggle < 0.04 nm)

UDF route の `clusterfix` (Steady velocity=0 atom 拘束) と機能等価。 freezegrps
(hard fix) は LINCS/SETTLE + Domain Decomposition で `determinant = -inf` で
abort するので採用せず、 harmonic posres を採用。

### Multi-component pair の thermostat / annealing 動作

複数 `ComponentSpec` の組み合わせ (mixed pair) と同名 component の重複
(pure-component pair) の両方で、`mdp` / `ndx` / `run_all.sh` が GROMACS の
per-group 仕様に整合するよう自動補正される:

#### Pure-component pair (同名 component 2 個、例: B_water-B_water)

`ComponentSpec(name="X")` を 2 個並べた場合 (A-B χ 計算で A=B にする
self-pair) は内部で dedup:

- `mdp` の `tc-grps` は `"X"` (1 group) → `ref-t` / `tau-t` も単一値で整合
- `system.ndx` の `[ X ]` group は両 component の atom を集約 (前 fix まで
  片方が overwrite されて atom 半分しか group 化されない fail mode あり)

#### Mixed-component pair (例: A_methanol-B_water)

異なる name の 2 component では dedup されず `tc-grps = "A_methanol B_water"`
(2 group) になるため、GROMACS の per-group 仕様に従って以下が自動で n 倍に
複製される:

- `tau-t` / `ref-t` (`_thermostat_block` 内): `"0.1"` → `"0.1 0.1"`、
  `"600.0"` → `"600.0 600.0"`
- `annealing` / `annealing-npoints` / `annealing-time` / `annealing-temp`
  (stage 4 の `generate_anneal_mdp`): 例として 2 group では
  `annealing = single single`、`annealing-temp = 600.0 300.0 600.0 300.0`
  のように複製。各 group は同じ 600 K → 300 K cooling schedule

これらは `tc-grps.split()` の token 数 (`n_groups`) で機械的に決まる。
1 group (legacy "System" or 単成分) では token 単一で従来出力と byte-identical。

#### MD ランチャー (`run_all.sh`)

default で `MDRUN_OPTS="${MDRUN_OPTS:--ntmpi 1}"` を使う (DD 無効化、OpenMP
は維持)。凝縮系で box が cooling 中に縮んだ際の
`box size ... is too small for a cut-off ... with 3 domain decomposition
cells` を回避するため。override したい場合は
`MDRUN_OPTS="-nt 8" bash run_all.sh` のように環境変数で渡す。

## CLI 引数

| 引数 | 説明 | デフォルト |
|------|------|-----------|
| `--mol` | SDF/MOL ファイル（複数指定） | - |
| `--smiles` | SMILES 文字列（複数指定） | - |
| `--name` | 分子名（省略時は自動） | - |
| `--n_mol` | 各成分の分子数 | - |
| `--weight_fraction` | 各成分の重量分率 | - |
| `--total_molecules` | 重量分率モード時の総分子数 | 200 |
| `--box` | ボックスサイズ [nm]（0=自動） | 0 (自動) |
| `--density` | 目標密度 [g/cm^3] | 0.8 |
| `--temperature` | 最終温度 [K] | 300 |
| `--T_high` | 高温 [K] | 600 |
| `--seed` | 乱数シード | None |
| `--forcefield` | OpenFF力場 | openff_unconstrained-2.1.0.offxml |
| `--packmol_tolerance` | Packmol 重なり距離 [A] | 2.0 |
| `--packmol_path` | Packmol バイナリパス | packmol |
| `--output_dir` | 出力ディレクトリ | . |
| `--config` | JSON設定ファイル | - |
| `-v` | 詳細ログ | off |

## 出力ファイル

```
{output_dir}/
├── input/
│   ├── config.json         # 再現用設定
│   ├── component_0.pdb     # 単分子 PDB
│   └── component_1.pdb
├── build/
│   ├── packmol.inp         # Packmol入力
│   ├── mixture.pdb         # Packmol出力
│   ├── system.gro          # GROMACS座標
│   ├── system.top          # GROMACS トポロジー
│   └── system.ndx          # インデックスグループ
└── md/
    ├── 01_em.mdp           # エネルギー最小化
    ├── 02_nvt_highT.mdp    # 高温NVT
    ├── 03_npt_highT.mdp    # 高温NPT
    ├── 04_anneal.mdp       # シミュレーテッドアニーリング
    ├── 05_npt_final.mdp    # 最終NPT平衡化
    ├── run_all.sh          # GROMACS実行スクリプト (5-stage 順次)
    ├── wrap_pbc.sh         # PBC ラップ後処理 (VMD 向け, 1.15.2+)
    └── gen_for_udf.sh    # OCTA viewer 用 energy.xvg + nojump gro 抽出 (1.30+)
```

ビルド後に MD を走らせ、`wrap_pbc.sh` を実行すると以下も生成されます:

```
md/
├── 02_nvt_highT_pbc.xtc    # 分子単位でラップした xtc (-pbc mol -ur compact)
├── 03_npt_highT_pbc.xtc
├── 04_anneal_pbc.xtc
├── 05_npt_final_pbc.xtc
└── 05_npt_final_pbc.gro    # 最終構造 (VMD 初期フレーム用)
```

`gen_for_udf.sh` を実行すると、OCTA viewer で読み込み可能な以下も生成されます:

```
md/
├── 05_npt_final_energy.xvg  # gmx energy: 全エネルギー term (seq 50 で 1-50 番選択)
└── 05_npt_final_nojump.gro  # gmx trjconv -pbc nojump: 分子を分断せず連続軌跡
```

`gmx trjconv -pbc nojump` は、`-pbc mol` (`wrap_pbc.sh`) と違って **分子を box
内に wrap せず、PBC を跨いで連続的に追跡**します。OCTA viewer (GOURMET) で軌跡を
再生する用途に適しています。

## MDP プロトコル（5段階）

| Stage | ファイル | アンサンブル | T | P | Steps | 目的 |
|-------|---------|-----------|---|---|-------|------|
| 1 | 01_em.mdp | - | - | - | 50000 | 急降下法 EM |
| 2 | 02_nvt_highT.mdp | NVT | 600K | - | 100000 | 高温 NVT (V-rescale) |
| 3 | 03_npt_highT.mdp | NPT | 600K | 1bar | 200000 | 高温 NPT (P-R) |
| 4 | 04_anneal.mdp | NPT | 600→300K | 1bar | 500000 | SA: annealing=single |
| 5 | 05_npt_final.mdp | NPT | 300K | 1bar | 500000 | 最終平衡化 |

## ビルド後のワークフロー

```bash
# 1. ビルド (SMILES または SDF 入力)
python build_amorphous.py --smiles "..." --name ... --n_mol ... --output_dir ./run1
# 2. MD 実行 (GROMACS 必須、CPU 8 コアで 1.3 ns 系なら ~10 分)
cd run1/md && bash run_all.sh
# 3. PBC ラップ (VMD 向けトラジェクトリ生成)
bash wrap_pbc.sh
# 4. VMD で可視化
vmd 05_npt_final_pbc.gro -xtc 05_npt_final_pbc.xtc
# アニーリング過程を見たい場合:
vmd 04_anneal.tpr -xtc 04_anneal_pbc.xtc
# 5. OCTA viewer 用 energy.xvg + nojump gro を抽出
bash gen_for_udf.sh
# → 05_npt_final_energy.xvg + 05_npt_final_nojump.gro を生成
# 6. abmptools.gro2udf で OCTA UDF/BDF に変換 (Step 7 参照)
```

`wrap_pbc.sh` は `gmx trjconv -pbc mol -ur compact` を各 xtc と最終 gro に適用し、
箱境界で分断された分子を修復したうえで compact box 表示に揃えます。

`gen_for_udf.sh` は OCTA viewer で読み込み可能な 2 種類の出力を生成します:

- `<stage>_energy.xvg` : `seq 50 | gmx energy -f <stage>.edr -o <stage>_energy.xvg`
  で全 energy term を一括取得 (1〜50 番、存在しない番号は gmx が無視)
- `<stage>_nojump.gro` : `echo 0 | gmx trjconv -f <stage>.{trr,xtc} -s <stage>.tpr
  -pbc nojump -o <stage>_nojump.gro` で分子を box 内に wrap せず、PBC を
  跨いで連続的に追跡した multi-frame gro

## OCTA UDF/BDF への変換 (`abmptools.gro2udf`)

GROMACS の出力 (top + gro + mdp + xvg) を OCTA cognac UDF に変換すると、
**OCTA viewer (GOURMET)** で topology + trajectory + energy plot を 1 ファイル
から再生できます。詳細リファレンスは [gro2udf.md](gro2udf.md) を参照。

### 標準フロー (multi-frame + xvg 全部入り)

```bash
cd run1   # build / md / input/ がある directory

# 全部入り UDF: topology + 101 frame trajectory + xvg energy + mdp Time params
python -m abmptools.gro2udf --from-top build/system.top md/05_npt_final.gro \
    --mdp md/05_npt_final.mdp \
    --trajectory md/05_npt_final_nojump.gro \
    --energy md/05_npt_final_energy.xvg \
    --out 05_full.udf
```

生成された `05_full.udf` に書込まれる内容:

| section | 内容 |
|---|---|
| `Set_of_Molecules.molecule[]` | 全分子の atom (element + 力場 type 名) / bond / angle / torsion |
| `Molecular_Attributes` | atom type / bond / angle / torsion potential、interaction site type |
| `Interactions.Pair_Interaction[]` | LJ σ/ε (per atom-type pair) |
| `Structure (record 0..N-1)` | 各 frame の atom 座標 + Cell + Time + Steps |
| `Statistics_Data` (record ごと) | Energy.{Bond, Angle, Torsion, Nonbonding, Electrostatic, Potential, Kinetic, Total, Hamiltonian} / Temperature / Pressure / Density / Volume、それぞれ Instantaneous / Batch_Average / Total_Average の 3 field |
| `Simulation_Conditions.Dynamics_Conditions.Time` | mdp 由来の `delta_T` (= dt) / `Total_Steps` (= nsteps) / `Output_Interval_Steps` (= nstxout-compressed) |

OCTA viewer (GOURMET) で `05_full.udf` を開くと:
- 3D viewer で trajectory アニメーション再生
- atom テーブルで `Atom_Name` (`C`/`H`/`O`) + `Atom_Type_Name` (`MOL0_X` / `c3` 等)
  + per-mol local `Atom_ID` 表示、元素ベース CPK 色描画
- 統計データプロットで全フレームの Total エネルギー / 温度 / 圧力 / 密度などを
  時系列表示

### topology のみ (skeleton UDF、reference 用)

trajectory / energy を別途 OCTA viewer で attach する流れにしたい場合は
`--topology-only` 単独で skeleton UDF を生成 (Structure record 0 件):

```bash
python -m abmptools.gro2udf --from-top build/system.top md/05_npt_final.gro \
    --mdp md/05_npt_final.mdp --topology-only --out 05_topology.udf
```

初期 1 frame も含めたい場合は `--initial-gro` を併用:

```bash
python -m abmptools.gro2udf --from-top build/system.top build/system.gro \
    --mdp md/05_npt_final.mdp --topology-only \
    --initial-gro md/05_npt_final.gro \
    --out 05_topology_with_initial.udf
```

### OCTA8.4 / J-OCTA-9.1-Student 環境向け

OCTA8.4 / J-OCTA-9.1-Student は cognac10.1 までしか同梱していないため、
default 出力 (cognac11.2 schema 要求) は読めません。`--cognac-version 101` を
付けて bundled の cognac10.1 互換 template を auto-select させます:

```bash
python -m abmptools.gro2udf --from-top build/system.top md/05_npt_final.gro \
    --mdp md/05_npt_final.mdp --cognac-version 101 \
    --trajectory md/05_npt_final_nojump.gro \
    --energy md/05_npt_final_energy.xvg \
    --out 05_full.udf
```

cognac10.1 / cognac11.2 とも gro2udf 出力データは **完全に同一** (Cell.a /
Statistics_Data 各 field 等)。違うのは `\include` する schema 定義 file 名のみ。

### 実機検証 (ketoprofen amorphous、50 mol × 33 atom × 101 frame)

```
totalRecord = 101
Simulation_Conditions.Dynamics_Conditions.Time:
  delta_T              = 0.0205 [tau]  (= mdp dt=0.001 ps)
  Total_Steps          = 500000
  Output_Interval_Steps = 5000

各 record で:
- 各 frame の atom 座標 (1650 atom)
- Cell.{a, b, c} (≈ 2.63 nm cube)
- Energy.{Bond, Angle, ..., Total}.{Inst, Batch, Total} (kJ/mol)
- Temperature ([K]) / Pressure ([bar]) / Density ([kg/m^3]) / Volume ([nm^3])
  全 3 field embed
```

詳細仕様 (CLI option / xvg→UDF mapping / unit alias / cognac engine 対応 /
トラブルシューティング) は [gro2udf.md](gro2udf.md) 参照。

## PubChem 自動ダウンロード (`abmptools.amorphous.pubchem`)

`abmptools.amorphous.pubchem` モジュールは PubChem PUG REST API を叩いて、
3D SDF (MMFF94 最適化済み、水素込み) や canonical SMILES を取得します。
`requests` などの追加依存なし (`urllib` 標準ライブラリのみ)。

```python
from abmptools.amorphous import fetch_3d_sdf, fetch_smiles, PubChemNo3DError

# CID から
sdf_text = fetch_3d_sdf(3825)
smi      = fetch_smiles(3825)

# 名前から
sdf_text = fetch_3d_sdf("ketoprofen", by="name")

# 3D が無いケース
try:
    fetch_3d_sdf("some_polymer", by="name")
except PubChemNo3DError as e:
    print(e)   # → 見つからない旨とメッセージ
```

CLI:

```bash
# SDF を DL
python -m abmptools.amorphous.pubchem --cid 3825 -o ketoprofen.sdf
python -m abmptools.amorphous.pubchem --name aspirin -o aspirin.sdf
# SMILES だけ取得 (標準出力)
python -m abmptools.amorphous.pubchem --cid 3825 --smiles-only
```

戻り値 / 例外:

| ケース | 動作 |
|---|---|
| 3D SDF あり | SDF テキストを返す / ファイルに保存 |
| 3D なし (or CID 不在) | `PubChemNo3DError` を送出 (CLI は exit code 3) |
| CID 形式不正 (例: 桁超過) | HTTP 400 → `PubChemError` (CLI exit 1) |
| ネットワーク不通 | `PubChemError` (URLError 起因) |

`build_amorphous.py` からは `--pubchem_cid` / `--pubchem_name` で直接利用でき、
解決済み SDF は `<output_dir>/input/` にキャッシュされます (`--pubchem_cache_dir`
で場所を上書き可)。

## 同梱サンプル

| パス | 入力形式 | 説明 |
|---|---|---|
| [`sample/amorphous/`](../sample/amorphous/) | SMILES | pentane / benzene 二成分混合 (`run_sample.sh`) |
| [`sample/amorphous/ketoprofen/`](../sample/amorphous/ketoprofen/) | SMILES | ケトプロフェン 50 分子の詳細手順書 (`README.md`) |
| [`sample/amorphous/ketoprofen_pubchem/`](../sample/amorphous/ketoprofen_pubchem/) | SDF | PubChem 3D SDF (CID 3825) を `--mol` で読み込む例 (`run_sample.sh` + 入力 SDF 同梱) |

## テスト

amorphous ビルダーには 2 層のテストが用意されています:

| テスト | 件数 | 依存 |
|---|---|---|
| `tests/test_builder_mocked.py` | 8 | なし (標準ライブラリのみ)。OpenFF/Packmol/Interchange をすべて mock して `AmorphousBuilder.build()` の 6-stage フロー、返り値 dict、`config.json` 書き出しを検証。CI で常時走らせる層 |
| `tests/test_builder_integration.py` | 12 | OpenFF + Packmol + AmberTools (sqm) + RDKit。methane ×10, box 2 nm の小系を実際に build 、gro/top/ndx/5 MDP/2 scripts/config.json の構造・整合を spot check。依存が欠けると自動 skip (`@pytest.mark.slow`) |

部品単体のテスト (`test_amorphous_models.py` / `test_density.py` / `test_mdp_protocol.py` / `test_ndx_writer.py` / `test_pubchem.py`) と合わせて 50+ 件で amorphous 機能をカバー。

```bash
# 常時走る軽量テストだけ
pytest tests/test_builder_mocked.py tests/test_pubchem.py -v

# 実ビルドまで含めた統合テスト (abmptoolsenv など依存込みの環境で)
pytest tests/test_builder_integration.py -v

# slow だけ除外 / 実行
pytest -m "not slow" tests/
pytest -m slow tests/
```

## JSON 設定ファイル例

```json
{
  "components": [
    {"smiles": "CCCCC", "name": "pentane", "n_mol": 200},
    {"smiles": "c1ccccc1", "name": "benzene", "n_mol": 50}
  ],
  "density_g_cm3": 0.8,
  "temperature": 300,
  "T_high": 600,
  "seed": 42,
  "output_dir": "./pentane_benzene"
}
```

## 依存関係

| パッケージ | 用途 | 必須 |
|-----------|------|------|
| openff-toolkit | 分子準備・力場適用 | Yes |
| openff-interchange | GROMACS出力 | Yes |
| openmm | バックエンド | Yes |
| numpy | 数値計算 | Yes |
| packmol | 初期配置 | Yes (外部バイナリ) |
| ambertools | AM1-BCC電荷 | No (nagl代替可) |
| openff-nagl | ML電荷 | No (ambertools代替可) |
| rdkit | SDF読み込み | 推奨 |

全て lazy import（`import abmptools.amorphous` 自体は依存不要）。
