# abmptools.amorphous — Amorphous Structure Builder

API + 短鎖ポリマー（オリゴマー）混合系のアモルファス構造ビルダー。
OpenFF でパラメータ化 → Packmol で初期配置 → GROMACS アニーリング MD 用ファイル一式を生成。

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
    └── wrap_pbc.sh         # PBC ラップ後処理 (VMD 向け, 1.15.2+)
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
```

`wrap_pbc.sh` は `gmx trjconv -pbc mol -ur compact` を各 xtc と最終 gro に適用し、
箱境界で分断された分子を修復したうえで compact box 表示に揃えます。

## 同梱サンプル

| パス | 入力形式 | 説明 |
|---|---|---|
| [`sample/amorphous/`](../sample/amorphous/) | SMILES | pentane / benzene 二成分混合 (`run_sample.sh`) |
| [`sample/amorphous/ketoprofen/`](../sample/amorphous/ketoprofen/) | SMILES | ケトプロフェン 50 分子の詳細手順書 (`README.md`) |
| [`sample/amorphous/ketoprofen_pubchem/`](../sample/amorphous/ketoprofen_pubchem/) | SDF | PubChem 3D SDF (CID 3825) を `--mol` で読み込む例 (`run_sample.sh` + 入力 SDF 同梱) |

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
