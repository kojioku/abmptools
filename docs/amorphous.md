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
    └── run_all.sh          # GROMACS実行スクリプト
```

## MDP プロトコル（5段階）

| Stage | ファイル | アンサンブル | T | P | Steps | 目的 |
|-------|---------|-----------|---|---|-------|------|
| 1 | 01_em.mdp | - | - | - | 50000 | 急降下法 EM |
| 2 | 02_nvt_highT.mdp | NVT | 600K | - | 100000 | 高温 NVT (V-rescale) |
| 3 | 03_npt_highT.mdp | NPT | 600K | 1bar | 200000 | 高温 NPT (P-R) |
| 4 | 04_anneal.mdp | NPT | 600→300K | 1bar | 500000 | SA: annealing=single |
| 5 | 05_npt_final.mdp | NPT | 300K | 1bar | 500000 | 最終平衡化 |

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
