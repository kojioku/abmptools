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

### Pure-component pair (同名 component 2 個) の動作

`ComponentSpec(name="X")` を 2 個並べた場合 (例: A-B χ 計算で A=B にする
self-pair)、内部で以下が自動 dedup される:

- `mdp` の `tc-grps` は `"X"` (1 group) → `ref-t` / `tau-t` の単一値と整合
- `system.ndx` の `[ X ]` group は **両方の component の atom を集約** (片方
  だけ書かれて atom 半分が抜ける fail mode は v?? で修正済)
- `run_all.sh` は default で `MDRUN_OPTS="${MDRUN_OPTS:--ntmpi 1}"` を使う
  (DD 無効化)。凝縮系で box が cooling 中に縮んだ際の `box size ... is too
  small for a cut-off ... with 3 domain decomposition cells` を回避。
  override したい場合は `MDRUN_OPTS="-nt 8" bash run_all.sh` のように環境
  変数で渡す

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
