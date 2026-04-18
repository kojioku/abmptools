# ケトプロフェン非晶質構造ビルド手順

`abmptools.amorphous` によるケトプロフェン非晶質系の作成手順と実行結果の記録。

## 結果サマリ

- 分子: ケトプロフェン (C16H14O3, MW=254.28)
- SMILES: `OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1`
- 分子数: 50
- 目標密度: 0.8 g/cm^3
- ボックスサイズ: 2.9772 nm (立方体, 密度から自動)
- 力場: OpenFF `openff_unconstrained-2.1.0.offxml`
- 電荷: AM1-BCC (AmberTools)
- 乱数シード: 42
- 総原子数: 1650 (50 mol × 33 atom)

## 実行環境

micromamba 環境 `abmptoolsenv` にインストールしたパッケージ:

| パッケージ | バージョン | 入手元 |
|---|---|---|
| python | 3.10 | conda-forge |
| openff-toolkit | 0.16.10 | conda-forge |
| openff-interchange | 0.4.2 | conda-forge |
| openmm | 8.5.1 | conda-forge |
| rdkit | 2023.09.6 | conda-forge |
| packmol | 21.2.1 | conda-forge |
| ambertools | — | conda-forge (AM1-BCC用) |
| setuptools | <81 (80.10.2) | pip (openff.amber_ff_ports互換) |

## セットアップ

```bash
# 依存パッケージ一括インストール
micromamba install -n abmptoolsenv -c conda-forge -y \
    openff-toolkit openff-interchange openmm rdkit packmol ambertools

# setuptools をダウングレード (pkg_resources 互換性のため)
~/.local/share/mamba/envs/abmptoolsenv/bin/pip install "setuptools<81"
```

## 実行コマンド

```bash
cd /home/okuwaki/llm-project/fcews-workspace/abmptools/sample/amorphous

PATH=~/.local/share/mamba/envs/abmptoolsenv/bin:$PATH \
~/.local/share/mamba/envs/abmptoolsenv/bin/python3 ../../build_amorphous.py \
    --smiles "OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1" \
    --name ketoprofen \
    --n_mol 50 \
    --density 0.8 \
    --temperature 300 \
    --seed 42 \
    --output_dir ./ketoprofen \
    -v
```

※ `PATH` に abmptoolsenv の bin を含めるのは、build スクリプト内の `subprocess` 経由で `packmol` を見つけさせるため。

## 6ステージ実行ログ

| Stage | 内容 | 所要時間 |
|---|---|---|
| 1 | Molecule preparation (SMILES → 3D, MW計算) | 数秒 |
| 2 | Density / box estimation (自動: 2.9772 nm) | <1秒 |
| 3 | Packmol packing (50分子配置) | 約6秒 |
| 4 | OpenFF parameterization (AM1-BCC電荷計算) | 約1分 |
| 5 | NDX index groups | <1秒 |
| 6 | MDP protocol 書き出し | <1秒 |

全工程で約80秒。

## 出力ファイル

```
ketoprofen/
├── input/
│   ├── component_0.pdb       # 単分子 PDB
│   └── config.json           # 再現用設定
├── build/
│   ├── packmol.inp           # Packmol入力
│   ├── packmol.log           # Packmolログ
│   ├── mixture.pdb           # Packmol出力 (50分子配置済)
│   ├── system.gro            # GROMACS座標 (1650原子)
│   ├── system.top            # GROMACSトポロジー
│   └── system.ndx            # ケトプロフェン用インデックス
└── md/
    ├── 01_em.mdp             # EM (50000 steps)
    ├── 02_nvt_highT.mdp      # NVT 600K (100000 steps = 100 ps)
    ├── 03_npt_highT.mdp      # NPT 600K 1bar (200000 steps = 200 ps)
    ├── 04_anneal.mdp         # SA 600→300K (500000 steps = 500 ps)
    ├── 05_npt_final.mdp      # NPT 300K 1bar (500000 steps = 500 ps)
    └── run_all.sh            # GROMACS実行スクリプト
```

## MD 実行 (別途 GROMACS が必要)

```bash
cd ketoprofen/md
bash run_all.sh   # 01_em → 02_nvt → 03_npt → 04_anneal → 05_npt_final
```

GROMACS は abmptoolsenv に同梱されていないため別途インストールが必要。
5 ステージで合計 1.3 ns 相当のアニーリング MD。

## つまずきポイント (修正済)

### 1. Packmol stdin シークエラー

`packmol 21.2.1` (conda-forge) は `subprocess.run(..., input=str)` で標準入力に
文字列を渡すと `Fortran runtime error: Illegal seek` が発生する。
これは Python が作るパイプが seekable でないため。

**修正**: `abmptools/amorphous/packing.py` で `stdin=open(inp_path, 'rb')` に変更し、
pdb_paths と output_pdb を絶対パスへ解決してから Packmol input を書き出すようにした。

### 2. setuptools 82+ の pkg_resources 削除

`openff.amber_ff_ports` が `from pkg_resources import ...` を使っており、
setuptools 82 以降では pkg_resources が別パッケージ化されたため import エラーになる。

**対処**: `pip install "setuptools<81"` で 80.x に固定。
(conda-forge の openff-toolkit が最新 setuptools に追随するまでの暫定対応)

### 3. PATH に packmol が通っていない

`abmptoolsenv` を activate せず、python3 バイナリだけフルパス指定で起動すると、
packmol バイナリが PATH 上に無く FileNotFoundError となる。

**対処**: 実行時に `PATH=~/.local/share/mamba/envs/abmptoolsenv/bin:$PATH` を付与。
または `--packmol_path` で絶対パス指定でも可。
