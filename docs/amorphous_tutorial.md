# amorphous builder チュートリアル

`abmptools.amorphous` を使って非晶質構造を組み立て、GROMACS アニーリング MD
まで走らせる **hands-on ガイド**。CLI オプション一覧や API 仕様は
[amorphous.md](amorphous.md) を参照。

本チュートリアルを順にやると:

1. pentane/benzene 混合系サンプルを 1 コマンドで動かせる
2. 自分の分子を SMILES / SDF / PubChem CID から構築できる
3. 出力ディレクトリの構造と各ファイルの役割がわかる
4. MD を走らせて最終構造を得られる
5. つまずきポイントを回避できる

---

## 0. このチュートリアルで作るもの

```
<output_dir>/
├── input/    # 単分子 3D (PDB) + config.json (再現用)
├── build/    # Packmol 配置結果 + GROMACS 入力 (.top/.gro/.ndx)
└── md/       # 5-stage アニーリング MDP + run_all.sh
```

5-stage MD: EM → NVT-600K → NPT-600K → SA 600→300K → NPT-300K
(合計 ≈ 1.3 ns、小系なら数十分〜1時間)

---

## 1. 環境構築

### 1-1. micromamba 環境 `abmptoolsenv` に依存を入れる

```bash
micromamba install -n abmptoolsenv -c conda-forge -y \
    openff-toolkit openff-interchange openmm rdkit packmol ambertools
```

これで通常は動きます。もし Stage 4 で `pkg_resources` 関連の ImportError が
出たら §8-2 を参照 (setuptools の固定)。

### 1-2. abmptools をインストール

```bash
cd ~/llm-project/fcews-workspace/abmptools
~/.local/share/mamba/envs/abmptoolsenv/bin/pip install -e .
```

### 1-3. 実行時に PATH に通す

`packmol` バイナリを `subprocess` から見つけるため:

```bash
export PATH=~/.local/share/mamba/envs/abmptoolsenv/bin:$PATH
```

(or 実行時に `--packmol_path` で絶対パス指定)

### 1-4. 動作確認

```bash
python -c "from abmptools.amorphous.cli import main; print('OK')"
# → OK
```

---

## 2. チュートリアル A: 同梱サンプル (pentane/benzene 混合)

一番手軽。`sample/amorphous/run_sample.sh` が全パラメータ指定済み。

```bash
cd ~/llm-project/fcews-workspace/abmptools/sample/amorphous
bash run_sample.sh
```

内部コマンド (参考):

```bash
python ../../build_amorphous.py \
    --smiles "CCCCC" "c1ccccc1" \
    --name pentane benzene \
    --n_mol 200 50 \
    --density 0.8 \
    --temperature 300 \
    --seed 42 \
    --output_dir ./pentane_benzene \
    -v
```

**期待される出力** (末尾):

```
=== Build complete: ./pentane_benzene ===

Output files:
./pentane_benzene/build/mixture.pdb
./pentane_benzene/build/packmol.inp
./pentane_benzene/build/packmol.log
./pentane_benzene/build/system.gro
./pentane_benzene/build/system.ndx
./pentane_benzene/build/system.top
./pentane_benzene/input/component_0.pdb
./pentane_benzene/input/component_1.pdb
./pentane_benzene/input/config.json
./pentane_benzene/md/01_em.mdp
./pentane_benzene/md/02_nvt_highT.mdp
./pentane_benzene/md/03_npt_highT.mdp
./pentane_benzene/md/04_anneal.mdp
./pentane_benzene/md/05_npt_final.mdp
./pentane_benzene/md/run_all.sh
./pentane_benzene/md/wrap_pbc.sh
```

所要時間: 200+50 = 250 分子、AM1-BCC 電荷計算含めて **約 2-3 分**。

---

## 3. チュートリアル B: 自分の分子を SMILES から作る (ketoprofen)

単成分 50 分子、密度 0.8 g/cm³ の非晶質を作る例。

```bash
cd ~/llm-project/fcews-workspace/abmptools

python build_amorphous.py \
    --smiles "OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1" \
    --name ketoprofen \
    --n_mol 50 \
    --density 0.8 \
    --temperature 300 \
    --seed 42 \
    --output_dir ./my_ketoprofen \
    -v
```

**6 ステージの進行ログ (参考、所要約 80 秒)**:

| Stage | 内容 | 目安時間 |
|---|---|---|
| 1 | Molecule preparation (SMILES → 3D conformer, MW 計算) | 数秒 |
| 2 | Density / box estimation (立方体サイズ自動算出) | <1 秒 |
| 3 | Packmol packing (50 分子配置) | 約 6 秒 |
| 4 | OpenFF parameterization (AM1-BCC 電荷) | 約 60 秒 |
| 5 | NDX index groups | <1 秒 |
| 6 | MDP protocol 書き出し | <1 秒 |

### 確認ポイント

```bash
cat ./my_ketoprofen/input/config.json
# 再現用に全パラメータが保存されている
```

```bash
head ./my_ketoprofen/build/system.gro
# 1 行目: title、2 行目: 原子数 (50×33=1650)、末尾: box 3 値 (nm)
```

---

## 4. チュートリアル C: PubChem から 3D SDF を取って使う

SMILES から生成した conformer が不自然な場合、PubChem の実験由来 3D SDF を
使う方が安定。`sample/amorphous/ketoprofen_pubchem/run_sample.sh` が一式を
やってくれる。

```bash
cd ~/llm-project/fcews-workspace/abmptools/sample/amorphous/ketoprofen_pubchem
bash run_sample.sh
```

内部:

```bash
# 1. PubChem REST API で SDF をダウンロード (CID 3825 = ketoprofen)
curl -sL "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/3825/SDF?record_type=3d" \
     -o input/ketoprofen_pubchem_cid3825.sdf

# 2. SDF 入力で build
python ../../../build_amorphous.py \
    --mol input/ketoprofen_pubchem_cid3825.sdf \
    --name ketoprofen \
    --n_mol 50 --density 0.8 --temperature 300 --seed 42 \
    --output_dir . -v
```

Python API 経由で PubChem 取得したいなら:

```python
from abmptools.amorphous import fetch_3d_sdf, PubChemNo3DError

try:
    sdf_path = fetch_3d_sdf(cid=3825, outdir="input")
except PubChemNo3DError:
    print("3D conformer 無し、SMILES fallback を検討")
```

### SMILES 版との比較 (参考値)

| 系 | 最終密度 (g/cm³) | 初期配座ソース |
|---|---|---|
| ketoprofen (SMILES) | 1.139 | OpenFF 単 conformer |
| ketoprofen (PubChem 3D) | 1.141 | MMFF94-optimized |

実測 1.26 g/cm³ に対して -10 % 程度だが、OpenFF 2.1.0 unconstrained の
典型精度内。**初期配座ソースは最終密度にほぼ影響しない**。

---

## 5. チュートリアル D: JSON 設定で多成分系

パラメータが多い場合は JSON にまとめる方が楽。

`sample/amorphous/mixture.json`:

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
  "forcefield": "openff_unconstrained-2.1.0.offxml",
  "output_dir": "./pentane_benzene"
}
```

実行:

```bash
cd ~/llm-project/fcews-workspace/abmptools
python build_amorphous.py --config sample/amorphous/mixture.json
```

SDF 入力もサポート: `"sdf": "path/to/mol.sdf"` を components に書く。

---

## 6. 出力ディレクトリの歩き方

```
<output_dir>/
├── input/
│   ├── component_0.pdb       # 成分 0 の単分子 3D (H 込み)
│   ├── component_1.pdb       # 成分 1 (多成分時)
│   └── config.json           # 再現用全パラメータ
│
├── build/
│   ├── packmol.inp           # Packmol 入力
│   ├── packmol.log           # Packmol ログ (配置成否)
│   ├── mixture.pdb           # Packmol 出力: 全分子配置済み
│   ├── system.gro            # GROMACS 座標 (全原子)
│   ├── system.top            # GROMACS トポロジー
│   │                         #   → #include で ff/molecule.itp を展開
│   └── system.ndx            # インデックス (System 以外の group 定義)
│
└── md/
    ├── 01_em.mdp             # Energy minimization (50000 steps)
    ├── 02_nvt_highT.mdp      # NVT 600K (100000 steps = 100 ps)
    ├── 03_npt_highT.mdp      # NPT 600K 1bar (200000 steps = 200 ps)
    ├── 04_anneal.mdp         # Simulated annealing 600→300K (500 ps)
    ├── 05_npt_final.mdp      # NPT 300K 1bar (500 ps)
    ├── run_all.sh            # 全 5 stage を grompp + mdrun で連鎖実行
    └── wrap_pbc.sh           # gmx trjconv -pbc mol -ur compact で
                              # <stage>_pbc.xtc / 05_npt_final_pbc.gro 生成 (VMD 用)
```

主要ファイルの中身確認:

```bash
# config.json: パラメータ再現用
cat input/config.json

# system.gro: 2 行目に原子数、末尾にセル
head -2 build/system.gro; tail -1 build/system.gro

# packmol.log の末尾: 配置失敗 (tolerance 未満) の警告を確認
tail build/packmol.log
```

---

## 7. MD 実行 (GROMACS 必須)

GROMACS (`gmx` コマンド) がインストールされている前提。abmptoolsenv には
**入っていない** ので別途用意:

```bash
# 例: conda-forge 版
conda install -c conda-forge gromacs
```

### 7-1. 5-stage 連鎖実行

```bash
cd <output_dir>/md
bash run_all.sh
```

各 stage の開始時に進捗メッセージが出る。CPU 8 コアで小系 (200-500 原子) なら
5 stage 合計 ≈ 10-30 分、大系 (数千〜数万原子) なら数時間。

### 7-2. PBC アンラップ (VMD / ParaView で見る前に)

`05_npt_final.xtc` の生座標は分子がセル境界で断裂している。分子単位で
wrap し直す:

```bash
bash wrap_pbc.sh
# → 05_npt_final_pbc.xtc / 05_npt_final_pbc.gro が生成される
```

```bash
# VMD
vmd 05_npt_final_pbc.gro 05_npt_final_pbc.xtc
```

### 7-3. 最終密度の確認

```bash
gmx energy -f 05_npt_final.edr -o density.xvg
# "Density" を選択 → 平均値 (g/L) を 1000 で割ると g/cm^3
```

---

## 8. つまずきポイントと対処

### 8-1. `Fortran runtime error: Illegal seek` (Packmol 21.2.1)

**症状**: `abmptools.amorphous.packing.run_packmol` で Stage 3 が落ちる。

**原因**: Packmol 21.2.1 (conda-forge) の仕様。`subprocess.run(input=str)`
で渡す stdin パイプが seekable でないと Fortran IO が拒否する。

**対処**: abmptools 1.15.1+ で修正済み (`stdin=open(inp_path, "rb")` +
絶対パス化)。**古い版を使っている場合はアップデート**:

```bash
pip install --upgrade abmptools
# または feature branch: git pull && pip install -e .
```

---

### 8-2. `ImportError: cannot import name 'escape' from 'pkg_resources'`

**症状**: Stage 4 (OpenFF parameterization) で import エラー。

**原因**: setuptools 82+ は `pkg_resources` を分離したが、`openff.amber_ff_ports`
が未対応。

**対処**:

```bash
~/.local/share/mamba/envs/abmptoolsenv/bin/pip install "setuptools<81"
```

---

### 8-3. `FileNotFoundError: [Errno 2] No such file or directory: 'packmol'`

**症状**: Stage 3 で Packmol 起動に失敗。

**原因**: `PATH` に abmptoolsenv の bin が無い。Python インタープリタだけ
フルパス指定で起動すると起きる。

**対処** (2 通り):

```bash
# 1. PATH を通す
export PATH=~/.local/share/mamba/envs/abmptoolsenv/bin:$PATH

# 2. 絶対パス指定
python build_amorphous.py ... \
    --packmol_path ~/.local/share/mamba/envs/abmptoolsenv/bin/packmol
```

---

### 8-4. `PubChemNo3DError: CID XXX has no 3D conformer`

**症状**: `ketoprofen_pubchem/run_sample.sh` (や `fetch_3d_sdf`) が例外。

**原因**: PubChem の一部 CID は 3D 構造データを持たない。

**対処**: SMILES 版 (`--smiles "..." --name ...`) にフォールバック。または
別 CID (同一化合物の別 tautomer 等) を試す。

---

### 8-5. Packmol 配置失敗警告

**症状**: `packmol.log` 末尾に `tolerance not met` 的な警告。

**原因**: 密度が高すぎる or tolerance が厳しすぎて配置失敗。

**対処**:

- `--density` を下げる (0.8 → 0.6)
- Packmol は配置結果を出す (`mixture.pdb` は生成される) ので、そのまま MD で
  緩和すれば多くの場合解決
- どうしても駄目なら `sample/amorphous/run_sample.sh` の pentane/benzene
  (density 0.8) で試して環境が壊れていないかだけ確認

---

### 8-6. WSL2 で GPU が使えない

**症状**: GROMACS が CPU のみで遅い。

**原因**: WSL2 では NVIDIA OpenCL ICD が無い (conda-forge の GROMACS OpenCL
ビルドは GPU 無効)。

**対処**: 小系 (数百原子) なら CPU 8 コアで 1.3 ns ≈ 10 分。実用上 OK。
大系なら Linux native 環境で CUDA 版 GROMACS。

---

## 9. 次のステップ

- **CLI 詳細**: [amorphous.md](amorphous.md) (オプション全一覧、JSON schema、Python API)
- **結果の永続化**: プライベート [md-archive](https://github.com/kojioku/md-archive) に
  `<system>_<kind>_<yyyymmdd>_<key params>/` 命名で格納
- **fcews-manybody との連携**: GROMACS ルート (Phase 2c-D) では fcewsmb が
  本 builder を呼び出し、amorphous build → MD → contact 抽出 → χ 計算まで
  pipeline 化。詳細は fcews-manybody の
  [integration_test.md](../../fcews-manybody/docs/integration_test.md)
- **FAQ**: [faq.md](faq.md)

---

## Appendix A: 手計算での依存関係チェックリスト

```bash
# すべて通れば環境 OK
python -c "import openff.toolkit"                # OpenFF toolkit
python -c "import openff.interchange"            # Interchange
python -c "import openmm"                        # OpenMM (簡易 MD 用)
python -c "import rdkit"                         # RDKit
which packmol                                    # Packmol
which sqm                                        # AmberTools
python -c "import abmptools.amorphous.cli"       # abmptools 本体
```

## Appendix B: 成功例 (reference)

| 系 | 分子数 | 密度 | box (nm) | build 時間 | MD 時間 |
|---|---|---|---|---|---|
| pentane + benzene | 200 + 50 | 0.8 | 3.0 立方体 | ≈ 2 分 | ≈ 10-15 分 (CPU) |
| ketoprofen | 50 | 0.8 | 2.98 立方体 | ≈ 80 秒 | ≈ 10 分 (CPU) |

いずれも `sample/amorphous/` の記録より。
