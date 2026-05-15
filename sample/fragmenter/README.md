# sample/fragmenter — abmptools.fragmenter 動作検証用サンプル

`abmptools.fragmenter` (v1.21.0+) の動作確認用に RDKit で生成した小規模 PDB。
`pdb2fmo` の入力としても流用できる構造で、`pip install -e '.[fragmenter]'` 後に
試せる。

## 内容

| ファイル | atoms | MW (g/mol) | 用途 |
|---|---|---|---|
| `pe_n20.pdb` | 62 (20 C + 42 H) | 282.6 | PE N=20 線形ポリエチレン。長い C-C 鎖の自動切断検証 |
| `pp_n10.pdb` | 92 (30 C + 62 H、主鎖 20 C + 分岐 10 C) | 422.8 | atactic polypropylene N=10。分岐含む鎖の主鎖検出 + 切断 |
| `propane5_acetone3.pdb` | 75 (combined: 5 propane + 2 acetone) | --- | 同一分子グループ化 + n_copies 展開検証 |
| `eude_block55.pdb` | 207 (90 heavy = MMA×5 + DMAEMA×5 block 共重合 + 117 H) | 1288.7 | 側鎖ありビニル系ポリマー (Eudragit E mimic)。MMA (ester 側鎖) + DMAEMA (amine 含む長い側鎖) block。主鎖 C-C 切断 + ester / amine 側鎖 + 多種 hetero atom 同居の検証 |

## 動作確認 (target_mw=200 default)

```bash
# PE N=20: 主鎖中央で 1 cut → 2 fragment
python -m abmptools.fragmenter suggest \
    --pdb sample/fragmenter/pe_n20.pdb \
    --output-dir ./review_pe \
    --target-mw 200
python -m abmptools.fragmenter apply \
    --pdb sample/fragmenter/pe_n20.pdb \
    --review-dir ./review_pe \
    --output ./review_pe/segment_data.dat

# PP N=10: 主鎖中央で 1 cut
python -m abmptools.fragmenter suggest \
    --pdb sample/fragmenter/pp_n10.pdb \
    --output-dir ./review_pp \
    --target-mw 200

# propane5_acetone3: 各 < 200 g/mol なので 0 cut、5+2=7 fragment 出力
python -m abmptools.fragmenter suggest \
    --pdb sample/fragmenter/propane5_acetone3.pdb \
    --output-dir ./review_mix \
    --target-mw 200
```

## 期待結果

| PDB | groups | n_copies | cut_sites/group | total fragments |
|---|---|---|---|---|
| pe_n20 | 1 | [1] | [1] | 2 |
| pp_n10 | 1 | [1] | [1] | 2 |
| propane5_acetone3 | 2 (CCC, CC(C)=O) | [5, 2] | [0, 0] | 7 |

## 補助検証 (リガンド付きタンパク質)

タンパク質 + リガンドが混在した PDB として、
`~/repos/abmptools-sample/sample/gbsa-genesis/input/3772L-rename.pdb` (FMODB データ、
CC-BY-SA、abmptools 同梱なし) を試せる。

```bash
python -m abmptools.fragmenter suggest \
    --pdb ~/repos/abmptools-sample/sample/gbsa-genesis/input/3772L-rename.pdb \
    --output-dir ./review_3772L \
    --target-mw 200
# -> 2 groups (protein chain, ligand) いずれも 0 cuts
#    (環内 / アミド / カルボニル 等で C-C 候補が排除される自然な挙動)
```

タンパク質残基ごとフラグメンテーションは **`abmptools.log2config`** で行うのが本筋
(本ツールは小分子・脂質・ポリマー専用)。

## サンプル PDB の生成方法

`pe_n20`, `pp_n10`, `propane5_acetone3` は RDKit から生成:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# PE N=20
mol = Chem.AddHs(Chem.MolFromSmiles("C" * 20))
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
Chem.MolToPDBFile(mol, "pe_n20.pdb")

# PP N=10 (atactic) — 主鎖 20 C + メチル分岐 10 C = 30 heavy + 62 H = 92 atoms
pp_smiles = "CC(C)" * 10
mol = Chem.AddHs(Chem.MolFromSmiles(pp_smiles))
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
Chem.MolToPDBFile(mol, "pp_n10.pdb")

# propane x5 + acetone x2
propane = Chem.AddHs(Chem.MolFromSmiles("CCC"))
acetone = Chem.AddHs(Chem.MolFromSmiles("CC(=O)C"))
AllChem.EmbedMolecule(propane, randomSeed=42)
AllChem.EmbedMolecule(acetone, randomSeed=43)
combined = propane
for _ in range(4):
    combined = Chem.CombineMols(combined, propane)
for _ in range(2):
    combined = Chem.CombineMols(combined, acetone)
Chem.MolToPDBFile(combined, "propane5_acetone3.pdb")

# EUD-E block 共重合体 (Eudragit E mimic): MMA×5 + DMAEMA×5
# 注: EUD-E は実際 random 共重合だが、simplified block で代表化学的特徴
# (ester + tertiary amine + 主鎖 C-C) を再現する
mma_unit = "CC(C)(C(=O)OC)"           # methyl methacrylate-like
dmaema_unit = "CC(C)(C(=O)OCCN(C)C)"  # DMAEMA-like (amine 側鎖)
eude_smiles = mma_unit * 5 + dmaema_unit * 5
mol = Chem.AddHs(Chem.MolFromSmiles(eude_smiles))
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
Chem.MolToPDBFile(mol, "eude_block55.pdb")
```

## 関連

- [`docs/fragmenter.md`](../../docs/fragmenter.md) — `abmptools.fragmenter` リファレンス
- [`abmptools/fragmenter/README.md`](../../abmptools/fragmenter/README.md) — モジュール構成早見表
- `tests/fragmenter/` — pytest による回帰テスト 14 件
