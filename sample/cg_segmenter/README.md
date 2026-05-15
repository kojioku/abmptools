# sample/cg_segmenter — `cg_segmenter` 動作検証用サンプル

`abmptools.fragmenter.cg_segmenter` (v1.24.0+) の検証用サンプル。

## 内容

| ファイル | atoms | MW (g/mol) | 用途 |
|---|---|---|---|
| `cholesterol_rdkit.pdb` | 74 (28 heavy = 27 C + 1 O + 46 H) | 386.7 | コレステロール (C27H46O、4 fused rings + iso-octyl tail + 3-OH)。fused ring の atom 共有検証 |

## 動作確認

```bash
python -m abmptools.fragmenter.cg_segmenter build \
    --pdb sample/cg_segmenter/cholesterol_rdkit.pdb \
    --output-dir ./cg_chol \
    --target-mw 200
```

期待結果:

| seg | kind | heavy | cap |
|---|---|---|---|
| 0 | ring_with_substituent | 6 | 3 H caps (3-OH を含む A ring) |
| 1 | ring | 6 | 5 H caps (B ring) |
| 2 | ring_with_substituent | 7 | 4 H caps (C ring + methyl) |
| 3 | ring_with_substituent | 7 | 3 H caps (D ring + methyl) |
| 4 | chain | 8 | 1 H cap (iso-octyl tail) |

合計 5 segments、約 99 atoms (cap 込み)。

**Shared atoms (6 pairs)**:

| 共有 atom | 属する segments |
|---|---|
| atom 11, 12 | seg 0 & seg 1 (A-B ring fusion) |
| atom 15, 16 | seg 1 & seg 2 (B-C ring fusion) |
| atom 19, 20 | seg 2 & seg 3 (C-D ring fusion) |

## サンプル生成方法 (本ディレクトリ内 PDB は以下で生成済)

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# cholesterol SMILES (canonical)
smi = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"
mol = Chem.AddHs(Chem.MolFromSmiles(smi))
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
Chem.MolToPDBFile(mol, "cholesterol_rdkit.pdb")
# atoms=74, heavy=28, MW=386.7, rings=4
```

## 別経路: PubChem CID 5997 から直接取得する場合

PubChem の cholesterol データ (CID 5997、CC-BY-NC ライセンス) を直接使うなら:

```bash
# SDF download (公式 REST API)
curl -o cholesterol_pubchem.sdf \
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5997/SDF?record_type=3d"

# SDF → PDB conversion (Open Babel 経由)
obabel -isdf cholesterol_pubchem.sdf -opdb -O cholesterol_pubchem.pdb
```

注: PubChem 由来ファイルは abmptools リポジトリには非同梱 (license の関係)。
本ディレクトリの `cholesterol_rdkit.pdb` は RDKit で SMILES から構築したもので、
同等の構造 (canonical SMILES が一致) を持つため検証用に十分。

## 関連

- [`docs/cg_segmenter.md`](../../docs/cg_segmenter.md) — `cg_segmenter` リファレンス
- [`abmptools/fragmenter/cg_segmenter/README.md`](../../abmptools/fragmenter/cg_segmenter/README.md) — モジュール構成
- `tests/cg_segmenter/` — pytest による回帰テスト 11 件
