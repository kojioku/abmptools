# amorphous builder samples

`abmptools.amorphous` の使用例 4 種。どれも同じディレクトリ構成:

```
<sample_dir>/
├── run_sample.sh          # 実行スクリプト (入力は script 内 or JSON)
├── README.md              # (ketoprofen*/ のみ) 実行ログ / 補足
└── (実行後) input/ build/ md/  # 出力 (.gitignore 済み)
```

どのサンプルも `bash <sample_dir>/run_sample.sh` で build が回る。
詳しい解説は [../../docs/amorphous_tutorial.md](../../docs/amorphous_tutorial.md)。

## サンプル一覧

| サンプル | 入力 | 系 | 目的 |
|---|---|---|---|
| [`pentane_benzene/`](pentane_benzene/) | SMILES 2 成分 | pentane (200) + benzene (50), 0.8 g/cm³ | 最小の多成分サンプル。builder が通るか確認用 |
| [`ketoprofen/`](ketoprofen/) | SMILES 単成分 | ketoprofen (50), 0.8 g/cm³ | API 単成分系、RDKit SMILES→3D 経路 |
| [`ketoprofen_pubchem/`](ketoprofen_pubchem/) | PubChem CID → SDF | ketoprofen (50), 0.8 g/cm³ | 外部 3D conformer (MMFF94) 入力経路。SMILES 版との密度比較も可 |
| [`mixture_json/`](mixture_json/) | JSON (SMILES 内包) | pentane + benzene (JSON 編集で自由変更) | JSON schema のサンプル。複雑な多成分系はこちら |

## 実行の流れ (全サンプル共通)

```bash
cd sample/amorphous/<sample_name>
bash run_sample.sh

# 出力は ./input ./build ./md
# MD 実行 (GROMACS 必要)
cd md
bash run_all.sh
bash wrap_pbc.sh
```

## 環境要件

[`docs/amorphous_tutorial.md#1-環境構築`](../../docs/amorphous_tutorial.md)
参照。要点:

```bash
micromamba install -n abmptoolsenv -c conda-forge -y \
    openff-toolkit openff-interchange openmm rdkit packmol ambertools
export PATH=~/.local/share/mamba/envs/abmptoolsenv/bin:$PATH
```
