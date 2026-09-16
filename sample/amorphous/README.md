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
python wrap_pbc.py
```

## 環境要件

[`docs/amorphous_tutorial.md#1-環境構築`](../../docs/amorphous_tutorial.md)
参照。要点:

```bash
micromamba install -n abmptoolsenv -c conda-forge -y \
    openff-toolkit openff-interchange openmm rdkit packmol ambertools
export PATH=~/.local/share/mamba/envs/abmptoolsenv/bin:$PATH
```

### Windows で組めない理由 (`.bat` を置いていない理由)

このサンプル群に `.bat` は無い。**Linux / macOS 専用**。他のサンプルが
`.sh` / `.bat` を対にしているのに、ここだけ `.sh` しかないのは意図的。

理由は **packmol** で、AM1-BCC ではない。

| 必要なもの | conda-forge win-64 | 備考 |
|---|---|---|
| `openff-toolkit` / `openff-nagl` | **あり** | 電荷は `--charge_method nagl` で解決済み |
| **`packmol`** | **無し** | 箱を詰める本体。代替経路なし |
| `ambertools` | 無し | `nagl` を使えば `sqm` は不要になる |

`--charge_method nagl` は「Windows でも動くように」入っている経路で、
`charge_from_molecules` を渡して **`sqm` を完全に回避する** (`sqm` に Windows
ビルドが無いため)。**そこは解決済み**。残るのが packmol で、
`abmptools.amorphous` は `shutil.which("packmol")` で探し、無ければ
`FileNotFoundError` で止まる。**詰める工程を飛ばす経路は無い。**

したがって Windows で組むなら WSL か Linux 機を使う。**組んだ後の解析は
Windows で動く** —— 出来上がった UDF / xtc を渡せば、`moldeck.hbond` を
Windows 側から `.bat` で回せる (`<moldeck>/sample/amorphous/` に対になった
スクリプトがある)。

> conda-forge の状況は 2026-09 に `micromamba search --platform win-64` で
> 確認した。`linux-64` では packmol / ambertools とも見つかるので、
> 検索自体は効いている。将来 win-64 ビルドが出たらこの節ごと見直すこと。
