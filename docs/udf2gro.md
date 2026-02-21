# udf2gro — COGNAC-UDF to GROMACS Converter

## 概要

`abmptools.udf2gro` は COGNAC が生成した UDF ファイルを GROMACS 形式
（`.gro` / `.top` / `.mdp` / `.ndx`）に変換するパッケージです。

---

## モジュール構成

```
abmptools/abmptools/udf2gro/
├── __init__.py         # Exporter / SystemModel を公開
├── __main__.py         # python -m abmptools.udf2gro 用
├── cli.py              # CLI エントリ（main()）
├── exporter.py         # 薄いオーケストレータ
├── system_model.py     # dataclass 中間表現（Gromacs 文字列を持たない）
├── udf_adapter.py      # UDF → SystemModel 変換（Adapter 層）
└── gromacs/
    └── writers/
        ├── gro_writer.py  # .gro 出力
        ├── top_writer.py  # .top 出力
        ├── mdp_writer.py  # .mdp 出力
        └── itp_writer.py  # .itp（スタブ）
```

### 設計方針（Aプラン）

| 層 | クラス/モジュール | 責務 |
|---|---|---|
| Adapter | `UdfAdapter` | UDFManager → `SystemModel` |
| Model | `SystemModel` (dataclasses) | 純粋データ（Gromacs文字列なし） |
| Writer | `GroWriter`, `TopWriter`, `MdpWriter` | 文字列生成・ファイル書き出し |
| Exporter | `Exporter` | Adapter → Writers の薄い協調役 |
| CLI | `cli.py` | 引数解析・`Exporter` 呼び出し |

---

## 使い方

### CLI
```bash
python -m abmptools.udf2gro input.udf output_prefix
```

### Python API

```python
from abmptools.udf2gro import Exporter

Exporter().export("system.udf", "output")
# → output.gro, output.top, output.mdp が生成される
```

---

## テスト手順

`udf2gro/test/` にリファレンス入出力が用意されています。

```bash
cd sample/udf2gro/input

# 1. 変換実行（出力はカレントディレクトリ）
python -m abmptools.udf2gro test.udf test

# 2. 期待値と比較
diff test.top  ../output/test.top
diff test.gro  ../output/test.gro
diff test.mdp  ../output/test.mdp

# 3. 全て差分なしであれば PASS
```


---

## 対応力場・アンサンブル

### 力場（Force Field）

| FF_TYPE | 名称 |
|---|---|
| 1 | AMBER |
| 2 | GAFF |
| 3 | GAFF2 |
| 4 | DREIDING |
| 5 | OPLS |
| 6 | OPLSUA |

### アンサンブル

| UDF アルゴリズム | GROMACS 設定 |
|---|---|
| NVE | integrator=md-vv, tcoupl=no, pcoupl=no |
| NVT_Nose_Hoover | integrator=md-vv, tcoupl=nose-hoover |
| NVT_Berendsen | integrator=md-vv, tcoupl=berendsen |
| NPT_Andersen_Nose_Hoover | integrator=md-vv, tcoupl=nose-hoover, pcoupl=MTTK |
| NPT_Parrinello_Rahman_Nose_Hoover | integrator=md, tcoupl=nose-hoover, pcoupl=Parrinello-Rahman |
| NPT_Berendsen | integrator=md, tcoupl=berendsen, pcoupl=berendsen |
| Kremer_Grest | integrator=sd |

---

## 変更点サマリ

### 変更ファイル一覧

| ファイル | 変更内容 | 理由 |
|---|---|---|
| `udf2gro/udf2gro.py` | `main()` を `abmptools.udf2gro.Exporter` への委譲に変更 | CLI 互換を維持しつつ実装を移管 |
| `abmptools/abmptools/__init__.py` | `Udf2groExporter` の条件付きインポートを追加 | パッケージ統合（UDFManager 未インストール時は無視） |
| `abmptools/docs/udf2gro.md` | 本ドキュメントを新規作成 | テスト手順・設計説明 |

### 新規ファイル一覧

`abmptools/abmptools/udf2gro/` 以下の全ファイル（12ファイル）が新規作成。
元の `exportGromacs()` の責務を以下に分割：

- **UDF読み取り** → `udf_adapter.py`
- **GRO書き出し** → `gromacs/writers/gro_writer.py`
- **TOP書き出し** → `gromacs/writers/top_writer.py`
- **MDP書き出し** → `gromacs/writers/mdp_writer.py`
- **データ定義（各種 dataclass）** → `system_model.py`
- **協調・NDX書き出し** → `exporter.py`
