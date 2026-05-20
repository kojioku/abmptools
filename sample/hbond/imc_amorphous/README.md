# Sample: IMC amorphous H-bond analysis

非晶質インドメタシン (IMC) の MD スナップショットに対する H-bond 検出/分類のサンプル。

## 入力

- BDF: `~/llm-project/SI/IMC_result450.0_out_rec900.bdf`
  - 125 分子 (GAFF2 atomtypes)、47.34 Å cubic box (initial)、record=900 単一スナップショット
  - T=450 K MD で平衡化された amorphous state

参考画像: `~/llm-project/SI/imc-bond-nmr.png`
13C CPMAS NMR で 3 種類のカルボニル環境 (~180/175/165 ppm) に分離する H-bond network。

## 実行

### CLI ワンライナー

```bash
bash run_cli.sh
```

### Jupyter notebook

```bash
jupyter notebook run_notebook.ipynb
```

## 出力 (typical results, T=450 K, per-functional-group, v1.27 候補)

| 量 | 値 |
|---|---|
| COOH dual (cyclic dimer 参加 COOH 数) | 10 (5 pair) |
| COOH single (amide に donate している COOH 数) | 49 |
| COOH free | 66 |
| amide accept (COOH から H-bond 受けている amide 数) | 49 |
| amide free | 76 |
| 総 H-bond 数 (cc / ca) | 31 / 50 |
| 検出基準 | Luzar-Chandler: d_DA ≤ 3.5 Å, ∠ ≥ 120° |

mol-level 代表 role (色付け用): 1 mol = 1 COOH なので per-COOH と同じ
(dual=10 / single=49 / free=66 mols)。

出力ファイル:
- `output/imc_hbond_summary.csv`: per-record の官能基単位集計 + mol-level legacy
- `output/imc_hbond_classification.csv`: 全 carboxyl / amide ごとの role table
- `output/imc_hbond_pairs.csv`: H-bond ペア一覧
- `output/imc_hbond.bdf`: Mol_Name 維持コピー (J-OCTA プリ描画用)
- `output/imc_hbond_colored.bdf`: 3 色グループ化された UDF (gourmet / OCTA ポスト描画用)
- `output/imc_hbond_count.png`: count vs record (1 record のため点)

## gourmet / J-OCTA ポスト描画での可視化

```bash
# OCTA gourmet 起動 (Windows or WSLg)
gourmet output/imc_hbond_colored.bdf
```

左パネル Python タブで `show.all("line", "mol", ...)` の 3 番目引数を **"molname"** に
書換えて Run。3 色グループに塗り分けられた分子描画が表示される。

J-OCTA プリ描画 (起動時の自動描画) は `output/imc_hbond.bdf` を使う (色付けなし、
分子形状確認のみ)。

3 色グループ:
- **Red**: IMC_DUAL (環状二量体に参加)
- **Blue**: IMC_SINGLE (アミドに H-bond)
- **Gray + transparent**: IMC_FREE (H-bond 未関与)

参考: GOURMET の Draw_Attributes.Molecule[] の color は select 型 (Red/Green/Blue/Magenta/Cyan/Yellow/White/Black/Gray) のみ。RGBA 任意色は不可。
