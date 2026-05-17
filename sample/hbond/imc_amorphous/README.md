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

## 出力 (typical results, T=450 K)

| 量 | 値 |
|---|---|
| dual COOH-COOH 二量体分子数 | 10 (5 pair) |
| single COOH→amide 関与分子数 | 73 |
| free 分子数 | 42 |
| 総 H-bond 数 (cc / ca) | 31 / 50 |
| 検出基準 | Luzar-Chandler: d_DA ≤ 3.5 Å, ∠ ≥ 120° |

出力ファイル:
- `output/imc_hbond_summary.csv`: per-record count
- `output/imc_hbond_pairs.csv`: H-bond ペア一覧
- `output/imc_hbond_colored.bdf`: 3 色グループ化された UDF (gourmet で可視化)
- `output/imc_hbond_count.png`: count vs record (1 record のため点)

## gourmet での可視化

```bash
# OCTA gourmet 起動 (Windows or WSLg)
gourmet output/imc_hbond_colored.bdf
```

3 色グループ:
- **Red**: IMC_DUAL (環状二量体に参加)
- **Blue**: IMC_SINGLE (アミドに H-bond)
- **Gray + transparent**: IMC_FREE (H-bond 未関与)

参考: GOURMET の Draw_Attributes.Molecule[] の color は select 型 (Red/Green/Blue/Magenta/Cyan/Yellow/White/Black/Gray) のみ。RGBA 任意色は不可。
