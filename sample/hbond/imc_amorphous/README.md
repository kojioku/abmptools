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
- `output/imc_hbond.bdf`: Mol_Name 維持コピー (J-OCTA プリ描画用、色付けなし)
- `output/imc_hbond_colored.bdf`: Mol_Name リネーム + Draw_Attributes (3 色、`molname` モード)
- `output/imc_hbond_action.bdf` + `imc_hbond_show.act`: Mol_Name 維持 + autorun action (`action` モード、**OCTA gourmet 用**、官能基単位色付け)
- `output/imc_hbond_show.py`: autorun ラッパーなしの Python script (`action` モード、**J-OCTA Viewer の Python パネル用**、対象は `imc_hbond.bdf`)
- `output/imc_hbond_count.png`: count vs record (1 record のため点)

`--colorize-mode {molname,action,both}` で経路選択 (default = molname、サンプル
`run_cli.sh` では `both` を渡して全部出している)。

## 可視化方法

### (1a) Python action 経路 — OCTA gourmet 推奨 (`<prefix>_action.bdf` + `<prefix>_show.act`)

```bash
gourmet output/imc_hbond_action.bdf
```

`autorun: showHbond()` が UDF を開いた瞬間に走り、carboxyl atoms (c/o/oh/ho) を
役割色 (red/blue/gray) + amide atoms (c/o/n) を accept=cyan / free=gray の
sphere overlay で塗り分ける。**1 分子内で COOH と amide が別役割の場合も
両方が独立色で描画**される。Mol_Name は変えないので J-OCTA プリ描画とも互換。

### (1b) Python パネル経路 — J-OCTA Viewer 推奨 (`<prefix>.bdf` + `<prefix>_show.py`)

J-OCTA Viewer は (1a) の autorun action 形式で落ちることがある。その場合は:

1. J-OCTA Viewer で `output/imc_hbond.bdf` (Mol_Name 維持 uncolored copy) を開く
2. 左下 Python パネルで `Load…` ボタンから `output/imc_hbond_show.py` を読み込み
   (または .py の内容を直接ペースト)
3. `Run` ボタンで実行 → 同じ役割色 overlay が描画される

(1a) と同じ描画ロジックだが、autorun ラッパー無しで Python パネル直接実行形式。
BDF も無改変なのでファイルの round-trip 安全。

### (2) Mol_Name リネーム経路 (`<prefix>_colored.bdf`、v1.25 legacy)

```bash
gourmet output/imc_hbond_colored.bdf
```

左パネル Python タブで `show.all("line", "mol", ...)` の 3 番目引数を **"molname"** に
書換えて Run。**分子代表 role** (priority dual > single > free) で 3 色塗分け。

### (3) J-OCTA プリ描画

`output/imc_hbond.bdf` (Mol_Name 維持、色付けなし) または
`output/imc_hbond_action.bdf` (Mol_Name 維持 + action 色付き) を使う。
`_colored.bdf` (Mol_Name リネーム) は J-OCTA プリ描画で空表示になる。

3 色グループ:
- **Red**: IMC_DUAL (環状二量体に参加)
- **Blue**: IMC_SINGLE (アミドに H-bond)
- **Gray + transparent**: IMC_FREE (H-bond 未関与)

参考: GOURMET の Draw_Attributes.Molecule[] の color は select 型 (Red/Green/Blue/Magenta/Cyan/Yellow/White/Black/Gray) のみ。RGBA 任意色は不可。
