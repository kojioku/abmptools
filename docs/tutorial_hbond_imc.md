# Tutorial: H-bond analysis on amorphous indomethacin (IMC)

非晶質インドメタシン (IMC) の MD trajectory に対して `abmptools.hbond` で
H-bond 検出 + 4-species 分類 (Yuan et al. 2015 NMR Table 1 と直接比較可能) +
可視化までを通すチュートリアル。

API / 設計の reference は [`hbond.md`](./hbond.md) を参照。NMR 帰属の背景は
論文 Yuan X. et al., *Mol. Pharmaceutics* 2015, 12, 4518–4528 (DOI
10.1021/acs.molpharmaceut.5b00705) を参照。

---

## 0. このチュートリアルで使うサンプルの構成

```
abmptools repo (軽量、~5 MB)
  sample/hbond/imc_amorphous/
    input/                           # ★ 入力ファイル (bundle 済み、本 tutorial の前提)
      IMC_result450.0_out_rec900.bdf # 1.5 MB、T=450 K, 125 mols, rec=900 単一スナップショット
      imc-bond-nmr.png               # Yuan 2015 Figure 5 (NMR 比較 plot 用)
    run_cli.sh                       # CLI ワンライナー
    run_notebook.ipynb               # Jupyter UI
    plot_nmr_comparison.py           # Yuan 2015 NMR vs MD 比較 plot
    README.md
    output/                          # 軽量出力 (CSV, .act, .py, png)
      imc_hbond_summary.csv
      imc_hbond_classification.csv
      imc_hbond_pairs.csv
      imc_hbond_show.act
      imc_hbond_show.py
      imc_hbond_count.png
      imc_hbond_nmr_comparison.png

abmptools-sample repo (重量 artifact、~5 MB、再生成可能)
  sample/hbond/imc_amorphous_20260516/
    output/
      imc_hbond.bdf                  # Mol_Name 維持 + Attributes (OCTA viewer pre-render)
      imc_hbond_action.bdf           # Action ヘッダ patch (gourmet autorun)
      imc_hbond_colored.bdf          # Mol_Name renamed (v1.25 legacy)
```

**入力ファイルは `sample/hbond/imc_amorphous/input/` に bundle 済み**なので、
本 tutorial の手順は abmptools repo を clone してすぐに走らせられる。別の
trajectory で解析したい場合は同じ `input/` 配下に置いて `BDF=...` 環境変数で
override する (詳細は 2 章)。

---

## 1. Environment

### 1.1 Python パッケージ

abmptools editable install 済みの conda / mamba 環境 (以下では `<myenv>` と
表記) を有効化:

```bash
mamba activate <myenv>            # または conda activate / venv activate
python -c "import abmptools.hbond; print('OK')"
```

**最小依存** (CLI ヘッドレス経路):

| パッケージ | 用途 | 入手 |
|---|---|---|
| `numpy` | 距離・角度計算 | abmptools 本体に含む |
| `matplotlib >= 3.5` | count vs record / NMR 比較 plot | `pip install abmptools[hbond]` |
| `UDFManager` | OCTA UDF / BDF I/O | **OCTA に同梱** (PyPI には無い) |

**Jupyter UI 経路を使うなら追加で**:

| パッケージ | 用途 |
|---|---|
| `ipywidgets >= 8.0` + `ipykernel` | Jupyter ipywidgets パネル |
| `rdkit-pypi >= 2022.09` | mol[0] 2D 構造プレビュー (carboxyl=赤 / amide=青 ハイライト) |

```bash
pip install "abmptools[hbond,jupyter,fragmenter]"
# fragmenter extras に rdkit-pypi が入る
```

### 1.2 UDFManager の入手 (OCTA / OCTA viewer)

UDFManager は PyPI に無く OCTA に同梱されている。Python から
`UDFManager` モジュールを通すには **OCTA の `python3/` ディレクトリ + UDF
定義 PATH** を通す必要がある。以下を `~/.bashrc` (または env の activate
script) に追記:

```bash
# OCTA85 の場合 (バージョンが違えば 85 を読み替え、install 先も適宜変更)
export OCTA85_HOME="$HOME/OCTA85"

# UDF 定義ファイル (engines 同梱の *.udf スキーマ)
export UDF_DEF_PATH="$OCTA85_HOME/ENGINES/udf:$UDF_DEF_PATH"

# Python module: GOURMET 側 (UDFManager 等) + ENGINES 側 (Cognac/COGNAC helpers 等)
export PYTHONPATH="$OCTA85_HOME/GOURMET/python3:$PYTHONPATH"
export PYTHONPATH="$OCTA85_HOME/ENGINES/python3:$PYTHONPATH"
```

設定後、別シェルを開きなおすか `source ~/.bashrc` で反映、以下で動作確認:

```bash
python -c "from UDFManager import UDFManager; print('UDFManager OK')"
```

env 起動時に自動で PATH が通っている運用もある。通っていない場合に
上記の 4 行を `~/.bashrc` に追加する。

### 1.3 可視化ツール (任意)

| ツール | 用途 |
|---|---|
| **OCTA gourmet** | `<prefix>_action.bdf` を開いて autorun overlay (推奨経路) |
| **OCTA viewer (GOURMET)** | `<prefix>.bdf` を開いて Python パネルで `<prefix>_show.py` を Load → Run、または Attribute フィルタで `hbond=Dual/Chain/Single/Free/Accept` 絞り込み |
| ブラウザ (Jupyter) | `run_notebook.ipynb` 経由で操作 (HTTP localhost) |

---

## 2. Smoke run: CLI 一発で全出力を生成

```bash
cd <abmptools>/sample/hbond/imc_amorphous

# CLI ワンライナー (約 10 秒) — bundle 済みの input/IMC_*.bdf を使う
bash run_cli.sh

# 別の BDF を解析したい場合
BDF=/path/to/my.bdf bash run_cli.sh
```

中身:

```bash
python -m abmptools.hbond input/IMC_result450.0_out_rec900.bdf \
    --out-prefix output/imc_hbond \
    --criteria luzar-chandler \
    --mol-name IMC \
    --colorize-mode both
```

実行ログ (要点):

```
Loaded input/IMC_result450.0_out_rec900.bdf
  125 molecules, 1 record(s)
  force field: GAFF2
  detected 125 carboxyls, 125 amides (125 tert), 0 N-H donors
  rec=0: COOH dual/chain/single/free=10/41/38/36 (8%/33%/30%/29%),
         amide accept/free=49/76 (39%/61%), hb_cc=31, hb_ca=50

Generated outputs:
  summary: output/imc_hbond_summary.csv
  classification: output/imc_hbond_classification.csv
  pairs: output/imc_hbond_pairs.csv
  uncolored: output/imc_hbond.bdf
  colored: output/imc_hbond_colored.bdf
  action_bdf: output/imc_hbond_action.bdf
  action_act: output/imc_hbond_show.act
  action_script: output/imc_hbond_show.py
  plot: output/imc_hbond_count.png
```

---

## 3. 出力ファイルの読み方

### 3.1 `summary.csv` (per-record の集計)

| column | 内容 |
|---|---|
| `record` | record index (0 始まり) |
| `n_carboxyls` / `n_amides` | 検出された官能基総数 (IMC は両方とも 125) |
| `n_carb_{dual,chain,single,free}` | per-COOH の 4 species カウント |
| `ratio_carb_*` | 同上の比率 (合計 1.0) |
| `n_amide_accept` / `n_amide_free` | per-amide の 2 species カウント |
| `n_hbonds_cc` / `n_hbonds_ca` | COOH-COOH / COOH-amide H-bond の総数 |
| `n_{dual,chain,single,free}_mols` | 分子代表 role の数 (色付け用 legacy) |

IMC ベースライン (T=450 K, rec=900):

| species | per-COOH | ratio |
|---|---|---|
| dual | 10 | 8.0% |
| chain | 41 | 32.8% |
| single | 38 | 30.4% |
| free | 36 | 28.8% |
| amide accept | 49 | 39.2% |
| amide free | 76 | 60.8% |

### 3.2 `classification.csv` (per-functional-group の役割表)

| column | 内容 |
|---|---|
| `record`, `group_type`, `mol_index`, `group_index` | 識別 |
| `role` | carboxyl: `dual` / `chain` / `single` / `free`; amide: `accept` / `free` |
| `partner_count` | 相手の数 |
| `partners` | `<mol>:<group_idx>;<mol>:<group_idx>;...` 形式 |

例: `mol_index=17` の COOH が `mol_index=42` の COOH と dual を組んでいる:

```
record,group_type,mol_index,group_index,role,partner_count,partners
0,carboxyl,17,0,dual,1,42:0
0,carboxyl,42,0,dual,1,17:0
```

### 3.3 `pairs.csv` (個別 H-bond)

```
record,kind,donor_mol,donor_d,donor_h,acceptor_mol,acceptor_a,d_da,d_ha,angle
0,cc,17,3,40,42,4,2.8123,1.7891,161.45
0,ca,0,3,40,47,2,2.8321,1.8308,163.50
...
```

`kind` は `cc` (COOH→COOH) または `ca` (COOH→amide)。`donor_d` は heavy donor
(OH oxygen)、`donor_h` は H、`acceptor_a` は heavy acceptor (=O)。距離は Å、
角度は度。

### 3.4 `count.png`

per-record の dual / chain / single / free mols の推移プロット (本 sample は
1 record なので点のみ表示)。

---

## 4. 可視化 (3 経路)

### 4.1 OCTA gourmet — autorun overlay (推奨)

```bash
gourmet output/imc_hbond_action.bdf
```

`autorun: showHbond()` が UDF オープン時に自動実行され、carboxyl atoms
(c/o/oh/ho) を **dual=red / chain=magenta / single=blue / free=gray**、amide
atoms (c/o/n) を **accept=cyan / free=gray** の sphere overlay で塗り分ける。

### 4.2 OCTA viewer (GOURMET) — Python panel (autorun crash 回避)

OCTA viewer (GOURMET) は autorun action 形式で crash することがある。その場合:

1. **File → Open** `output/imc_hbond.bdf` (Mol_Name 維持コピー、プリ描画 OK)
2. 左下 **Python panel** で `Load…` ボタンから `output/imc_hbond_show.py`
   を読み込む (or 中身をペースト)
3. **Run** ボタンで実行 → 4.1 と同じ overlay が描画

### 4.3 OCTA viewer Attribute フィルタ — 色なしカテゴリ可視化

`output/imc_hbond.bdf` の各 functional-group atom の `Attributes[]` に
`hbond=Dual/Chain/Single/Free/Accept` が append されているので、OCTA viewer の
Attribute フィルタで:

- `hbond=Dual` の atom のみ表示 → cyclic dimer の場所が分かる
- `hbond=Chain` の atom のみ表示 → disordered chain の場所
- `hbond=Accept` の atom のみ表示 → amide acceptor

等のカテゴリ別絞り込みが可能。

### 4.4 Mol_Name リネーム経路 (v1.25 legacy)

`output/imc_hbond_colored.bdf` は `Mol_Name` を `IMC_DUAL/CHAIN/SINGLE/FREE` に
書き換えて `Draw_Attributes.Molecule[]` に Red/Magenta/Blue/Gray を書く形式。
gourmet で開いてから Python panel の `show.all('line','mol','molname',...)`
を実行する必要がある (OCTA viewer プリ描画では空表示になるので注意)。

---

## 5. Jupyter UI で操作 (任意)

`sample/hbond/imc_amorphous/run_notebook.ipynb` は IMC 入力 BDF を bundle
した default で起動するが、**環境変数 `HBOND_BDF` を渡せば任意の BDF / UDF**
(PVA / peptide / アルコール / 混合系等) に切替できる。

### 5.1 起動コマンド (OS 別)

`abmptools` editable install 済みの env を有効化したシェルから:

#### Linux / macOS (bash / zsh)

```bash
cd <abmptools>/sample/hbond/imc_amorphous

# (a) default (bundle 済み IMC)
jupyter notebook --no-browser --port=8888 run_notebook.ipynb

# (b) 別 BDF / UDF に切替 (環境変数で 1 行 inline 指定)
HBOND_BDF=/path/to/pva.udf jupyter notebook --port=8888 run_notebook.ipynb

# (c) 環境変数を export 後に起動
export HBOND_BDF=/path/to/pva.udf
jupyter notebook --port=8888 run_notebook.ipynb
```

#### Windows PowerShell

```powershell
cd <abmptools>\sample\hbond\imc_amorphous

# (a) default
jupyter notebook --no-browser --port=8888 run_notebook.ipynb

# (b) 別 BDF / UDF に切替
$env:HBOND_BDF = "C:\path\to\pva.udf"
jupyter notebook --port=8888 run_notebook.ipynb

# 後で default に戻す
Remove-Item Env:HBOND_BDF
```

#### Windows cmd

```cmd
cd <abmptools>\sample\hbond\imc_amorphous

REM (a) default
jupyter notebook --no-browser --port=8888 run_notebook.ipynb

REM (b) 別 BDF / UDF に切替
set HBOND_BDF=C:\path\to\pva.udf
jupyter notebook --port=8888 run_notebook.ipynb

REM default に戻す
set HBOND_BDF=
```

(port 8888 が使用中なら `--port=8890` 等に変更)

### 5.2 パネル展開後の widget

ブラウザで notebook を開き、最初のセルを実行すると最終的に
`from abmptools.hbond import open_panel; panel = open_panel(BDF)` で
ipywidgets GUI パネルが展開:

| Widget | 内容 |
|---|---|
| info HTML | n_molecules / n_records / 検出 force field / 官能基カウント |
| **mol[0] 2D 構造図** | RDKit SVG (carboxyl=赤 / amide=青 ハイライト) |
| **Mode** dropdown | `imc` (4-species、本チュートリアル) / `generic` (PVA 等の任意系) |
| **mol prefix** text | default `IMC`、別系名 (`PVA`, `peptide` 等) に編集可 |
| **Donor groups** checkboxes | carboxyl / amide_donor / amine_donor / hydroxyl を複数選択可 |
| **Acceptor groups** checkboxes | carboxyl_O / amide_O / hydroxyl_O / ether_O を複数選択可 |
| 検出基準 dropdown | luzar-chandler / strict / custom (custom 時 d_DA / d_HA / 角度 slider 有効) |
| Lifetime checkbox + gap tolerance + dt | multi-record の時 |
| **Run** ボタン | 解析実行 → stats 表 + count plot inline 表示 |

### 5.3 別系での使い方例

PVA (アルコール、generic mode) を解析する場合:

```bash
# Linux
HBOND_BDF=/path/to/pva_traj.udf jupyter notebook run_notebook.ipynb
```

ブラウザのパネルで:

1. **Mode dropdown → `generic`** に切替
2. **mol prefix → `PVA`** に編集
3. **Donor groups → `hydroxyl` のみチェック** (COOH は無いので外す)
4. **Acceptor groups → `hydroxyl_O` のみチェック**
5. 検出基準は luzar-chandler のまま
6. **Run** をクリック

→ donor-type × acceptor-type pair の per-record 集計が表示される
(出力 CSV: `<prefix>_pair_stats.csv`)。

詳細は `sample/hbond/imc_amorphous/run_notebook.ipynb` の cell[0] markdown
にも書かれている。

---

## 6. NMR vs MD 比較 plot (Yuan 2015)

```bash
cd <abmptools>/sample/hbond/imc_amorphous
python plot_nmr_comparison.py
```

`output/imc_hbond_nmr_comparison.png` (3 段: Yuan Figure 5 画像 / Yuan Table 1
deconvolution bars / MD per-COOH bars、共通 ppm 軸) が生成される。

Yuan 2015 Table 1 (Mol. Pharm. 12, 4518) と MD の比較:

| ppm | species | NMR (Yuan) | MD (本 sample) | 差 |
|---|---|---|---|---|
| 179.3 | dual (cyclic dimer) | **58.5%** | 8.0% | **−50.5pt** |
| 176.3 | chain (disordered) | 15.2% | 32.8% | +17.6pt |
| 172.4 | single (COOH-amide) | 18.9% | 30.4% | +11.5pt |
| 170.4 | free | 7.5% | 28.8% | +21.3pt |

**観察**: MD で cyclic dimer 形成が著しく過小、chain と free が過大。
要因候補:

- T=450 K が高すぎ (Yuan 実験は室温〜200 °C、~373 K でも cyclic 優位)
- equilibration 不足 (cyclic dimer 形成は遅い過程の可能性)
- GAFF2 の COOH-COOH H-bond パラメータ精度
- 急冷の冷却速度が速すぎて dimer 形成が間に合わない

---

## 7. Trouble shooting

### 7.1 `ModuleNotFoundError: No module named 'UDFManager'`

section 1.2 の手順に従って OCTA の Python module PATH + UDF 定義 PATH を
通す (4 行とも `~/.bashrc` に追記):

```bash
export OCTA85_HOME="$HOME/OCTA85"
export UDF_DEF_PATH="$OCTA85_HOME/ENGINES/udf:$UDF_DEF_PATH"
export PYTHONPATH="$OCTA85_HOME/GOURMET/python3:$PYTHONPATH"
export PYTHONPATH="$OCTA85_HOME/ENGINES/python3:$PYTHONPATH"
```

### 7.2 Jupyter で mol[0] 2D 構造が "image" alt text しか表示されない

abmptools commit `23f6a88` 以降で fix 済み。古いバージョンを使っている場合
は editable install を更新:

```bash
cd <abmptools>
git pull
pip install -e .
```

それでも出る場合: Jupyter Kernel restart (`Kernel → Restart Kernel`) して
セル再実行。

### 7.3 Jupyter で count plot が "Image" alt text しか表示されない

同上、`6c0eb89` 以降で fix 済み。Kernel restart で反映。

### 7.4 OCTA viewer (GOURMET) で `<prefix>_action.bdf` を開いた瞬間 crash

OCTA viewer は autorun action 形式で落ちることがある。`<prefix>.bdf` (Mol_Name
維持コピー) を開いてから Python panel で `<prefix>_show.py` を Load → Run
する経路 (4.2) を使う。

### 7.5 OCTA viewer (GOURMET) プリ描画で分子が空表示

`<prefix>_colored.bdf` は Mol_Name が `IMC_DUAL/CHAIN/SINGLE/FREE` にリネーム
されており、OCTA viewer 内部の Mol_Name lookup が壊れる。`<prefix>.bdf` (Mol_Name
維持) を使う。

### 7.6 検出された carboxyls / amides が 0

入力 UDF の `Atom_Type_Name` が `MOL0_X` 等 per-atom unique で FF mapping が
hit しない場合 (OpenFF SMIRNOFF UDF の典型)。v1.28+ では default で element +
bond-graph fallback が自動 tag するので動くはず。古いバージョンの場合は
abmptools を更新するか、antechamber で GAFF type を assign して
`Atom_Type_Name` を書き換える ([`hbond.md`](./hbond.md) の 「Element + bond-graph
fallback」section 参照)。

### 7.7 `Cannot find \end{header}` エラー

gro2udf で生成した UDF (header section 無し) を `--colorize-mode action/both`
で処理した時に出る。v1.28+ で fallback が入っているので、abmptools を更新。

---

## 8. 次のステップ

- 別系での解析: COOH を持たない系 (PVA / アルコール / ペプチド主鎖)
  → `--classify-mode generic` で donor-type × acceptor-type pair 統計。
  例: `sample/amorphous/pva_amorphous/` を参照
- multi-record trajectory: `--gap-tolerance N --dt PS` で lifetime +
  Luzar-Chandler 自己相関 + τ_HB を出す (本 sample は 1 record のみなので
  lifetime は無効化されている)
- カスタム H-bond criteria: `--criteria custom --d-da-max 3.2 --angle-min 130`
  で閾値を変更
- API 直叩き: `from abmptools.hbond import Analyzer, AnalyzerConfig`、
  `AnalyzerConfig(bdf_path=..., classify_mode='imc', donor_groups=['carboxyl'],
  acceptor_groups=['carboxyl_O', 'amide_O'])` で組み立て、`run()` →
  `write_outputs()` 直叩きも可

---

## 関連ドキュメント

- [`hbond.md`](./hbond.md) — API / 設計 / CLI option 一覧
- [`architecture.md`](./architecture.md) — abmptools.hbond の位置づけ
- [Yuan et al. 2015](https://doi.org/10.1021/acs.molpharmaceut.5b00705)
  — *Mol. Pharm.* 12, 4518 (本 sample の NMR ref 元)
