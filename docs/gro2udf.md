# gro2udf — GROMACS .gro to COGNAC-UDF Converter

## 概要

`abmptools.gro2udf` は GROMACS が生成した `.gro` ファイル（位置・速度・セル情報）を
COGNAC-UDF の Structure レコードに書き戻すパッケージです。

逆方向（UDF → GRO）の変換は `abmptools.udf2gro` が担います。

### まず動かしてみたい人へ

**入力が一式そろったサンプルが同梱してあります。** 読む前に 1 回通すのが
早いです。

```bash
cd sample/gro2udf
bash run.sh
```

3 本走ります —— 座標の差し替え、`.top` からの新規作成、そして
**軌跡 + エネルギー込み**の変換。`.edr` と `.xtc` も同梱してあるので、
エネルギーを入れた変換もそのまま試せます。詳しくは
[`sample/gro2udf/README.md`](../sample/gro2udf/README.md)。

自分の系で打つコマンドは **§ 使い方** にあります。

---

## 2 つの使い方 — どちらを選ぶか

**「UDF が既にあるか」で決まります。**

| | **udf-and-gro モード** | **`--from-top` モード** |
|---|---|---|
| 何をするか | **既にある UDF の座標だけを差し替える** | **UDF を新しく作る** |
| 入力 | `.udf` ＋ `.gro` | `.top` ＋ `.gro` |
| UDF の分子構造・力場 | **元の UDF のものをそのまま使う** | `.top` から組み立てる |
| 使う場面 | OCTA で系を組んだ後、GROMACS で MD を回して**構造を戻したい** | GROMACS で組んだ系を**初めて OCTA に持ち込む** |
| コマンド | `python -m abmptools.gro2udf in.udf out.gro` | `python -m abmptools.gro2udf --from-top system.top out.gro` |

### udf-and-gro モード — 座標の差し替え

もとの UDF が持っている **Set_of_Molecules（分子構造）と Interactions（力場）は
一切触らず**、Structure レコードの座標・速度・セルだけを `.gro` の値で上書きします。

```
既存 UDF  ─┬─ 分子構造・力場  ──────────────→  そのまま
           └─ 座標・速度・セル  ←── .gro から差し替え
```

**OCTA で系を組み、GROMACS で MD を回し、結果を OCTA に戻す**、という
往復のうち「戻り」にあたります。UDF 側に何も足せないので、GROMACS 側で
分子を増減させた場合は使えません。

> **実測 (2026-09-07)**: IMC + PVP の 30 分子系で確認。`Set_of_Molecules` の
> ブロックは 726,489 文字が**バイト単位で一致**し、先頭原子の座標だけが
> `[9.67, 24.80, 1.15]` → `[9.66, 24.82, 1.00]` に置き換わった。

### `--from-top` モード — UDF を新規作成

GROMACS の topology（`.top` / `.itp`）を読んで、**分子構造・結合・力場・電荷を
すべて UDF に書き起こします**。`.gro` からは座標とセルを取ります。

```
.top  ──→ 分子構造・結合・角度・二面角・LJ・電荷
.gro  ──→ 座標・セル                              ─→ 新しい UDF
.mdp  ──→ 温度・カットオフ（任意、下記）
```

元になる UDF が要らないので、**GROMACS だけで組んだ系を OCTA に持ち込める**
のが利点です。UDF のスキーマ（どの項目を持つか）はテンプレートから取るので、
`--template` を省略すると同梱の `default_template.udf` が使われます。

`--mdp` を渡すと、温度・`tau_t`・`rcoulomb` から **Nose-Hoover の Q と Ewald の
カットオフを自動計算**します。省略した場合はテンプレートの値のままです。

> **どちらを使うか迷ったら**
> 手元に「元になる UDF」があるなら udf-and-gro、無いなら `--from-top`。

### 入力はどこから来るか

`--from-top` に渡す `.top` / `.gro` は、多くの場合 **`abmptools.amorphous` で
組んだ非晶質セルを GROMACS で回した結果**です。

```
abmptools.amorphous  ─→  build/system.top    ─┐
                         build/system.gro     │
GROMACS (5 段階 MD)  ─→  md/05_npt_final.gro ─┼─→  gro2udf --from-top  ─→  UDF
                         md/05_npt_final.mdp  │
                         md/..._energy.xvg   ─┘
```

構築から MD までの手順は
[amorphous_tutorial.md](amorphous_tutorial.md)、amorphous からの標準的な変換
フロー (trajectory と energy を一緒に入れる例) は
[amorphous.md の「OCTA UDF/BDF への変換」](amorphous.md) にあります。

---

## モジュール構成

### udf-and-gro モードの実装

```
abmptools/abmptools/gro2udf/
├── __init__.py        # Exporter / TopExporter を公開
├── __main__.py        # python -m abmptools.gro2udf 用
├── cli.py             # CLI エントリ（main()）
├── exporter.py        # 薄いオーケストレータ
├── gro_parser.py      # GRO ファイル → GROFrame（Parser 層）
├── gro_adapter.py     # GROFrame → AtomPosition / CellGeometry（Adapter 層）
└── udf_writer.py      # AtomPosition / CellGeometry → UDF レコード（Writer 層）
```

### `--from-top` モードの実装

```
abmptools/abmptools/gro2udf/
├── top_parser.py          # TOP/ITP → TopRawData（Parser 層）
├── top_adapter.py         # TopRawData + GRO → TopModel（Adapter 層）
├── top_model.py           # 中間表現 dataclass 群
├── top_exporter.py        # TopModel → COGNAC UDF（Writer 層）
├── mdp_parser.py          # MDP → MdpParams（シミュレーション条件）
├── guard.py               # 変換できない .top を検出して止める
└── default_template.udf   # テンプレート未指定時のデフォルト UDF
```

### ★ 変換できない `.top` は通さない（`guard.py`）

`--from-top` が読むセクションと、書ける COGNAC ポテンシャル型は固定である。
それ以外は**メッセージ無しで捨てられていた**ので、変換は成功したように見えて
中身が静かに間違った UDF が出ていた。

きっかけは Martini 2 の topology である。非結合パラメータが全部
`[ nonbond_params ]` にあり `[ atomtypes ]` は `c6 = c12 = 0` なので、
**`Pair_Interaction` が 0 本の UDF が無警告で生成されていた**
（`Interaction_Site_Type[].Range` は `sigma * 1.5` なのでこれも 0 になる）。

重大度は 2 つに分かれる。**復旧可能性が違うため**である。

| 重大度 | 何が起きるか | 既定の挙動 |
|---|---|---|
| **fatal** | 項が**誤った関数形・誤った数値で書かれる**。UDF は完成品に見えるので下流で気づけない | **エラーで停止** |
| warning | 項が**落ちる**。不在として見えるうえ、全原子 topology はずっとこの挙動だった | 警告のみ。変換は続行 |

**fatal**

| 検出 | 失われるもの |
|---|---|
| `[ defaults ]` comb-rule 1 | `[ atomtypes ]` 末尾 2 列は `c6`/`c12` だが σ/ε として読まれる。全 LJ が誤り |
| `[ nonbond_params ]` に行がある | 明示された非結合ペア全部。Martini はここに全 LJ を置く |
| `[ constraints ]` に行がある | 拘束結合。力の定数を持たないので、原子が何にも繋がらなくなる |
| bond funct ≠ 1 | funct に関わらず COGNAC `Harmonic` として書かれる |
| angle funct ≠ 1 | funct に関わらず COGNAC `Theta` として書かれる。Martini の G96 角度 (funct 2) は**力の定数の単位すら違う** (kJ/mol vs kJ/mol/rad²) |

**warning**

`[ pairtypes ]` / `[ settles ]` / `[ exclusions ]` / `[ virtual_sites* ]` /
`[ dummies* ]` / `[ cmaptypes ]` / `[ cmap ]`、および dihedral の funct が
1 / 3 / 4 / 9 以外の場合。

`--allow-unsupported` を付けると fatal も警告に落として変換する。

```bash
python -m abmptools.gro2udf --from-top martini.top conf.gro --out out.udf
# UnsupportedTopFeatureError: this .top uses features gro2udf writes incorrectly:
#   - [ defaults ] declares comb-rule 1, ...
#   - [ nonbond_params ] has 10 entries but is not read; ...
#   - angle funct 2 is written as COGNAC Theta regardless ...
```

#### 検出しないもの・その理由

- **`[ pairs ]`** は損失ではない。COGNAC は `Scale_1_4_Pair`（`fudgeLJ` から設定
  済み）で 1-4 を再現する。`[ pairtypes ]` による型ごとの上書きだけが、この
  単一のスケール係数では表せないので警告する
- **`[ position_restraints ]` 等の restraints 系**は挙げていない。`#ifdef` は
  gro2udf のどこでも評価されないため、`#ifdef POSRES` の中にある拘束を「有効」
  として警告すると、当たるより外れる方が多くなる。
  `[ settles ]` は `#else` 側＝実際に有効な枝なので残してある
- **dihedral の funct は生の行から数えている**。非対応 funct の行は parser の
  「型参照のみ」分岐に落ち、**funct が torsion 型の番号として保存される**ため、
  パース後には funct が残らない（しかも存在しない型を指す）

**guard は検出しかしない。** これらを実際に変換できるようにするのは別の作業である。

### 設計方針

変換パイプラインは以下の 3 段階:

```
GROParser  →  GROFrame
GROAdapter →  List[AtomPosition], CellGeometry   ← SystemModel の dataclass を再利用
UDFWriter  →  UDF レコード（Steps, Time, Position, Velocity, Unit_Cell）
Exporter   →  上記の薄い協調役（フィルタ・ループ管理）
```

| 層 | クラス/モジュール | 責務 |
|---|---|---|
| Parser  | `GROParser`  | GRO ファイルをフレーム単位で読み込む |
| Adapter | `GROAdapter` | GROFrame → `AtomPosition` / `CellGeometry` |
| Writer  | `UDFWriter`  | UDF に Steps/Time/Position/Velocity/Cell を書く |
| Exporter| `Exporter`   | Parser→Adapter→Writer を協調させる |
| CLI     | `cli.py`     | 引数解析・`Exporter` 呼び出し |

**中間表現**: `abmptools.udf2gro.system_model` の `AtomPosition` / `CellGeometry`
を再利用します（新たな中間表現クラスは作りません）。

`gro2udf` 文脈での `AtomPosition` フィールドの解釈:
- `mol_id`  → 0-based UDF mol インデックス（`Structure.Position.mol[m]` の m）
- `atom_id` → 0-based atom インデックス（`atom[a]` の a）
- `vx/vy/vz` → [nm/ps]（UDFWriter が × 1000 して [m/s] で put）

### ★ 座標が入るのは動的レコードだけ（静的 Structure はテンプレート由来）

UDF は **静的セクション（共通部）** と **動的レコード（フレームごと）** に分かれる。
`gro2udf` が静的セクションへ書くのは次の 4 つだけである。

| 静的セクションに書くもの | 書かないもの |
|---|---|
| `Set_of_Molecules` | **`Structure`**（座標・セル）|
| `Simulation_Conditions` | **`Initial_Structure`** |
| `Molecular_Attributes` | |
| `Interactions` | |

**座標とセルは `newRecord()` した先、つまり動的レコードにしか書かれない。**
静的 `Structure` と `Initial_Structure` は `--template` に渡した UDF の内容が
**そのまま残る**。

同梱テンプレート（`default_template.udf` / `default_template_cognac101.udf`）は
この部分が空なので、何も主張しない。問題になるのは **実在の UDF を
`--template` に渡したとき**で、COGNAC の系を往復させるときは自然に
そうなる。このとき出来上がる UDF は:

- 静的 `Structure` … **変換前**の座標と箱（テンプレートのもの）
- 動的レコード … **変換後**の座標と箱（GROMACS から読んだもの）

の 2 つを同時に持つ。**どちらも形式としては正しいので、`gmx` も COGNAC も
OCTA viewer も文句を言わない。** 静的側を初期構造として読む経路に乗せると、
変換したはずの構造ではなく元の構造で走る。

テンプレートの静的 `Structure` が空でない場合は **変換時に警告を出す**
（変換自体は止めない）:

```
template foo.udf already holds positions for 1024 molecule(s) in its static
Structure block. gro2udf writes coordinates and the cell into dynamic records
only, so that block is carried over unchanged and keeps the template's
structure and box. ...
```

静的側にも変換後の構造を反映したい場合は、後処理で書き込むか、空のテンプレート
（同梱の `default_template.udf`）を使って `Set_of_Molecules` 以下を作り直す。

---

## ★ 下流の GROMACS 変換器へ渡すときに要るもの

`gro2udf` の UDF を下流の GROMACS コンバータに通す
場合、**UDF が形式として正しいだけでは足りない**。下流は UDF スキーマが
「空でも妥当」としている場所を規約として読んでいる。**変換は成功し、下流の
検証も通り、NVT では何も起きず、NPT でだけ落ちる**という形になるので、
気付きにくい (2026-09 に実機で 5 件見つかった)。

### Nose-Hoover の `Q` と `Cell_Mass`

`Q` は Nose-Hoover の**熱浴質量**で、GROMACS の `tau_t` に対応する量。
両者の関係はどちらのマニュアルにも書かれている。

```
COGNAC マニュアル 式 2.6      dζ/dt = ( Σ pᵢ²/mᵢ − g·k_B·T ) / Q
GROMACS リファレンス 式 51    Q = τ_T² · N_f · k · T₀ / (4π²)
  → τ_T について解くと         tau_t = 2π · √( Q_d / (g·k_B·T) )
```

`Q` は `g·k_B·T` と組で現れるので、**時間に直すには `g·k_B·T` で割る**。
`Q_d` は `Q · Unit_Parameter.Mass · Unit_Parameter.Length²`。詳細は
[udf2gro.md](udf2gro.md)。

`Q` は **アルゴリズムごとに別のフィールド**に置かれる。`NVT_Nose_Hoover.Q`
だけ書いて NPT 側を空にすると、**NPT に切り替えた瞬間に `Q = 0` が読まれ `tau_t = 0` に
なり、Nose-Hoover が 0 除算して箱が `nan` に飛ぶ**。`gro2udf` は NVT / NPT の
両系統に同じ `Q` を書く。

`Cell_Mass` (バロスタット質量) も同様で、未記入だと `tau_p` の式が 0 になり
下限の 2.0 に丸められる。`CognacSystemUtil.setCellMass` は
**系の全質量**を入れるので、`gro2udf` もそれに合わせる。

### ★ `--nh-dof` — 自由度の数え方

`Q = g * k_B * T * tau^2` の `g` をどう数えるか。**既定は `3N-3`**。

| | 何が変わるか |
|---|---|
| **`3N-3`** (既定) | `Q = (3N-3)·k_B·T·τ²` / mdp は **`comm-mode = Linear`** (`nstcomm = 100`) |
| `3N` | `Q = 3N·k_B·T·τ²` / mdp は **`comm-mode = None`**。**既存の UDF と同一**になる |

**既定を `3N-3` にしている理由。** GROMACS 自身の既定が `comm-mode = Linear`
で、重心のドリフトを除くのが MD の通常の作法。`None` のままだと grompp が
毎回「運動エネルギーが重心に溜まる」と警告する。出力の行き先は GROMACS
なので、その流儀に合わせている。**既存の UDF と揃えたいときは
`--nh-dof 3N`**。

**`comm-mode` も一緒に切り替わる。** `Q` の数え方だけ変えて mdp が
`comm-mode = None` のままだと、GROMACS は 3N で積分してしまい**オプションが
説明どおりに動かない**。`gro2udf` は
`Simulation_Conditions.Dynamics_Conditions.Moment` の
`Calc_Moment` / `Stop_Translation` を立てて下流に `Linear` を書かせる。

確認は GROMACS の報告で:

```
--nh-dof 3N     Number of degrees of freedom in T-Coupling group System is 9150.00
--nh-dof 3N-3   Number of degrees of freedom in T-Coupling group System is 9147.00
```

```bash
python -m abmptools.gro2udf --from-top system.top md.gro --nh-dof 3N-3 --out out.udf
```

根拠は実測。既存の UDF 4 系 (原子数 23 / 80 / 110 / 3050) の `Q` は
いずれも `3N * k_B * T * (0.1 ps)^2` にぴったり乗る (逆算した tau が 4 系とも
0.099999972 ps)。**`3N-3` では乗らない** —— 80 原子で 0.1006 とずれる。
GROMACS 自身もその mdp に対し `degrees of freedom ... is 9150.00`
(= 3 x 3050) と報告した。

差は **3050 原子で 0.03%、80 原子で 0.6%** 程度。`Q` は熱浴の応答の速さを
決めるだけでサンプリングされるアンサンブルは変えないので、**大きい系では
どちらでも実害はない**。小さい系ほど効く。

> **`Q` を書いたら、戻した先の `tau_t` を確かめる。** `Q` は自由度と `k_B` を
> 含む熱浴質量なので、`tau_t` に直すには `g·k_B·T` で割る必要がある
> ([udf2gro.md の「`Q` → `tau_t` の変換」](udf2gro.md))。この系の `Q` なら
> `tau_t = 0.628 ps` に対応する。
>
> `.mdp` は雛形と考え、**本番では `tau_t` / `tau_p` を実験条件に合わせて
> 見直すこと。** `tau_p >= 2 * tau_t` (共鳴を避ける) も確認する。

### 力場 ID と 1-4 スケーリング

- **`Unit_Parameter.Comment` の `FF=n`** — 下流はここで力場を判定する。
  空だと「力場が分からない」扱いになり書き出しが通らない。`--ff` で指定
  (既定 `gaff` = `FF=2`)。番号は力場定義ファイルの並び
- **`User_Torsion.Parameters[]` の `SCNB` / `SCEE`** — 1-4 のスケーリング。
  空だと下流が 1-4 の扱いを決められず NPT が流せない。`.top` の
  `fudgeLJ` / `fudgeQQ` をそのまま入れる

### 静的セルとテンプレート

- 静的 `Structure.Unit_Cell` と `Initial_Structure.Initial_Unit_Cell` には
  **最初のフレームの箱**を書く。書かないとテンプレートの値 (同梱テンプレート
  なら 20 A 立方と 100 A 立方) が残り、座標と箱が別の系の UDF になる
- **`<top と同名>.udf` を自動でテンプレートに採用しない。** 隣にあるのは
  たいてい MD 前の UDF で、その箱を引き継いでしまう。使うなら `--template` で
  明示する

## 使い方

### ★ どれを打てばよいか — 3 段階

想定は「GROMACS で回し終えて、その結果を OCTA viewer (GOURMET) で見たい」人。
**下に行くほど 1 つずつ足すだけ**で、前の段を打ち直す必要はありません。

#### 1. 構造だけ (gmx は要らない)

```bash
python -m abmptools.gro2udf --from-top system.top md/prod.gro \
    --mdp md/prod.mdp \
    --out prod.udf
```

`--mdp` は任意ですが、**渡すと Nose-Hoover の `Q` と Ewald カットオフが
その `.mdp` の値から決まります**。渡さないと既定値が入るので、下流で
GROMACS に戻す予定があるなら渡してください (**§ Nose-Hoover 熱浴質量 Q**)。

#### 2. + 軌跡 (gmx が要る)

```bash
python -m abmptools.gro2udf --from-top system.top md/prod.gro \
    --mdp md/prod.mdp \
    --trajectory md/prod.xtc \
    --max-frames 100 \
    --out prod.udf
```

- `--trajectory` は **生の `.xtc` をそのまま渡して構いません**。`gro2udf` が
  先に `gmx trjconv -pbc nojump` を通すので、周期境界で割れた分子が
  そのまま UDF に入ることはありません (下の **§ `--trajectory` は gmx を呼ぶ**)
- **`--max-frames` は付けてください。** 付けないと**全フレーム**が入ります
  (下の **§ `--max-frames` と `--frame-step`**)
- `.xtc` を読むには `MDAnalysis` が要ります。多フレーム `.gro` なら不要

#### 3. + エネルギー (gmx が要る)

```bash
python -m abmptools.gro2udf --from-top system.top md/prod.gro \
    --mdp md/prod.mdp \
    --trajectory md/prod.xtc \
    --max-frames 100 \
    --edr md/prod.edr \
    --out prod.udf
```

`--edr` を渡すと `gmx energy` がその場で走り、全エネルギー項が各フレームの
`Statistics_Data` に入ります。OCTA viewer で開けば**軌跡とエネルギープロットが
同じ UDF から再生されます**。

> **`.edr` は軌跡と同じ run のものを使ってください。** エネルギーは
> **フレームの時刻に最も近い行**を割り当てる仕組みなので、別の run の `.edr`
> を当てても**それらしい値が入り、出力を見ても気付けません**。

---

### `--from-top` のオプション

| オプション | 既定 | 何をするか | gmx |
|---|---|---|---|
| `--out PATH` | `<gro の stem>_fromtop.udf` | 出力先 | |
| `--mdp PATH` | 使わない | `ref_t` / `tau_t` / `rcoulomb` から `Q` と Ewald を決める | |
| `--template PATH` | 同梱テンプレート | スキーマの元にする UDF | |
| `--cognac-version N` | `112` | テンプレートの `cognac<N>.udf` include を差し替える (OCTA8.4 なら `101`) | |
| `--ff NAME` | `gaff` | `Unit_Parameter.Comment` に書く力場 ID | |
| `--nh-dof {3N,3N-3}` | `3N-3` | `Q` の自由度の数え方 | |
| `--topology-only` | off | Structure レコードを**書かない** (骨格だけ) | |
| `--initial-gro PATH` | — | `--topology-only` に 1 フレームだけ足す | |
| `--trajectory PATH` | — | 軌跡 (`.xtc` / 多フレーム `.gro`) を全フレーム埋め込む | **要** |
| `--max-frames N` | 全部 | **合計 N 枚**まで間引く | |
| `--frame-step N` | 1 | N 本に 1 本 (ストライド) | |
| `--edr PATH` | — | `.edr` から xvg をその場で作って埋め込む | **要** |
| `--energy PATH` | — | **既にある** `.xvg` を埋め込む | |
| `--skip-nojump` | off | `--trajectory` に `-pbc nojump` を掛けない | |
| `--tpr PATH` | 位置引数の `.gro` | `-pbc nojump` の参照構造 | |
| `--gmx PATH` | `gmx` | 使う gmx の実行ファイル | |
| `--allow-unsupported` | off | 書き出しが正しくない項があっても続行する | |

`--help` で同じ一覧が出ます。

```bash
python -m abmptools.gro2udf --from-top --help
```

### ★ `--trajectory` は gmx を呼ぶ (`-pbc nojump` が既定)

**`--trajectory` を渡すと、変換の前に `gmx trjconv -pbc nojump` が走ります。**
`gro2udf` はもともと gmx を一切呼ばないファイル変換器なので、**ここだけ前提が
変わります**。

そうしてあるのは、**生の `.xtc` をそのまま渡した人が、何も指定せずに正しい
UDF を得られるようにする**ためです。周期境界をまたいだ分子は `.xtc` の中では
割れていて、それが UDF に入ると **OCTA viewer で開いて初めて分かります**。

| 指定 | `-pbc nojump` | gmx |
|---|---|---|
| `--trajectory` のみ (既定) | **走る** | **要る** |
| `--trajectory` + `--skip-nojump` | 走らない | 要らない |
| `--trajectory` 無し | 走らない | **要らない** |

- **`nojump` は冪等**です。既に nojump 済みの軌跡に掛け直しても座標は
  変わりません (PVA 30 分子 2250 原子 × 6 frame で実測、max |Δr| = 0.0000 Å /
  移動した原子 0 個)。ですから `--skip-nojump` は「結果を変えるため」ではなく
  **「gmx を呼ばせないため」**のものです。
  (**ファイルは byte 一致にはなりません。** 同じ実測で 13500 行のうち 5 行が
  `-0.000` と `0.000` の符号だけ違いました。`.gro` の 3 桁表記でのゼロの
  書き方の差で、数値は同じです)
- フラグが `--skip-nojump` (**動作だけ**を述べる) なのは、**nojump 済みでなくても
  意図して飛ばすこと**があるためです —— gmx が無い、別の後処理で通す、
  割れたままの座標を見たい
- 参照構造は `--tpr` があればそれ、無ければ位置引数の `.gro`。`-pbc nojump` は
  結合情報を読まないので `.gro` で成立します。古い gmx が新しい `.tpr` を
  読めないとき (tpx の版違い) は `.gro` に退避して、そう表示します
- `--trajectory` を渡していないのに `--skip-nojump` と書くとエラーになります
  (`--trajectory` の書き忘れを黙って通すと、**topology だけの UDF が「成功」と
  して出ます**)

### `--edr` と `--energy` の違い

**どちらも省けば、エネルギーは一切読みません。** 軌跡だけの UDF になります。

| | 何を渡すか | gmx | いつ使うか |
|---|---|---|---|
| `--edr md.edr` | GROMACS の `.edr` | **要る** | 普通はこちら。`gmx energy` をその場で回して xvg を作る |
| `--energy e.xvg` | `gmx energy` が出した `.xvg` | 要らない | 既に xvg がある / 項を自分で選びたい |

```bash
# --edr: 1 コマンドで済む
python -m abmptools.gro2udf --from-top system.top md/prod.gro \
    --trajectory md/prod.xtc --edr md/prod.edr --out prod.udf

# --energy: gmx energy を自分で呼ぶ形。上と同じ結果になる
gmx energy -f md/prod.edr -o energy.xvg
python -m abmptools.gro2udf --from-top system.top md/prod.gro \
    --trajectory md/prod.xtc --energy energy.xvg --out prod.udf
```

両方渡すとエラーです (どちらがエネルギーの出どころか決められないため)。

> **`--edr` が作る `.xvg` は `.edr` の隣に出ます** (`<edr の stem>_energy.xvg`)。
> 出力先は選べないので、入力ディレクトリを汚したくなければ `.edr` を作業用の
> 場所へ写してから渡してください。

エネルギーがどのフィールドに入るかは
**§ Multi-frame trajectory + energy を 1 UDF に embed** を参照。

### ★ `--max-frames` と `--frame-step` — 間引き

**意味が違います。** 同時には渡せません (エラー)。

| | 意味 | 例 (10000 フレームの `.xtc`) |
|---|---|---|
| `--max-frames 100` | **合計 100 枚**にする | 100 本に 1 本 → 100 枚 |
| `--frame-step 100` | **100 本に 1 本**にする | 100 枚 |

下流で効くのは**枚数**なので、普段は `--max-frames` を使ってください。
どちらも **gmx を必要としません** (読み込み時に落とすだけ)。

**付けないと全フレームが入ります。** UDF はテキストなので素直に効きます ——
実測 (240 原子、同梱サンプル):

| フレーム数 | UDF |
|---|---|
| 0 (`--topology-only`) | 75 KB |
| 1 | 89 KB |
| 11 | 227 KB |

**1 フレームあたり約 13.8 KB = 原子 1 個あたり約 57 バイト**。この比例は
系が大きくなっても変わらず、2250 原子 × 101 フレームの UDF は **12.3 MB**
でした (原子あたり 56.6 バイト)。したがって:

| 系 | フレーム数 | UDF |
|---|---|---|
| 2250 原子 | 101 | **12.3 MB** (実測) |
| 2250 原子 | 10001 (100 ns を 10 ps ごと) | **約 1.3 GB** |
| 2250 原子 | 100 (`--max-frames 100`) | **約 12 MB** |

**目安は 100 フレーム前後**です。それ以上は OCTA viewer 側でも重くなるので、
`--topology-only` で骨格だけ作って軌跡は viewer に読ませる手もあります。

### gmx が PATH に無いとき

`--gmx` で実行ファイルを直接指します。

```bash
python -m abmptools.gro2udf --from-top system.top md/prod.gro \
    --trajectory md/prod.xtc --edr md/prod.edr \
    --gmx /opt/gromacs/bin/gmx --out prod.udf
```

- **MD 環境付属の gmx** (conda 環境、J-OCTA 同梱など) は PATH に出ていない
  ことが多いので、そのときはここで指してください
- **Windows** は J-OCTA 同梱の GROMACS が使えます
  (`C:\J-OCTA-<版>\additional\GROMACS\bin\gmx.exe`)。`-s` に `.gro` を
  渡す経路では `GMXLIB` も要るので、同じ木の `share\top` を指してください ——
  指していないと `residuetypes.dat not found` で止まります

  ```bat
  set "PATH=C:\J-OCTA-12.0\additional\GROMACS\bin;%PATH%"
  set "GMXLIB=C:\J-OCTA-12.0\additional\GROMACS\share\top"
  ```

- gmx が見つからないときのメッセージには、**`--skip-nojump` で止められる**
  ことも併記されます

### テンプレートの決まり方

`--template` を渡さなければ、**同梱の `default_template.udf` を使います**。

```
abmptools/gro2udf/default_template.udf              (cognac11.2 / OCTA85)
abmptools/gro2udf/default_template_cognac101.udf    (--cognac-version 100/101/102)
```

> **`.top` と同名の `.udf` が隣にあっても、自動では使いません。**
> 存在すれば「ある。ただし使っていない」と表示するだけです。そこにあるのは
> たいてい **MD 前の UDF** で、テンプレートは静的構造と箱を供給するため、
> 黙って採ると**MD 前の箱を引き継いだ UDF**ができてしまいます。使いたいなら
> `--template` で明示してください。

### udf-and-gro モード（後方互換）

既にある UDF の**座標とセルだけ**を `.gro` で差し替えます。

```bash
python -m abmptools.gro2udf test.udf output.gro
# -> test_groout.udf  (カレントディレクトリ。名前と場所は選べません)
```

`.top` は要りません。力場もトポロジもテンプレート UDF のものがそのまま残ります。

### Python API

```python
import os.path
import abmptools.gro2udf
from abmptools.gro2udf import Exporter, TopExporter

# udf-and-gro モード
Exporter().export("test.udf", "output.gro")

# --from-top モード
#
# ★ template_path は必須です (None は通りません)。CLI が既定で使うのと
#   同じものを渡すなら、同梱テンプレートのパスを自分で組みます。
template = os.path.join(os.path.dirname(abmptools.gro2udf.__file__),
                        "default_template.udf")

TopExporter().export(
    "system.top", "md/prod.gro",
    template_path=template,
    out_path="prod.udf",
    mdp_path="md/prod.mdp",
    trajectory_path="md/prod_nojump.xtc",
    energy_path="energy.xvg",
    max_frames=100,
)
```

`max_frames` / `frame_step` / `topology_only` / `initial_gro_path` /
`cognac_version` / `force_field` / `nh_dof` / `allow_unsupported` も
同名の引数で受けます。**CLI のオプションとほぼ 1 対 1** です。

**ただし `TopExporter` は gmx を呼びません。** CLI の `--edr` と
`-pbc nojump` は `cli.py` が `abmptools.trajectory` に委譲しているので、
API から同じことをしたいときは自分で呼んでください。

```python
from abmptools.trajectory import gmx_energy, nojump_with_fallback

# -pbc nojump。tpr が古い gmx で読めなければ fallback の .gro に退避する
traj, used_fallback = nojump_with_fallback(
    trajectory="md/prod.xtc",
    reference="md/prod.tpr",
    fallback="md/prod.gro",
)

# .edr -> .xvg (CLI の --edr と同じもの)
xvg = gmx_energy(edr="md/prod.edr", output="energy.xvg", terms=range(1, 51))
```

つまり **API 側は「gmx を呼ぶかどうか」を呼び出し側が決める**構造で、CLI は
その組み合わせを 1 コマンドにまとめたものです。

---

## 動かして確かめる

### サンプルを通す

```bash
cd sample/gro2udf
bash run.sh
```

3 本とも通れば環境は揃っています (UDFManager / MDAnalysis / gmx)。
結果は `run_out/` に出ます —— 詳しくは
[`sample/gro2udf/README.md`](../sample/gro2udf/README.md)。

### 回帰テスト

出力が以前と変わっていないことは pytest が見ます。

```bash
python -m pytest tests/test_regression.py -k gro2udf
```

**UDFManager が無い環境では skip されます** (`importorskip`)。2 件とも通る
ことが「両モードの出力が参照値と一致している」ことの確認です。

> **`sample/gro2udf/*/output/*.udf` を比較対象にしないでください。**
> あれは古い版で作られたスナップショットで、**現在の出力とは一致しません**
> (Nose-Hoover まわりの field が当時のまま)。回帰の参照値はテストが別に
> 持っています。`run.sh` がそれらを上書きしないのもこのためです。

---

## 変換仕様（--from-top モード）

### Nose-Hoover 熱浴質量 Q

UDF フィールド: `Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q`

**単位**: amu·Å²

**計算式**:

```
Q = g · k_B · T · τ²
```

| 変数 | 意味 | 値・出典 |
|---|---|---|
| `g` | 自由度 | `3·N − 3`（N = 全原子数、並進 COM 3 自由度を除く） |
| `k_B` | ボルツマン定数（COGNAC 内部単位） | `0.83144626 amu·Å²/(ps²·K)` |
| `T` | 設定温度 [K] | MDP の `ref_t`（未指定時 300.0 K） |
| `τ` | 緩和時間 [ps] | MDP の `tau_t`（未指定時 0.1 ps） |

**k_B の単位変換**:

```
k_B [amu·Å²/(ps²·K)]
    = 1.38064852×10⁻²³ J/K
      / (1.66054×10⁻²³ J/(amu·Å²/ps²))
    = 0.83144626
```

> MDP を指定しない場合は `T=300.0 K, τ=0.1 ps` のデフォルト値が使用されます。
> 正確な Q を得るには `--mdp` オプションで MDP を渡してください。

---

### Ewald 静電相互作用デフォルト値

`--from-top` モードは以下の Ewald パラメータを UDF に自動設定します。

| UDF パス | 設定値 | 備考 |
|---|---|---|
| `Interactions.Electrostatic_Interaction[0].Name` | `"POINT_CHARGE"` | 固定 |
| `Interactions.Electrostatic_Interaction[0].Algorithm` | `"Ewald"` | 固定 |
| `Interactions.Electrostatic_Interaction[0].Scale_1_4_Pair` | `0.83333333333333` | 5/6（AMBER 慣例） |
| `Interactions.Electrostatic_Interaction[0].Ewald.Dielectric_Constant` | `0.0` | 固定 |
| `Interactions.Electrostatic_Interaction[0].Ewald.R_cutoff` | 下記フォーミュラで計算 **[Å]** | COGNAC 内部単位 |
| `Interactions.Electrostatic_Interaction[0].Ewald.Ewald_Parameters` | `"Auto"` | 固定 |

`Simulation_Conditions.Calc_Potential_Flags` は AMBER/GAFF の慣例に合わせて明示的に設定されます
（テンプレートの既定値は `Angle` / `Torsion` / `Non_Bonding_1_4` が `0` なので、そのままだと
結合角・二面角が書き出されているのに無効、1-4 対も `Scale_1_4_Pair` を設定しているのに落ちる）:

| フラグ | 値 |
|---|---|
| `Bond` / `Angle` / `Torsion` | `1` |
| `Non_Bonding_Interchain` / `Non_Bonding_Intrachain` | `1` |
| `Non_Bonding_1_3` | `0`（除外） |
| `Non_Bonding_1_4` | `1`（`Scale_1_4_Pair` で縮小） |
| `Electrostatic` | `1` |

#### Ewald R_cutoff の計算式

Deserno & Holm (J. Chem. Phys. **109**, 7678, 1998) の最適化条件（実空間誤差 ≈ 逆空間誤差）に基づく式:

```
R_cutoff = sqrt(11.5) / α × 10   [Å]   (COGNAC 内部長さ単位)

α = sqrt(π) × (5.5 N / V²)^(1/6) [nm⁻¹]   (V は nm³)
```

| 変数 | 意味 | 出典 |
|---|---|---|
| N | 全原子数 | TopModel.n_atoms_total |
| V | シミュレーションボックス体積 [nm³] | GRO 第1フレームのセル a×b×c |

GRO ファイルがない場合（空フレーム）は定数フォールバック値 `11.303883305209` Å を使用します。

---

## 変換仕様（udf-and-gro モード）

### フレームフィルタ

`Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps` が N の場合、
GRO ファイル内の `step % N == 0` のフレームのみ UDF レコードとして書き込まれます。

### 単位変換

| 量 | GRO単位 | UDF put 単位 |
|---|---|---|
| 位置 x/y/z | nm | `[nm]` |
| 速度 vx/vy/vz | nm/ps | `[m/s]`（×1000 変換） |
| セルサイズ a/b/c | nm | `[nm]` |
| 時刻 time | ps（`t=` フィールド） | `[ps]` |

### セル形式

- 3値（直交系）: a=box[0], b=box[1], c=box[2], α=β=γ=90°
- 9値（三斜晶）: GROMACS ベクトル表現から a, b, c, α, β, γ を計算

---

## 出力ファイル名

| モード | デフォルト出力ファイル名 |
|---|---|
| udf-and-gro | `<udf_stem>_groout.udf`（カレントディレクトリ） |
| --from-top | `<gro_stem>_fromtop.udf`（カレントディレクトリ） |

## Multi-frame trajectory + energy を 1 UDF に embed

`gen_for_udf.py` で出力される multi-frame `.gro` (or `.xtc`) + `.xvg` を
直接 gro2udf に渡すと、全フレームの座標 + 各フレームのエネルギー値を
1 つの UDF に書き込める。OCTA viewer (GOURMET) で開けば trajectory も energy plot
も同じ UDF から再生される:

```bash
# 1. (前提) MD 実行 + gen_for_udf.py で trajectory + energy を生成
bash run_all.sh
python gen_for_udf.py
# → md/05_npt_final_nojump.gro  (multi-frame、gmx trjconv -pbc nojump)
# → md/05_npt_final_energy.xvg  (gmx energy 全 term)

# 2. gro2udf で全部入りの UDF を生成
#    gen_for_udf.py の出力は既に nojump 済みなので --skip-nojump
#    (付けなくても結果は同じ。付けると gmx を呼ばない)
python -m abmptools.gro2udf --from-top build/system.top md/05_npt_final.gro \
    --mdp md/05_npt_final.mdp \
    --trajectory md/05_npt_final_nojump.gro --skip-nojump \
    --energy md/05_npt_final_energy.xvg \
    --out 05_full.udf

# 3. OCTA viewer (GOURMET) で:
#    File -> Open 05_full.udf
#    → trajectory + energy plot がそのまま見える
```

`--trajectory` 指定時の動作:

- **まず `gmx trjconv -pbc nojump` を通す** (下の「★ `--trajectory` は gmx を
  呼ぶ」)
- `.gro` 拡張子: multi-frame gro を pure-Python parser で読む (各 frame の
  title 行 `t=<ps> step=<N>` から時刻 / step を抽出、coord は fixed-column)
- `.xtc` 拡張子: 既存の MDAnalysis 経由 `frames_from_xtc()` 経路 (要
  `MDAnalysis` install)

> **`--trajectory` は `gmx trjconv -pbc nojump` を先に通します** (既定)。
> 止めるのは `--skip-nojump`、参照は `--tpr`、gmx の場所は `--gmx`。
> 詳しくは **§ 使い方 → `--trajectory` は gmx を呼ぶ**。

`--energy` 指定時の挙動:

- xvg の `@ sN legend "..."` から column 名を取得
- xvg の time grid と trajectory frame time を nearest-neighbour で照合
- 各 frame の `Statistics_Data.Energy.Instantaneous` に以下を書込み:

  | xvg legend | UDF field |
  |---|---|
  | Bond | Bond |
  | Angle | Angle |
  | Proper Dih. + Improper Dih. | Torsion (合算) |
  | LJ-14 + LJ (SR) + Disper. corr. | Nonbonding (合算) |
  | Coulomb-14 + Coulomb (SR) + Coul. recip. | Electrostatic (合算) |
  | Potential | Potential |
  | Kinetic En. | Kinetic |
  | Total Energy | Total |

  schema 上に対応 field が無い場合は silently skip (cognac10.1 で field 数が
  少ない場合)。

- `Hamiltonian` は xvg に直接無いので 0.0 default。

加えて以下も書込まれる (xvg にあれば自動):

| xvg legend | UDF Statistics_Data path | native unit |
|---|---|---|
| Temperature | `Temperature.Instantaneous` | `[T]` = epsilon/R ≈ 503 K |
| Pressure | `Pressure.Instantaneous` | `[P]` = epsilon/(Av·sigma³) |
| Density | `Density.Instantaneous` | `[mass/sigma³]` ≈ 1660 kg/m³ |
| Volume | `Volume.Instantaneous` | `[sigma³]` = 0.001 nm³ |

### Simulation_Conditions.Dynamics_Conditions.Time の MDP 同期

`--mdp` を渡すと、`Simulation_Conditions.Dynamics_Conditions.Time` も MDP から
自動で書き込まれる:

| UDF field | MDP key | 単位 / 意味 |
|---|---|---|
| `delta_T` | `dt` | 1 step の時間 [ps]、UDF schema は [tau] |
| `Total_Steps` | `nsteps` | 全 step 数 |
| `Output_Interval_Steps` | `nstxout-compressed` (fallback: `nstenergy`) | trajectory 出力 step 間隔 |

これで cognac UDF を **後から COGNAC で MD 再開** したい場合の dt / nsteps
情報も保持される。

実機検証 (ketoprofen amorphous、50 mol × 33 atom × 101 frame):
```
totalRecord = 101
Simulation_Conditions.Dynamics_Conditions.Time:
  delta_T              = 0.0205 [tau]  ( = mdp dt=0.001 ps )
  Total_Steps          = 500000
  Output_Interval_Steps = 5000

rec=0:   Bond=381.6, Total=1357.6, T=0.62, P=-0.01, ρ=0.7000, V=18163.0  [native units]
rec=50:  Bond=390.4, Total=1241.8, T=0.57, P=-0.01, ρ=0.6809, V=18672.7
rec=100: Bond=353.6, Total=1361.4, T=0.60, P=-0.00, ρ=0.6694, V=18992.4
```

Time / Cell.a / Density 等は COGNAC native unit (`[tau]` / `[sigma]` /
`[mass/sigma³]`) で書かれるが、OCTA viewer (GOURMET) は `Unit_Parameter` 経由で
[ps] / [nm] / [g/cm³] 表記に戻して描画する。例えば `ρ=0.7 [mass/sigma³]`
は `Unit_Parameter.{Length=0.1, Mass=1}` 経由で `0.7 × (1 amu / (0.1 nm)³) =
0.7 × 1660 kg/m³ ≈ 1162 kg/m³` ≈ 1.16 g/cm³ (ketoprofen の expected density
1.26 g/cm³ 付近)。

## OCTA viewer (GOURMET) 連携: topology-only UDF + 後付け trajectory / energy

OCTA viewer (GOURMET) は **topology だけ含む UDF** を開き、別途 `.gro` (trajectory)
や `.xvg` (energy time-series) を attach するワークフローをサポートする。
gro2udf を `--topology-only` で呼ぶと、Structure record を **0 件** に
した skeleton UDF が出力される (OCTA viewer で trajectory を attach 経由で扱う):

```bash
# (a) topology だけの skeleton UDF (initial frame なし) — Structure record 0 件
python -m abmptools.gro2udf --from-top build/system.top md/05_npt_final.gro \
    --mdp md/05_npt_final.mdp --topology-only --out 05_topology.udf

# (b) 初期 1 frame も含めたい場合は --initial-gro を併用
python -m abmptools.gro2udf --from-top build/system.top build/system.gro \
    --mdp md/05_npt_final.mdp --topology-only \
    --initial-gro md/05_npt_final.gro \
    --out 05_topology_with_initial.udf
```

その後 OCTA viewer で trajectory + energy を attach:

```bash
# trajectory + energy は gen_for_udf.py で生成
python gen_for_udf.py
# → md/05_npt_final_nojump.gro  (multi-frame trajectory)
# → md/05_npt_final_energy.xvg  (energy time-series)

# OCTA viewer (GOURMET) で:
#    File -> Open 05_topology.udf
#    Trajectory -> Load 05_npt_final_nojump.gro
#    Energy plot -> Load 05_npt_final_energy.xvg
```

`--topology-only` 指定時の動作:

- 出力 UDF の `totalRecord` =
  - **0** (`--initial-gro` 省略時、skeleton モード、OCTA viewer で trajectory 全体を attach)
  - **1** (`--initial-gro <path>` 指定時、その gro の 1 frame だけ書く)
- `Set_of_Molecules.molecule[]` には全 mol の atom / bond / angle / torsion
  が書かれる (描画 + topology 認識用)
- `Molecular_Attributes` および `Interactions` も通常通り書かれる
- `Simulation_Conditions.Dynamics_Conditions.Time.{delta_T, Total_Steps,
  Output_Interval_Steps}` も mdp から書かれる

`gen_for_udf.py` の `*_nojump.gro` (`gmx trjconv -pbc nojump`) と組み合わせる
と、PBC を跨いだ連続軌跡が OCTA viewer で滑らかに再生される。OCTA viewer の reference
比較や energy plot 機能を使う際は、本 mode の skeleton UDF を **reference**
として読み込み、別途実トラジェクトリ / xvg を attach する流れ。

## トラブルシューティング

### 多成分系で `Atom_Type_Name` が衝突する (2.13.7 で修正)

OpenFF / interchange の型名 `<moleculetype>_<index>` の index は**分子内**の
通し番号で、分子種ごとに 0 から振り直される。`<元素><index>` に書き換えると
その分子種を区別する情報が落ちるので、**元素と位置が一致する型どうしが衝突する**。

実例 (アリピプラゾール + 乳酸): `APZ_2` と `LAC_2` はどちらも index 2 の O で、
両方 `O2` になっていた。UDF は同名の `Atom_Type` を 2 つ持ち、名前で型を引く経路
(`interaction_Site` / `Pair_Interaction` / 読み戻す全コンバータ) は先頭を拾うので、
**2 成分目が 1 成分目のパラメータで静かに走る**。型 69 個に対し名前 68 種、
往復後の LJ-14 が 0.9 % ずれた。

**単一成分では起きない**ので長く気付かれなかった。

2.13.7 からは topology 全体で一度に名前を割り当てる:

1. 先に現れた型が `<元素><index>` をそのまま取る (**単一成分の出力は従来とバイト一致**)
2. 衝突した型は、その元素で既に使われている最大 index の次を取る
3. 力場固有の型名 (GAFF `c3` 等) は素通しで、先に予約される

衝突を解消したときは警告が出る:

```
atom type 'LAC_2' wanted the display name 'O2', which is already taken by
another moleculetype; using 'O4' instead. The index in an interchange type
name is per molecule, so it repeats across components.
```

**2.13.6 以前が書いた多成分 UDF は作り直すこと。** 名前が潰れている以上、
後から復元することはできない (どちらの型のパラメータが残っているか、
UDF だけからは判別できない)。

### `UDFExportError: failed while writing section ...`

`top_exporter.py` は UDFManager 経由で UDF を書き出す各 stage を
`UDFExportError` で wrap している。OCTA のバージョン違い (例: OCTA84 vs
OCTA85) で template UDF の schema に存在しない field を書こうとした時など、
UDFManager の cryptic な `RuntimeError` をそのまま投げる代わりに、

- **どの section** で失敗したか (`Set_of_Molecules`, `Structure[record=12]`,
  `Molecular_Attributes`, `Interactions`, `default_condition` 等)
- 使った **template UDF / 出力 UDF のパス**
- **underlying exception** (UDFManager の元エラー)
- **対処 hint** (異なる OCTA version の場合は当該 OCTA から template を
  再生成するか、bundled の `abmptools/gro2udf/default_template.udf` を使う)

を含む診断メッセージで再 raise する。エラーが出た時は **section 名 +
underlying** をまず確認し、当該 OCTA version で当該 field が定義されているか
を `<OCTA>/ENGINES/udf/*.udf` (UDF schema 定義ファイル) で確認する。

例 (OCTA84 で `Interactions.Pair_Interaction[].Lennard_Jones.sigma` が無い等):

```
UDFExportError: gro2udf: failed while writing section 'Interactions'.
  template UDF: /opt/OCTA84/ENGINES/cognac/sample/cognac.udf
  output  UDF: ./05_npt_final.udf
  underlying  : RuntimeError: UDFManager: undefined node ...
  hint        : this often means the template UDF schema does not contain a
                field this section needs. If you are using a different OCTA
                version (e.g. OCTA84 vs OCTA85), try regenerating the template
                with that OCTA's `udfreader` / `udfdef.py`, or use the bundled
                template at `abmptools/gro2udf/default_template.udf`.
```

### OCTA8.4 / OCTA8.4 で `file not found:cognac112.udf` エラー

bundled template (`abmptools/gro2udf/default_template.udf`) は OCTA85 の
cognac11.2 schema を前提に `\include{"cognac112.udf"}` を要求する。
**OCTA8.4 / OCTA8.4 には cognac11.2 schema 定義が無く** (cognac101
までしか同梱されていない)、UDFManager が template-open 時点で失敗する:

```
RuntimeError: file not found:cognac112.udf.
[error 1] line 7 near ":": No data definition[Simulation_Conditions].
```

加えて、cognac10.1 → 11.2 で `Structure` セクションの data structure
(`Molecular_Coord` や `Energy_element` 等の field 数) にも変更があり、
**bundled template の data section を OCTA8.4 で parse することはできない**。
よって `--cognac-version 101` で include 行を書き換えても、`Structure` 書込み
時点で別エラーが出る (OCTA85 環境で検証済み)。

### 推奨対処: `--cognac-version 101` で bundled cognac10.1 template を使う

`abmptools` には cognac10.1 schema 互換の bundled template
(`abmptools/gro2udf/default_template_cognac101.udf`) が同梱されています。
`--cognac-version 101` を指定すれば、`--template` 省略時に自動でこの template
が選択されます:

```cmd
python -m abmptools.gro2udf --from-top build\system.top md\05_npt_final.gro ^
    --mdp md\05_npt_final.mdp --cognac-version 101 --out 05_output.udf
```

bundled cognac10.1 template は `Unit_Parameter:{"","",1.0,4.184,0.1}` を含み、
cognac10.1 でも `[nm]` / `[ps]` / `[kJ/mol]` の unit alias が解決されるため、
cognac11.2 を使った場合と **出力データは完全に同一** (Cell.a=26.6805 [sigma]
= 2.66805 nm 等)。

### 代替: OCTA8.4 の GOURMET で保存したミニマル COGNAC UDF を `--template` で渡す

1. **OCTA8.4 の GOURMET を起動**し、OCTA viewer 同梱のミニマル COGNAC sample UDF を
   読み込む。例えば以下のいずれか:
   - `C:\OCTA8.4\bin\win64\ENGINES\cognac\sample\*.udf`
   - `C:\OCTA8.4\GOURMET\sample\cognac*.udf`
   - 既存の自分の COGNAC 計算 input UDF (小さいもの)

2. **File → Save As** で別名保存 (例: `octa84_template.udf`)。GOURMET が
   data section を OCTA8.4 互換に書き出すため、parse 失敗しない。

3. その template を `--template` で渡して gro2udf を実行:

```cmd
python -m abmptools.gro2udf --from-top build\system.top md\05_npt_final.gro ^
    --template octa84_template.udf ^
    --mdp md\05_npt_final.mdp --out 05_output.udf
```

この方法だと OCTA8.4 + cognac101 schema 完全互換の data structure が確保される
ため、Structure / Molecular_Attributes / Interactions 全 section が書ける。

### その他の option (実用上動かないことが多い)

- `--cognac-version 110` / `--cognac-version 101` だけで bundled template の
  include 行を書き換える方法は **data section 構造差で失敗する** ことを確認済み
  (OCTA85 + cognac101 で `Structure[record=0]` 書込み時に
  `RuntimeError: ArgumentError: put data.`)。bundled template が cognac11.2 の
  data 構造を含むため、古い schema では parse できない。
- OCTA8.4 の `udfdef.py` で skeleton UDF を生成する方法もあるが、
  GOURMET の Save As で済むため通常不要。

### `UDFExportError: UDFManager module is required but could not be imported`

OCTA の `python3/` ディレクトリへ `PYTHONPATH` を通していない。section 1.2 の
`OCTA85_HOME` / `UDF_DEF_PATH` / `PYTHONPATH` × 2 (GOURMET + ENGINES)
の 4 行を `~/.bashrc` に追記してから `source ~/.bashrc` で反映する。
