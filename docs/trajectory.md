# abmptools.trajectory — GROMACS trajectory 後処理ラッパー

`gmx trjconv` / `gmx energy` を **クロスプラットフォーム(Windows native 対応)の
Python API + CLI** で叩く後処理ユーティリティ。従来 sample / builder に散在していた
bash script(`trajectory_thin_nojump.sh`、amorphous の `wrap_pbc.sh`、`gen_for_udf.sh`
等)を Python に統一したもの。GROMACS 本体は同梱せず `subprocess` で呼ぶだけ。

## CLI

```
python -m abmptools.trajectory <subcommand> --traj <xtc> --tpr <tpr> [options]
```

| subcommand | 内容 | 対応する gmx |
|---|---|---|
| `thin_nojump` | `--skip N` で間引き + `-pbc nojump`(分子を割らない連続座標) | `trjconv -skip N -pbc nojump` |
| `nojump` | `-pbc nojump` のみ | `trjconv -pbc nojump` |
| `thin` | `--skip N` で間引きのみ | `trjconv -skip N` |
| `wrap_pbc` | `-pbc mol -ur compact`(+ 任意 `--center`)。VMD 向け compact 表示 | `trjconv -pbc mol -ur compact` |
| `energy` | `.edr` → `.xvg`(エネルギー項抽出) | `energy` |
| `gen_for_udf` | OCTA / gro2udf 用に `energy` + `nojump` をまとめて実行。stage は自動判定 | `energy` + `trjconv -pbc nojump` |

### 共通オプション(座標系サブコマンド)

| オプション | 既定 | 内容 |
|---|---|---|
| `--traj` | (必須) | 入力 trajectory(`.xtc` 等) |
| `--tpr` | (必須) | reference 構造(`.tpr` 推奨) |
| `--out` | 自動命名 | 出力パス(省略時は `<traj>_<op>_skip<N>.xtc` 等) |
| `--group` | `System` | trjconv の出力グループ |
| `--ndx` | なし | index file |
| `--gmx` | `gmx` | 使用する gmx 実行ファイル(下記 gotcha 参照) |

`energy` サブコマンドは `--edr` / `--out`(`.xvg`)/ `--terms-max`(既定 50)/ `--gmx`。

### 例

```bash
# 間引き + nojump(prod.xtc を 1/10 に、分子を割らずに連続化)
python -m abmptools.trajectory thin_nojump --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10

# VMD 向け compact wrap
python -m abmptools.trajectory wrap_pbc --traj 05_npt_final.xtc --tpr 05_npt_final.tpr \
    --out 05_npt_final_pbc.xtc --ur compact

# エネルギー項を xvg に
python -m abmptools.trajectory energy --edr prod/prod.edr --out prod_energy.xvg

# OCTA / gro2udf に渡す 2 点セット (energy.xvg + nojump.gro) をまとめて
python -m abmptools.trajectory gen_for_udf
```

### `gen_for_udf` — OCTA へ渡す 2 点セット

`<stage>_energy.xvg`(全エネルギー項)と `<stage>_nojump.gro`(PBC を跨いで
連続な軌跡)を 1 コマンドで作ります。amorphous の `md/gen_for_udf.py` の中身は
これです。

**stage 名は決め打ちしません。** カレント(または `--dir`)の中で
`<name>.tpr` と `<name>.edr` / `.xtc` / `.trr` が揃っているものを stage として
拾うので、amorphous の `05_npt_final` でも **Tg 計算後の構造でも同じ呼び方**で
通ります。候補が複数あるときは `05_npt_final` / `prod` / `production` を優先し、
それでも決まらなければ候補を列挙して停止します(黙って 1 つ選ばない)。

| オプション | 既定 | 内容 |
|---|---|---|
| `--stage` | 自動判定 | `<stage>.edr` / `.tpr` / `.xtc` の basename |
| `--dir` | カレント | stage ファイルのあるディレクトリ |
| `--ndx` | 自動探索 | index file。既定は `../build/system.ndx` → `system.ndx` の順に探す |
| `--no-ndx` | — | index を使わない(自動探索も止める) |
| `--terms-max` | 50 | energy term 番号の上限 |
| `--group` | `0` | trjconv の group(0 = System) |

```bash
# amorphous の md/ で (stage も index も自動)
python -m abmptools.trajectory gen_for_udf

# Tg 計算の出力 (index file が無い run)
python -m abmptools.trajectory gen_for_udf --stage prod --no-ndx

# 別ディレクトリを指定して
python -m abmptools.trajectory gen_for_udf --dir run1/md
```

#### 引数なしのとき、何が読まれるか

`gen_for_udf` は必須引数がありません。引数なしで実行すると、カレント
ディレクトリ(`--dir` があればそちら)の中身から読むファイルを決めます。

**1. どの stage か** —— `*.tpr` を列挙し、**同じ basename の `.edr` /
`.xtc` / `.trr` が 1 つでもある**ものだけを候補にします(`.tpr` だけの
grompp 残骸は無視)。

- 候補が 1 つ → それを使う(名前は何でもよい)
- 複数 → `05_npt_final` → `prod` → `production` の順で優先
- それでも決まらない → **候補を列挙して RC 1 で停止**

**2. 決まった stage から読むファイル**

| 読むもの | 用途 | 無い場合 |
|---|---|---|
| `<stage>.edr` | `gmx energy` → `<stage>_energy.xvg` | energy をスキップ (明示表示) |
| `<stage>.trr` or `.xtc` | `gmx trjconv -pbc nojump` の入力 | 軌跡をスキップ (明示表示) |
| `<stage>.tpr` | 上の reference 構造 | 軌跡をスキップ |

`.trr` と `.xtc` が両方あれば、**速度も持つ `.trr` を優先**します。
両方ともスキップになれば RC 1 です(何も作らずに「完了」と言わないため)。

**3. index file** —— `../build/system.ndx` → `system.ndx` の順に探し、
見つかれば使い、無ければ使いません。

出力は**読んだファイルと同じディレクトリ**に置かれます。
どの stage を選んだかは必ず 1 行目に出るので、実行後にそこを見れば
意図と合っているか確認できます:

```
stage: 05_npt_final
index: /path/to/build/system.ndx
  /path/to/md/05_npt_final_energy.xvg  (gmx energy, 0.3 MB)
  /path/to/md/05_npt_final_nojump.gro  (trjconv -pbc nojump, 4.3 MB)
```

#### `--no-ndx` を使うとき

index file (`.ndx`) は原子を group にまとめた定義で、`gmx trjconv` の `-n`
に渡すものです。`--no-ndx` は**自動探索ごと止めて `-n` を付けずに実行**します。

なぜ要るかというと、**group 番号の意味が `.ndx` の有無で変わる**からです。

| | group 0 の中身 |
|---|---|
| `.ndx` なし | tpr の既定 group。**0 = System** |
| `.ndx` あり | **その `.ndx` の最初の group** (System とは限らない) |

`gen_for_udf` は group 0 を出力するので、隣に**別の系・別の目的で作った
`.ndx`** が置いてあると、意図せず一部の原子だけを切り出した `.gro` が
**エラーなしで**出来てしまいます。Tg 計算のように index を使わない run の
ディレクトリで、たまたま `system.ndx` が同居しているようなときに付けて
ください。

逆に「別の `.ndx` の、特定の group を出したい」場合は
`--ndx <file> --group <名前か番号>` を明示します。index が元から無い run
では、付けても付けなくても結果は同じです。

## Python API

```python
from abmptools.trajectory import (
    thin_and_nojump, nojump, thin, wrap_pbc, gmx_energy, gen_for_udf,
    find_stage, find_ndx, run_trjconv, GmxError,
)

out = thin_and_nojump(trajectory="prod/prod.xtc", tpr="prod/prod.tpr", skip=10)
# -> PosixPath('prod/prod_nojump_skip10.xtc')
```

| 関数 | 内容 |
|---|---|
| `thin_and_nojump(trajectory, tpr, skip=..., ...)` | 間引き + nojump。出力 Path を返す |
| `nojump(...)` / `thin(...)` | それぞれ単独 |
| `wrap_pbc(..., ur="compact", center=None)` | `-pbc mol -ur compact` |
| `gmx_energy(edr, out, terms_max=50, gmx="gmx")` | `.edr` → `.xvg` |
| `gen_for_udf(stage=None, directory=".", ndx=None, ...)` | energy + nojump をまとめて実行。`{"stage", "energy", "trajectory", "ndx"}` の dict を返す |
| `find_stage(directory)` / `find_ndx(directory)` | stage 名 / index file の自動判定(単体でも使える) |
| `run_trjconv(...)` | 低レベル `trjconv` ラッパー(任意フラグ) |
| `GmxError` | gmx 実行失敗時に送出される例外 |

## 注意 (gotcha)

- **tpr の版と gmx の版が合わないと失敗する**。新しい GROMACS で書いた `.tpr`
  (例: 2026 系の tpr v138)は古い `gmx` では読めない。使いたい `gmx` を
  **`--gmx /path/to/gmx` で明示**して、tpr を書いた版と揃える。
- グループ選択は `--group`(既定 `System`)。溶質だけ等にしたい場合は `--ndx` +
  グループ名を指定。

## 関連

- 上流の MD 生成: [`amorphous.md`](./amorphous.md)
  (どちらも `wrap_pbc` 相当の後処理を内部で使う)。
- 生成した `.xtc` を可視化・解析へ: [`gro2udf.md`](./gro2udf.md)(COGNAC UDF 化)。
