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

> **軌跡だけなら、これを通さなくても構いません。** `gro2udf --trajectory` は
> 渡された軌跡に自分で `-pbc nojump` を掛けるので、生の `.xtc` を直接渡せます
> ([`gro2udf.md`](./gro2udf.md))。`gen_for_udf` の値打ちは **energy.xvg と
> 軌跡を stage 判定込みで揃えて出す**ところにあります。ここを通したものを
> `gro2udf` に渡すときは `--already-nojump` を付けると gmx を呼ばずに済みます
> (付けなくても結果は同じ)。

**stage 名は決め打ちしません。** カレント(または `--dir`)の中で
`<name>.tpr` と `<name>.edr` / `.xtc` / `.trr` が揃っているものを stage として
拾うので、amorphous の `05_npt_final` でも **Tg 計算後の構造でも同じ呼び方**で
通ります。候補が複数あるときは `05_npt_final` / `prod` / `production` を優先し、
それでも決まらなければ候補を列挙して停止します(黙って 1 つ選ばない)。

| オプション | 既定 | 内容 |
|---|---|---|
| `--stage` | 自動判定 | `<stage>.edr` / `.tpr` / `.xtc` の basename |
| `--dir` | カレント | stage ファイルのあるディレクトリ |
| `--ndx` | 使わない | index file。系の一部だけを UDF にするときだけ `--group` とセットで指定 |
| `--ref` | `<stage>.tpr` | `trjconv -s` に渡す構造。古い gmx が tpr を読めないときは `<stage>.gro` へ自動退避 |
| `--nojump-format` | `gro` | `gro` か `xtc`。`xtc` は 10 倍ほど小さいが、読むのに MDAnalysis が要る |
| `--terms-max` | 50 | energy term 番号の上限 |
| `--group` | `0` | trjconv の group(0 = System) |

```bash
# amorphous の md/ で (stage も index も自動)
python -m abmptools.trajectory gen_for_udf

# stage を名指し (例: Tg 計算の出力)
python -m abmptools.trajectory gen_for_udf --stage prod

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

**3. index file** —— 既定では使いません (理由は次節)。

出力は**読んだファイルと同じディレクトリ**に置かれます。
どの stage を選んだかは必ず 1 行目に出るので、実行後にそこを見れば
意図と合っているか確認できます:

```
stage: 05_npt_final
  /path/to/md/05_npt_final_energy.xvg  (gmx energy, 0.3 MB)
  /path/to/md/05_npt_final_nojump.gro  (trjconv -pbc nojump, 4.3 MB)
```

#### index file (`.ndx`) が要るのはどんなときか

**既定では使いません。** `gen_for_udf` が作るのは系全体を写したものなので、
group 0 = tpr の System (全原子) がそのまま欲しいものです。

`--ndx` を渡すと **group 0 の意味が変わります**。

| | group 0 の中身 |
|---|---|
| index なし (既定) | tpr の System = **全原子** |
| index あり | **その index file の最初の group** (System とは限らない) |

つまり無関係な `.ndx` を渡すと、**原子の一部だけを切り出した `.gro` が
エラーなしで**出来ます。`gmx energy` のほうは index を一切使いません。

**`build/system.ndx` は grompp のためのものです。** amorphous の
`run_all.sh` が `grompp -n ../build/system.ndx` として使うのは、mdp の
`tc-grps` が成分ごとの group (`IMC` など) を参照するからで、**軌跡の
切り出し用ではありません**。abmptools が書く `system.ndx` は先頭が必ず
`[ System ]` なので、渡しても渡さなくても結果は同じでした (実測で
`.xvg` / `.gro` とも一致)。効果が無い一方で、先頭 group が System でない
`.ndx` を拾うと黙って壊れます。既定で使わないのはこのためです。

要るのは **「系の一部だけを UDF にしたい」場合だけ**です。その場合は
`--ndx <file> --group <名前か番号>` を明示し、**下流の `.top` も同じ部分系に
揃えてください** —— `gro2udf --from-top` は `.top` と軌跡が同じ系である前提で
組み立てます。

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

  ```
  Fatal error:
  reading tpx file (prod.tpr) version 138 with version 119 program
  ```

  **MD 環境が自前の GROMACS を持っている場合に、これを踏みます。** 手元の
  MD を新しい GROMACS で流していると、付属の古い `gmx` では tpr を読めません。

  `-pbc nojump` は結合情報を使わないので、**reference を `.gro` に替えれば
  通ります**。`gen_for_udf` は上のエラーを検出したとき `<stage>.gro` に自動で
  退避し、警告を出します(`--ref` で明示も可)。

  ただし**同じ結果にはなりません**。nojump は reference から積み上げるので、
  **分子まるごとが別の周期イメージに置かれることがあります**(実測で最大
  23 Å = 約 1 箱の差)。**分子が割れることはありません** — 同じ系で 1 分子の
  最大の広がりは tpr 参照・gro 参照とも 13.58 Å(箱は 22.6 Å)で一致しました。

  なお**数値そのものは gmx の版に依存しません**。同じ `.gro` を reference に
  して GROMACS 2026.3 と 2020.4 で処理した結果は**完全に一致**
  (最大差 0.000e+00 Å)しました。`.edr` も 2020.4 で問題なく読めます。

- **`.gro` を `-s` に渡すときは `GMXLIB` が要ることがある**。設定しないと
  `residuetypes.dat not found` で止まります(`.tpr` はこのデータを自分で
  持っているので出ません)。使っている GROMACS の `share/top` を指します。

  ```cmd
  set "GMXLIB=<gromacs>\share\top"
  ```

  ```bash
  export GMXLIB=<gromacs>/share/top
  ```

## 関連

- 上流の MD 生成: [`amorphous.md`](./amorphous.md)
  (どちらも `wrap_pbc` 相当の後処理を内部で使う)。
- 生成した `.xtc` を可視化・解析へ: [`gro2udf.md`](./gro2udf.md)(COGNAC UDF 化)。
