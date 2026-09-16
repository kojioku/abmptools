# `abmptools.gro2udf` サンプル — GROMACS → COGNAC UDF

```bash
bash run.sh
```

| ディレクトリ | 何を見せるか |
|---|---|
| `udf_and_gro_mode/` | 既存 UDF の座標を `.gro` で差し替える (後方互換モード) |
| `gro_top_mode/` | `.top` + `.gro` から UDF を新規に組む (`--from-top`) |

系はどちらも**ベンゼン 20 分子 / 240 原子**。`gro_top_mode/input/output.gro`
は 11 frame の multi-frame `.gro` で、その 1 frame 目 (箱 1.48 nm、約
0.8 g/cm³ の液体) が構造として使われます。

## エネルギー込みの変換 (`--edr` + `--trajectory`)

`run.sh` の 3 本目です。

```bash
python -m abmptools.gro2udf --from-top \
    gro_top_mode/input/test.top gro_top_mode/input/output.gro \
    --template gro_top_mode/input/template_empty.udf \
    --mdp gro_top_mode/input/md.mdp \
    --trajectory gro_top_mode/input/md.xtc \
    --edr gro_top_mode/input/md.edr \
    --out run_out/md_full.udf
```

**ここだけ gmx が要ります。** `--edr` が `gmx energy` を、`--trajectory` が
`gmx trjconv -pbc nojump` を走らせます (後者は既定。止めるなら
`--skip-nojump`)。`.xtc` を読むにはさらに MDAnalysis が要りますが、
`run.sh` は無ければ `gmx` で `.gro` に直してから渡します。

### ★ energy は 5 frame に 1 回。trajectory は 1 回

`md.mdp` は出力の頻度を**わざと変えて**あります。

| | 値 | 結果 |
|---|---|---|
| `nstenergy` | 100 | **51** rows (0.0 〜 5.0 ps、0.1 ps 刻み) |
| `nstxout-compressed` | 500 | **11** frames (0.0 〜 5.0 ps、0.5 ps 刻み) |

**これが実際の姿**です。amorphous の本番設定は `1000` / `5000` で、やはり
**5:1**。エネルギーは軌跡より細かく取るのが普通なので、1:1 のサンプルでは
**必ず起きる対応付けを一度も確かめられません**。

`gro2udf` は **frame の時刻に最も近い xvg の行**を割り当てます (行番号では
なく時刻で照合)。実データでの対応は次のとおり:

| UDF の record | frame 時刻 | `Instantaneous` | `Batch_Average` | `Total_Average` |
|---|---|---|---|---|
| 0 | 0.0 ps | xvg row 0 | row 0 | row 0 |
| 1 | 0.5 ps | xvg row **5** | row 1〜5 の平均 | row 0〜5 の平均 |
| 2 | 1.0 ps | xvg row **10** | row 6〜10 の平均 | row 0〜10 の平均 |
| … | | | | |
| 10 | 5.0 ps | xvg row **50** | row 46〜50 の平均 | row 0〜50 の平均 |

書き込まれる先は **`Statistics_Data.<Class>.<Avg>.<項目>`** です。
平均の種類が項目より**上**に来ることに注意してください:

```
Statistics_Data.Energy.Instantaneous.Potential
Statistics_Data.Energy.Batch_Average.Potential
Statistics_Data.Energy.Total_Average.Potential
```

Potential の実測値 (この `md.edr`、`gmx energy` の xvg そのまま):

```
frame  t[ps]  Instantaneous  Batch_Average  Total_Average   [kJ/mol]
    0    0.0       1046.452       1046.452       1046.452
    1    0.5        374.611        484.682        578.310
    2    1.0        381.107        289.972        447.247
   ...
   10    5.0        421.669        431.276        488.470
```

最後の `Total_Average` 488.470 は `gmx energy` が報告する全区間平均
(488.47 kJ/mol) と一致します。

> **UDF の中の数字は 4.184 分の 1 になります。** UDF の native 単位が
> kcal/mol だからです。上の 1046.452 は UDF では `250.108`、374.611 は
> `89.534` として入っています (実際に書き出した `run_out/md_full.udf` を
> `UDFManager` で読んで確認した値)。`Time` も同様に native 単位なので、
> 0.5 ps は `10.2` と出ます。

`Temperature` / `Pressure` / `Density` / `Volume` も xvg にあれば同じ規則で
入ります (詳細は [`docs/gro2udf.md`](../../docs/gro2udf.md))。

## `md.edr` / `md.xtc` の作り方

`md.mdp` に全部書いてあります。同じものを作り直すなら:

```bash
# 作業用ディレクトリで (input/ を中間ファイルで汚さないため)
mkdir -p /tmp/gro2udf_md && cd /tmp/gro2udf_md
IN=<abmptools>/sample/gro2udf/gro_top_mode/input
head -243 "$IN/output.gro" > conf.gro     # 1 frame 目だけ取り出す
gmx grompp -f "$IN/md.mdp" -c conf.gro -p "$IN/test.top" -o md.tpr
gmx mdrun -deffnm md -ntmpi 1 -ntomp 1
# md.edr / md.xtc を input/ へ戻す
```

- **カットオフは 0.7 nm** で、`test.mdp` の 0.85 nm とは違います。箱が
  1.48 nm なので半分の 0.74 nm を超えられません (`test.mdp` の 0.85 nm は
  `output.gro` の**最後の** frame、6.17 nm 箱のほうに対応します)
- **バロスタットは使っていません。** 軌跡を作るのが目的で、NVT で足ります
- 速度は `gen_vel = yes` で生成。`gen_seed` を固定してあるので、同じ
  GROMACS なら同じ軌跡が出ます

## 出力

**`run.sh` が書くのは `run_out/` だけ**です (`.gitignore` 済み)。

| 場所 | 追跡 | 内容 |
|---|---|---|
| `run_out/test_groout.udf` | なし | 1 本目の結果 |
| `run_out/output_fromtop.udf` | なし | 2 本目の結果 |
| `run_out/md_full.udf` | なし | 3 本目の結果 (11 frame + energy) |
| `*/output/*.udf` | あり | **参照用のスナップショット。`run.sh` は触りません** |

分けてあるのは、**サンプルを流しただけで作業ツリーが汚れないようにする**ため
です。`--edr` が作る `.xvg` は `.edr` の隣に出て場所を選べないので、`.edr` も
`run_out/` に写してから渡しています。

> `*/output/` のスナップショットは**現在の出力と一致しません**。古い版で
> 作られたもので、Nose-Hoover まわりの field が当時のままです。回帰の比較は
> これらではなく `tests/test_regression.py` が持っている参照値で行います。

## 見えるかもしれない警告

```
Nose-Hoover Q = 422669.4 corresponds to tau_t = 9.659 ps (g = 717).
```

`md.mdp` の `tau_t` は 1.537 ps で、9.659 はその **2π 倍**です。Q と `tau_t`
の換算には 2π を含める流儀と含めない流儀があり、この警告は**下流の
`udf2gro` に渡すときに気付けるように**出しています。変換そのものは正常です
(詳細は [`docs/udf2gro.md`](../../docs/udf2gro.md))。

```
template ... declares a static cell of 20.0000 A, but the data is 14.8011 A.
```

`template_empty.udf` は箱 2.0 nm の空テンプレートで、座標は 1.48 nm の系です。
静的セルは 1 frame 目から取り直されるので、これも正常です。
