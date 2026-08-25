# md-fmo: GROMACS MD → 液滴切り出し → FMO 前処理

タンパク質間相互作用 (PPI) 等の GROMACS MD トラジェクトリから
フレームを切り出し、局所最適化 → 液滴 (水和殻) 抽出 → FMO 入力生成
までを行うスクリプト群。

## 手順

| スクリプト | 役割 |
|---|---|
| `0_parmed.sh` | GROMACS `.top` → Amber `.prmtop` (ParmEd `gromber`)。後段 cpptraj 用 |
| `1_trajsep.sh` | トラジェクトリから指定区間・間隔でフレームを `.gro` に切り出し |
| `2_optmask-frame.sh` | 各フレームを最小化 → PBC 再構成 (cpptraj `autoimage`) → 液滴を切り出し |
| `2b_recover-from-center1.sh` | 既に旧版を流した後の**救済用**。最小化をやり直さず液滴だけ作り直す |
| `check_contact.py` | 切り出した液滴が「複合体が接触した状態」かを機械検証 |
| `3_fmosetup.sh` | 液滴 PDB から ajf (FMO 入力) を生成 |
| `sample/` | 本問題 (複合体の分離・水和殻の穴) を溶媒つき toy 系・数秒で再現するデモ ([sample/README.md](sample/README.md)) |

```bash
bash 0_parmed.sh <system>.top
bash 1_trajsep.sh -f <traj>.xtc -s <system>.gro
cd gmxpdbs-foropt
bash ../2_optmask-frame.sh -n index.solute.ndx -p <top> -f <traj>
bash ../3_fmosetup.sh
```

既に旧版 (`-pbc cluster` を使っていた頃の `2_optmask-frame.sh`) を流し終えている場合は、`*_fmo_center1.pdb` が残っていれば
最小化をやり直さずに液滴だけ作り直せる (数分):

```bash
bash 2b_recover-from-center1.sh
```

## なぜ cpptraj `autoimage` を使うか

### 症状
複合体 (タンパク A + タンパク B + リガンド等) の系で、
旧版の PBC 再構成 (`-pbc cluster` を挟む 4 段) が、
**もともと正しく接触していた複合体を引き剥がす**ことがある。
壊れたフレームでは複合体が数十 Å 離れ、水和殻に穴が空く。
見た目では気づきにくく、そのまま FMO に投入すると物理的に無意味な結果になる。

200 ns 系で 152.5–200 ns から 2.5 ns 間隔・20 フレームを切り出した実測:

| 段階 | A–B 接触が壊れたフレーム数 | 液滴の水分子数 | 水和殻の穴 |
|---|---|---|---|
| 抽出直後の生 `.gro` | **0 / 20 (全て正常)** | — | — |
| 旧 `-pbc cluster` | **11 / 20** | 5865–7442 (欠損あり) | 最大 28 Å の空洞 |
| v2 `autoimage` | **0 / 20** | 7149–7441 (安定) | 0 |

### 原因
1. Amber で構築した複合体を acpype で GROMACS に変換すると、
   溶質 (タンパク A+B+リガンド) が **1 つの `[moleculetype]`** に潰れる
   (Amber prmtop に moleculetype 再利用の概念が無いため)。
2. その結果、GROMACS の安全な PBC 手段が使えない:
   - `-pbc whole` … 結合をたどるだけでフラグメント間の相対位置は直せない
   - `-pbc mol` … 溶質全体を 1 分子扱いするので直せない
   - `-pbc cluster` … 唯一相対位置を動かせるが、反復アルゴリズムが
     **初期配置依存**で、正しい複合体を壊すギャンブルになる
3. index (`.ndx`) は原子番号の集合にすぎず、分子のつながり情報を持たない。
   index を細かく分けても `-pbc` は複数グループを剛体として扱えないので解決しない。

### 解決
cpptraj `autoimage` は **prmtop の結合情報から各フラグメントを別分子と認識**し、
「anchor に対する最近接イメージ」を反復なしで決定論的に選ぶため破綻しない。
どうせ後段の液滴切り出しで prmtop を使う (`0_parmed.sh`) ので、
PBC 再構成もそちらに寄せると一貫する。

```
# 旧: -pbc whole -> -pbc cluster -> -pbc mol -ur compact -center -> -fit rot+trans
# 新: -pbc whole -> cpptraj autoimage anchor :A fixed :B+lig mobile :WAT -> mask
```

## 系ごとの指定

系ごとの値は **スクリプト冒頭を編集する**。実行はそのまま:

```bash
bash ../2_optmask-frame.sh -n index.solute.ndx -p system.top -f traj.xtc
```

| 冒頭の変数 | 意味 |
|---|---|
| `anchormask` | 箱の中心に固定する分子 (通常タンパク片方) |
| `fixedmask` | anchor の最近接イメージへ寄せる分子 (相手タンパク + リガンド) |
| `mobilemask` | 箱に詰め直す溶媒 (既定 `:WAT`) |
| `maskinfo` | 液滴に残す溶質の残基レンジ |
| `stripdist` | 溶質から何 Å 以内の水を残すか (既定 6.0) |
| `pdb_snum` / `pdb_enum` / `pdb_interval` | 処理するフレーム番号 |
| `minscript` | minimize の条件ファイル |

**編集するのは冒頭のブロックだけ。** 以前は末尾の `--frag` 行も直す必要が
あり、480 行離れた 2 か所を揃えなければならなかった。今は下記のとおり
自動導出する。

**コマンドラインで渡した値は、冒頭の値より必ず優先される。** 引数の解釈は
冒頭のブロックより後で行われるので、一時的に別の系・別のフレームで流したい
ときはスクリプトを書き換えずに済む。解決後の値は実行時に表示される:

```
--- 設定 ---
  anchor : :2            <- --anchor :2 が冒頭の :1 を上書きした
  fixed  : :3
  solute : 1-2  (液滴に残す水: 溶質から 9.9 A)
  frames : 1..3 step 2
  mdp    : min_quick.mdp
  verify : --residues anchor:2-2 --residues fixed:3-3
  threads: 4
```

対応は以下のとおり:

| オプション | 対応する変数 |
|---|---|
| `--anchor` `--fixed` `--mobile` | `anchormask` / `fixedmask` / `mobilemask` |
| `--solute` `--strip` | `maskinfo` / `stripdist` |
| `--frames s:e[:i]` | `pdb_snum` / `pdb_enum` / `pdb_interval` |
| `--mdp` | `minscript` |
| `--fit` | `fitmask` (指定すると `dofit=1`) |

実行時に解決後の設定を表示するので、編集漏れはその場で分かる。

> **検証用の残基レンジは書かなくてよい。** `check_contact.py` に渡す
> `--residues` は `--anchor` / `--fixed` から自動で作る
> (`anchor:1-323` / `fixed:324-324`)。同じレンジを二度書く必要はない。
> 以前は冒頭のマスクと末尾の `--frag` が 480 行離れて別々に置かれており、
> 片方だけ直して気づかない形になっていた。マスクが単純なレンジでない場合だけ
> `--residues <名前>:<開始>-<終了>` を明示する。

上のオプションを省いたときの既定値はスクリプト冒頭にある:

- `anchormask` / `fixedmask` … 複合体フラグメントの残基レンジ
  (`anchor` = 中心に固定する分子、`fixed` = それに寄せる分子)
- `maskinfo` … 液滴に残す溶質の残基レンジ
- `stripdist` … 溶質から何 Å 以内の水を残すか

**編集し忘れは自動で止まる。** 起動時に 2 段の検査をする。

1. **範囲**: prmtop の `%FLAG POINTERS` から残基数を読み、マスクが超えていないか
2. **溶媒混入**: `%FLAG RESIDUE_LABEL` を読み、`anchormask` / `fixedmask` に
   `WAT` / `HOH` / `SOL` / `TIP3` / `T3P` が入っていないか

> **2 段目が要る理由。** 溶媒つきの系では NRES が水込みで数万になるので、
> 溶質しか無い前提で書いたマスクが**範囲内に収まって素通りする**。実測:
> EGFR (NRES=22215) に同梱の例マスク `:1-840` をかけると 1 段目は通り、
> 水 516 残基が anchor に入る。autoimage の基準が壊れるが、エラーにも
> 警告にもならない。2 段目はこれを
> `ERROR: anchormask=":1-840" に溶媒が 516 残基入っている` で止める。

> **なぜ検査が要るか。** cpptraj は範囲を超えたレンジを**無警告で丸める**。
> 334 残基の系に `:1-840` と書くと全 1055 原子 (= 溶媒込み) が選ばれ、警告は
> 一切出ない。anchor に水が全部入った状態で autoimage が回り、静かに間違った
> 液滴が出来る。逆に開始側が範囲外の `:841-1185` は警告こそ出るが選択ゼロで
> 走り続ける。どちらも実行時には気づけない。

## ジョブ投入スクリプト (3_fmosetup.sh)

`3_fmosetup.sh` は最後に投入スクリプト (`setupv1dd2024-bulk.sh` 等) を使う。
**手でコピーする必要はない。** 見つからなければ

1. カレント
2. `3_fmosetup.sh` と同じ場所
3. 上へ 3 階層

の順に探し、見つけたらカレントへ複製する。名前が冒頭の `runsh=` と違っていても、
候補が 1 つに決まるならそれを使い、読み替えたことを表示する。

```
[copy] ../setupv1dd2024-bulk.sh を使う (runsh を setupv1dd2024-bulk.sh に読み替えた)
runsh = setupv1dd2024-bulk.sh
```

どこにも無い場合は、探した場所を示して止まる。

> **以前は「カレントに無ければ即終了」だった。** 実際に `3_fmosetup.sh` を流すのは
> `-optedpdb` の中で、投入スクリプトは配布物と一緒に上の階層にある。手で `cp` する
> 手順が 1 つ増えるだけで、忘れると理由の分かりにくいエラーになっていた。

## 通し確認を短時間で済ませる

本番の `min_aftermd.mdp` は `emtol=1255` まで落とすので、70950 原子の系だと
1 フレーム 5-7 分かかる (収束 step は構造次第で 294-711 と 2.4 倍ばらつく)。
段の繋がり・再開・検証だけを確かめたいときは `min_quick.mdp` を使う。

```bash
bash ../2_optmask-frame.sh -n index.solute.ndx -p system.top -f traj.xtc \
    --mdp min_quick.mdp
```

`nsteps=15` で打ち切り、カットオフも 1.2 nm に縮めてある。同じ EGFR 系で
**5 フレーム 36 秒**。

> **min_quick.mdp で作った液滴を FMO に投入しないこと。** 力が全く落ちて
> いないので構造として使えない。パイプラインの配線を見るためだけのもの。

## スレッド数とログインノード

`minimize` の `gmx mdrun` は `-ntmpi 1 -ntomp <n>` で回し、`OMP_NUM_THREADS`
も同じ値にそろえる。既定は 4 で、`-t` か環境変数 `OPT_THREADS` で変えられる。

```bash
bash 2_optmask-frame.sh -n index.solute.ndx -p system.top -f traj.xtc -t 8
```

- **`-nt` だけを渡さないこと。** 環境に `OMP_NUM_THREADS` が居ると
  `The total number of threads requested (12) is not divisible by the number
  of OpenMP threads requested (8)` で mdrun が即死する。HPC のジョブ
  スクリプトはたいてい `OMP_NUM_THREADS` を撒くので、環境任せにしない。
- **共有のログインノードで回すときは既定のまま (4) にしておく。** コアを
  埋めると他の利用者の作業が止まる。大きくするのは計算ノードにジョブとして
  流すときだけにする。

## 再実行と中断

- **途中で失敗したらそこで止まる** (`set -euo pipefail` + 段ごとの出力検査)。
  以前は `gmx mdrun` が落ちても走り続け、真因の数十行下で cpptraj の空実行 →
  `mv: cannot stat` → `check_contact` の `FileNotFoundError` と別のエラーが
  4 つ出て、原因が追いにくかった。
- **同じディレクトリで再実行できる。** 済んだ段は `[skip]` で飛ばすので、
  途中で失敗しても直して流し直せばよい。

### 「済んだ」の判断は完了マーカーだけで行う

再開でいちばん危ないのは、**完成していないものを完成とみなすこと**である。
判断の根拠をファイルの存在に置くと、次の 2 つを取り違える。

- `mkdir` は作業の**前**に走る。**ディレクトリがあっても中身は無い**
- 殺されて切り詰められた `.gro` / `.pdb` も「存在する」し「空でもない」。
  `[ -s ファイル ]` は通ってしまう

そこで、フレームごとに段の完了マーカー (`.done.minimize` / `.done.arrange`) を
**中身を検めたあと、最後に**書き、判断はこれだけで行う。

- `.gro` は 2 行目の原子数と実際の行数 (N+3) が一致し、最終行が箱ベクトルで
  あることまで見る。書きかけはここで落ちる
- `.pdb` は原子行があり `END` で終わっている (または最終行が座標欄まで
  書けた完全な原子行である) ことを見る
- 液滴は **一時名に書かせてから rename** する。同一ファイルシステム内の
  rename は原子的なので、殺されても半端な `_mask.pdb` が残らない
- マーカーが無いのに成果物がある場合 (前回ここで落ちた、あるいはマーカー
  導入前の実行) は、**完全なら採用してマーカーだけ付け、壊れていれば
  park してから作り直す**

```bash
# 何がどこまで済んでいるか (何も実行しない。prmtop が無くても見られる)
bash 2_optmask-frame.sh -n index.ndx -p system.top -f traj.xtc --check

frame 0 | minimize: 済 | arrange: 済
frame 1 | minimize: 済 | arrange: 途中(成果物あり・マーカー無し)
frame 2 | minimize: 未 | arrange: 未
frame 3 | minimize: 途中(成果物が壊れている) | arrange: 未(ディレクトリのみ)

# やり直したいとき
bash 2_optmask-frame.sh ... --redo              # 全部
bash 2_optmask-frame.sh ... --redo-from arrange # arrange から
```

### 上書きで消さない

作り直しが要る場面では、既存を `mv -f` で潰さず `<名前>.<日時>` にどける
(`[park]` と表示される)。`-optedpdb/` へのコピーも同じ扱いなので、**再実行で
前回の液滴が黙って消えることはない**。`_whole.pdb` / `_arranged.pdb` は同じ段で
作り直す中間物で、前者は gmx が `#file#` として自動退避する。

## 検証 (必須)

液滴を FMO に投入する前に必ず通し、NG フレームは使わないこと。

```bash
python3 check_contact.py \
    --residues A:1-840 --residues B:841-1184 --residues lig:1185 \
    <head>-optedpdb/*_fmo_mask.pdb
```

- 指定した部分どうしの実座標最小距離が `--contact` (既定 5 Å) 未満 → 接触 OK
- 水和殻の最大空隙が `--hole` (既定 8 Å) 未満 → イメージング健全
- 依存: `numpy`, `scipy`。**python は自動で探す** —— 有効な venv、`~/fmoenv`、
  `python3`、`python` の順に試し、`numpy` と `scipy` を import できたものを使う。
  明示したいときだけ `PYTHON=/path/to/python` を渡す
- 見つからない場合は**液滴を残したまま検証だけ保留**し、後から流すコマンドを
  出して exit 2 で終わる。「検証したつもり」にはしない。富岳のシステム
  `python3` は 3.6.8 で scipy が無いので、python 環境を作る前に流すとこうなる
- `amber.sh` が `PYTHONPATH` に入れる Amber の py3.9 site-packages は
  **スクリプト側で外す**。残したまま venv の python から `parmed` を import すると
  numpy 2 系で `No module named 'numpy.compat'` になる。利用者が自分で設定した
  `PYTHONPATH` は残す
- NG が 1 つでもあれば exit 1 を返すのでバッチに組み込める
- **測れなかったものは合格にしない。** 指定した残基レンジに該当が無い、
  比較できた組が 1 つも無い、液滴なのに水が 0 — いずれも NG にして理由を出す
  (以前は素通りしていた。原子 50 個・水 0 に切り詰めた液滴が
  `water= 0 max_gap= nan OK` と報告された)。液滴化前の構造を検証する等で
  水が無くて当然の場合は `--allow-dry`

> **`--residues` に渡すのは残基番号であって、FMO のフラグメント番号ではない。**
> FMO の分割は次の `3_fmosetup.sh` (ajf 生成) で決まるので、この時点では
> まだ存在しない。溶質の範囲では一致することが多いが保証は無く、実データ
> (EGFR) の `AUTOMATIC FRAGMENTATION` 表では残基通番 335 が FMO
> フラグメント 325 に対応する (液滴に残らなかった水の分だけずれる)。
> cpptraj の `:N-M` と同じ番号を渡すこと。
> 旧名 `--frag` も当面受け付けるが、この取り違えを招くので `--residues` を使う。

**接触距離だけを見ないこと。** `-pbc cluster` は「接触は戻すのに水和殻を置き去りにする
(溶質が脱水する)」失敗をする。溶媒つき toy 系でも再現でき (空隙 3.7 Å → 13.9 Å、
接触距離は正常のまま)、実データの破綻フレームでは空隙が 19–28 Å に達した。
`max_gap` を必ず併せて確認する。

## NG が出たときの切り分け

1. **抽出直後の生 `.gro` を先に調べる**。
   ```bash
   python3 check_contact.py --residues ... gmxpdbs-foropt/<head>_<i>.gro
   ```
   ここで既に NG なら **PBC の問題ではなく、実際に解離している**可能性が高い
   (長時間 MD では起こりうる)。生が正常なら壊しているのは後処理側。
   実測例では生 20/20 が正常で、破綻は全て後処理由来だった。
2. **`anchor` を変えてみる** — 相手タンパクを anchor にすると通ることがある。
3. NG フレームは**除外する**。無理に直すより使わない方が安全。

## よくある落とし穴 (道具側)

1. **conda/miniconda を PATH 先頭に置くと `cpptraj` が AmberTools 版でなくなる**。
   同梱の別ビルドを掴む。python はフルパス変数で持つのが安全:
   ```bash
   source <gromacs>/bin/GMXRC
   source <amber>/amber.sh
   PY=<miniconda>/bin/python3      # PATH には入れない
   ```
2. **`gmx trjconv -pbc cluster` に `.pdb` を渡すと無言で空ファイルを吐く**。
   中間ファイルは `.gro` で繋ぐこと。`-pbc whole` は `.pdb` でも通るので気づきにくい。
3. **`gmx make_ndx` は `keep 0` を先に打たないとグループ番号がずれる**。
   ずれたまま cluster にかけると
   `Fatal error: Molecule 1 marked for clustering but not atom 1 in it` で落ちる。
   ```
   keep 0
   a 1-<溶質末尾の原子番号>
   name 1 solute
   q
   ```
4. **cpptraj `autoimage` の `mobile` マスクは対象が空だとセットアップで落ちる**
   (`Could not set up image list` → `Could not allocate mobile list`)。
   無溶媒系では `mobile` を**省く**。
5. `mask ... maskpdb` の出力は連番 `.1` が付くので `mv` が要る。
