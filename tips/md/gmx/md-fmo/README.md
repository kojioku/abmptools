# md-fmo: GROMACS MD → 液滴切り出し → FMO 前処理

タンパク質間相互作用 (PPI) 等の GROMACS MD トラジェクトリから
フレームを切り出し、局所最適化 → 液滴 (水和殻) 抽出 → FMO 入力生成
までを行うスクリプト群。

## 手順

| スクリプト | 役割 |
|---|---|
| `0_parmed.sh` | GROMACS `.top` → Amber `.prmtop` (ParmEd `gromber`)。後段 cpptraj 用 |
| `1_trajsep.sh` | トラジェクトリから指定区間・間隔でフレームを `.gro` に切り出し |
| `2_optmask-frame.sh` | 各フレームを最小化 → PBC 再構成 → 液滴を切り出し (**旧版**) |
| `2_optmask-frame_v2.sh` | 上の PBC 再構成を cpptraj `autoimage` 化した **推奨版** |
| `2b_recover-from-center1.sh` | 既に旧版を流した後の**救済用**。最小化をやり直さず液滴だけ作り直す |
| `check_contact.py` | 切り出した液滴が「複合体が接触した状態」かを機械検証 |
| `3_fmosetup.sh` | 液滴 PDB から ajf (FMO 入力) を生成 |
| `sample/` | 本問題 (複合体の分離・水和殻の穴) を溶媒つき toy 系・数秒で再現するデモ ([sample/README.md](sample/README.md)) |

```bash
bash 0_parmed.sh <system>.top
bash 1_trajsep.sh -f <traj>.xtc -s <system>.gro
cd gmxpdbs-foropt
bash ../2_optmask-frame_v2.sh -n index.solute.ndx -p <top> -f <traj>
bash ../3_fmosetup.sh
```

既に旧 `2_optmask-frame.sh` を流し終えている場合は、`*_fmo_center1.pdb` が残っていれば
最小化をやり直さずに液滴だけ作り直せる (数分):

```bash
bash 2b_recover-from-center1.sh
```

## なぜ v2 (autoimage) を推奨するか

### 症状
複合体 (タンパク A + タンパク B + リガンド等) の系で、
旧 `2_optmask-frame.sh` の PBC 再構成が、
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

## 系ごとの編集ポイント

`2_optmask-frame_v2.sh` 冒頭の以下を系に合わせて変更する:

- `anchormask` / `fixedmask` … 複合体フラグメントの残基レンジ
  (`anchor` = 中心に固定する分子、`fixed` = それに寄せる分子)
- `maskinfo` … 液滴に残す溶質の残基レンジ
- `stripdist` … 溶質から何 Å 以内の水を残すか
- `check_contact.py` の `--frag` … 末尾の検証行も同じレンジに合わせる

**編集し忘れは自動で止まる。** スクリプトは prmtop の `%FLAG POINTERS` から
残基数を読み、上のマスクが範囲内かを起動時に検査する。範囲外なら
`ERROR: anchormask=":1-840" は残基 840 を指しているが, ... は 334 残基しかない`
で中断する。

> **なぜ検査が要るか。** cpptraj は範囲を超えたレンジを**無警告で丸める**。
> 334 残基の系に `:1-840` と書くと全 1055 原子 (= 溶媒込み) が選ばれ、警告は
> 一切出ない。anchor に水が全部入った状態で autoimage が回り、静かに間違った
> 液滴が出来る。逆に開始側が範囲外の `:841-1185` は警告こそ出るが選択ゼロで
> 走り続ける。どちらも実行時には気づけない。

## スレッド数とログインノード

`minimize` の `gmx mdrun` は `-ntmpi 1 -ntomp <n>` で回し、`OMP_NUM_THREADS`
も同じ値にそろえる。既定は 4 で、`-t` か環境変数 `OPT_THREADS` で変えられる。

```bash
bash 2_optmask-frame_v2.sh -n index.solute.ndx -p system.top -f traj.xtc -t 8
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
- **同じディレクトリで再実行できる。** `${head}_ref.tpr` と minimize 済みの
  フレームは `[skip]` で飛ばすので、途中で失敗しても直して流し直せばよい。

## 検証 (必須)

液滴を FMO に投入する前に必ず通し、NG フレームは使わないこと。

```bash
python3 check_contact.py \
    --frag A:1-840 --frag B:841-1184 --frag lig:1185 \
    <head>-optedpdb/*_fmo_mask.pdb
```

- フラグメント間の実座標最小距離が `--contact` (既定 5 Å) 未満 → 接触 OK
- 水和殻の最大空隙が `--hole` (既定 8 Å) 未満 → イメージング健全
- 依存: `numpy`, `scipy`。システムの `python3` に無い環境では
  `PYTHON=/path/to/venv/bin/python` で差し替える (富岳のログインノードは
  `python3` が 3.6.8 で scipy 無し)
- NG が 1 つでもあれば exit 1 を返すのでバッチに組み込める

**接触距離だけを見ないこと。** `-pbc cluster` は「接触は戻すのに水和殻を置き去りにする
(溶質が脱水する)」失敗をする。溶媒つき toy 系でも再現でき (空隙 3.7 Å → 13.9 Å、
接触距離は正常のまま)、実データの破綻フレームでは空隙が 19–28 Å に達した。
`max_gap` を必ず併せて確認する。

## NG が出たときの切り分け

1. **抽出直後の生 `.gro` を先に調べる**。
   ```bash
   python3 check_contact.py --frag ... gmxpdbs-foropt/<head>_<i>.gro
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
