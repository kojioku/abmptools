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
| `check_contact.py` | 切り出した液滴が「複合体が接触した状態」かを機械検証 |
| `3_fmosetup.sh` | 液滴 PDB から ajf (FMO 入力) を生成 |
| `sample/` | 本問題を数十 KB・数秒で再現する極小 toy 系デモ ([sample/README.md](sample/README.md)) |

```bash
bash 0_parmed.sh CR8_AF3_MD_cpx_solv_noIon.top
bash 1_trajsep.sh -f CR8_AF3_MD_pr.xtc -s CR8_AF3_MD_pr.gro
cd gmxpdbs-foropt
bash ../2_optmask-frame_v2.sh -n index.solute.ndx -p <top> -f <traj>
bash ../3_fmosetup.sh
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

## 検証 (必須)

液滴を FMO に投入する前に必ず通し、NG フレームは使わないこと。

```bash
python3 check_contact.py \
    --frag A:1-840 --frag B:841-1184 --frag lig:1185 \
    <head>-optedpdb/*_fmo_mask.pdb
```

- フラグメント間の実座標最小距離が `--contact` (既定 5 Å) 未満 → 接触 OK
- 水和殻の最大空隙が `--hole` (既定 8 Å) 未満 → イメージング健全
- 依存: `numpy`, `scipy`
