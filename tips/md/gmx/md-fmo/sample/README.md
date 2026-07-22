# md-fmo PBC 再構成 デモ (極小 toy 系)

親ディレクトリ [`../README.md`](../README.md) で説明した
「複合体を引き剥がさずに液滴を切り出す」問題を、
**数十 KB・数秒で再現できる合成系**で体験するサンプル。
実タンパクの座標は使わない (公開データの機微を持ち込まない)。

## 何を示すか

タンパク A + タンパク B + リガンドが **1 つの GROMACS moleculetype に同居し、
フラグメント間に結合が無い**系 (Amber→acpype→GROMACS で溶質が 1 分子に潰れた
状況) を 18 原子で模す。B が周期境界を越えた「生構造」を各手法で組み直す:

| 手法 | 結果 | 理由 |
|---|---|---|
| raw (そのまま) | **NG** | B が境界の向こうのイメージにある |
| `gmx trjconv -pbc whole` | **NG** | A,B,L は単一 moleculetype・相互に無結合。結合をたどる whole では再結合できない |
| `gmx trjconv -pbc cluster` | OK | この小さな系では修復できる (下記の注意) |
| cpptraj `autoimage` | **OK** | prmtop の分子情報から 3 フラグメントを別分子と認識し決定論的に組む |

要点は **「`-pbc whole`/`-pbc mol` では原理的に直らない」** ことと
**「`autoimage` は決定論的に直る」** こと。この 2 つを toy で厳密に確認できる。

### `-pbc cluster` の注意
この極小系では cluster も修復できてしまうが、**実系 (多体 + 溶媒) では
cluster の反復アルゴリズムが初期配置依存になり、正しく接触した複合体を
壊すことがある**。実測 (200 ns 系・20 フレーム) では 11/20 が破綻し、
水和殻に最大 28 Å の空洞が生じた ([`../README.md`](../README.md) の表)。
`autoimage` はこの初期配置依存が無い。

## 実行

```bash
# 必要: gmx, cpptraj (AmberTools), python3 + numpy + scipy + parmed
bash run_demo.sh
```

`make_toy.py` が `toy.gro` (生構造) / `toy_good.gro` (正解) / `toy.top` を生成し、
`run_demo.sh` が prmtop 生成 → 4 手法で組み直し → `../check_contact.py` で判定する。

## 期待される出力

```
###### 生構造 (B が周期境界の向こう) を組み直した結果 ######
W_raw.pdb        A-B= 40.94 A-L=  2.98 B-L= 43.84   *** NG ***
W_whole.pdb      A-B= 40.94 A-L=  2.98 B-L= 43.84   *** NG ***
W_cluster.pdb    A-B=  2.94 A-L=  2.98 B-L=  2.97   OK
W_autoimage.pdb  A-B=  2.94 A-L=  2.98 B-L=  2.97   OK
```

## ファイル

| ファイル | 内容 |
|---|---|
| `make_toy.py` | toy 系 (gro/top) を生成。numpy のみ依存 |
| `run_demo.sh` | prmtop 生成 → 4 手法比較 → 判定 |
| 生成物 (`toy*.gro`, `toy.top`, `toy.prmtop`, `*.pdb`, `*.tpr` 等) | 実行時に作られる。コミットしない |
