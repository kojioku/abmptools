# md-fmo PBC 再構成 デモ (溶媒つき toy 系)

親ディレクトリ [`../README.md`](../README.md) で説明した
「複合体を引き剥がさず・水和殻を壊さずに液滴を切り出す」問題を、
**約 1000 原子・数秒で再現できる合成系**で体験するサンプル。
実タンパクの座標は使わない (公開データの機微を持ち込まない)。

## 何を示すか

タンパク A + タンパク B + リガンドが **1 つの GROMACS moleculetype に同居し、
フラグメント間に結合が無い**系 (Amber→acpype→GROMACS で溶質が 1 分子に潰れた
状況) を、水 (WAT) で溶媒和して模す。箱は複合体に対しやや小さめ (tight box) で、
複合体が周期境界をまたぐ「生構造」を、実系の `2_optmask-frame.sh` と同じ 4 段
(`whole → cluster → mol/compact`) と cpptraj `autoimage` で組み直して比較する。

`run_demo.sh` の出力例:

```
good input : A-B=  3.20 A-L=  1.67 B-L=  1.67   water= 331 max_gap=  3.67   OK
whole      : A-B= 18.68 A-L= 31.66 B-L= 22.21   water= 331 max_gap=  3.55   *** NG ***
cluster    : A-B=  3.20 A-L=  1.67 B-L=  1.67   water= 331 max_gap= 13.90   *** NG ***
mol/compact: A-B=  3.20 A-L=  1.67 B-L=  1.67   water= 331 max_gap=  3.67   OK
autoimage  : A-B=  3.20 A-L=  1.67 B-L=  1.67   water= 331 max_gap=  3.67   OK
```

| 手法 | 接触 | 水和殻 (max_gap) | 判定 | 理由 |
|---|---|---|---|---|
| raw / `-pbc whole` | **NG** (18–32 Å) | — | NG | A,B,L は単一 moleculetype・相互に無結合。結合をたどる whole では再結合できない |
| `-pbc cluster` | OK (2–3 Å) | **NG (13.9 Å の穴)** | **NG** | 接触は戻すが**水和殻を置き去りにする (溶質が脱水)** ← 溶媒系特有の失敗 |
| `-pbc mol/compact` | OK | OK | OK | 後段でようやく水を戻す (**全段が成功すればの話**) |
| cpptraj `autoimage` | OK | OK | **OK** | prmtop の分子情報から 1 手・決定論的に接触も水和殻も回復 |

### 読み取れること
1. `-pbc whole`/`-pbc mol` だけでは複合体を組み直せない (単一 moleculetype)。
2. **`-pbc cluster` は接触を戻せても水和殻に穴を開ける** (溶質が脱水)。
   正しい液滴になるかは後段の `-pbc mol/compact/center` が成功するかに依存する
   多段依存の危うい経路。
3. **実系 (多体 + 大きなタンパク) では `-pbc cluster` が接触自体も壊す**。
   実測 (200 ns 系・20 フレーム) で 11/20 が破綻し、水和殻に最大 28 Å の空洞が
   生じた ([`../README.md`](../README.md) の表)。この toy でも cluster 段で
   同種の脱水が起きることが確認できる。
4. `autoimage` は接触も水和殻も **1 手で決定論的に**回復する。

## 実行

```bash
# 必要: gmx, cpptraj (AmberTools), python3 + numpy + scipy + parmed
bash run_demo.sh
```

`make_toy.py` が `toy.gro` (生構造) / `toy_good.gro` (正解) / `toy.top` を生成し、
`run_demo.sh` が prmtop 生成 → 実 4 段パイプライン + autoimage → `../check_contact.py`
で各段を判定する。

## ファイル

| ファイル | 内容 |
|---|---|
| `make_toy.py` | 溶媒つき 3 フラグメント複合体 (A,B,リガンド + 水) を生成。numpy のみ依存 |
| `run_demo.sh` | prmtop 生成 → whole/cluster/mol-compact + autoimage を比較 → 判定 |
| 生成物 (`toy*.gro`, `toy.top`, `toy.prmtop`, `*.pdb`, `*.tpr`, `*.ndx` 等) | 実行時に作られる。コミットしない |
