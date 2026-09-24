# QMOptimizerPySCF — PySCF QM Geometry Optimiser

`abmptools.geomopt.QMOptimizerPySCF` は PySCF を使った DFT geometry optimiser です。
MACE・OpenFF optimizer と同等のインターフェースで使えます。

## 依存ライブラリ

| パッケージ | 必須/任意 | ライセンス | バージョン（動作確認済） | インストール |
|---|---|---|---|---|
| `pyscf` | **必須** | Apache 2.0 | 2.12.1 | `pip install pyscf` |
| `geometric` | 必須（推奨） | BSD 3-Clause | 1.1 | `pip install geometric` |
| `pyberny` | 任意 | MPL-2.0 | 0.6.3 | `pip install pyberny`\* |
| `simple-dftd3` | 任意 | MIT | — | `pip install simple-dftd3` |
| `dftd3` | 任意 | LGPL-3.0 | — | `pip install dftd3` |

\* pyberny 0.6.3 は setuptools 82+ 環境で `pkg_resources` インポートエラーが発生する
既知の問題があります（2026-02 確認）。`geometric` を推奨します。

pyscf と geometric（または pyberny）は必須です。

分散補正ライブラリが入っていない場合は警告を出して dispersion なしで実行されます。
**この分岐は結果の見た目を変えません**（計算は最後まで走り、終了コードも 0 です）ので、
D3 を効かせたい場合は次のどちらかで確認してください。

- `opt_results.jsonl` の **`dispersion_applied`** —— `false` なら D3 は付いていません
- 出力 xyz の 2 行目 —— D3 が付いたときだけ `B3LYP-D3BJ/def2-SVP` のように
  汎関数名にハイフンで続きます。付かなかったときは `B3LYP/def2-SVP` と素のまま出ます

`simple-dftd3` を入れる場合、abmptools が使う入口は `dftd3.pyscf.energy` です
（`DFTD3Model` という名前のクラスは存在しません）。

## 溶媒中の最適化

`solvent` に `"pcm"` / `"cpcm"` / `"iefpcm"` / `"ddcosmo"` / `"smd"` を渡すと
PySCF の対応する連続誘電体モデルを被せます（既定は `"none"` = 気相）。
誘電率は `solvent_eps` に**数値**または**溶媒名**（`water`, `methanol`,
`ethanol`, `acetone`, `dichloromethane`, `thf`, `chloroform`, `toluene`,
`benzene`, `cyclohexane`, `hexane`）で指定します。省略するとモデル既定値です。

```python
opt = QMOptimizerPySCF(functional="B3LYP", basis="def2-SVP",
                       dispersion="d3bj",
                       solvent="cpcm", solvent_eps="water")
```

**イオン・双性イオンでは気相最適化が構造を壊します。** 実測（B3LYP-D3BJ/def2-SVP）
では、グリシン双性イオンを気相で最適化すると N から O へ陽子が移り、
最短 O–H が 1.640 Å → **0.991 Å**（O–H 結合が生成）になって中性型へ崩壊しました。
CPCM(water) では 1.712 Å で双性イオンのまま保たれます。多価アルコールが
分子内水素結合で折り畳む問題も同じ性質のものです。

分散補正と同じく、**モデルを付けられなかった場合は警告を出して気相で続行します**。
結果の `solvent_applied` と、出力 xyz の 2 行目（付いたときだけ
`B3LYP-D3BJ/def2-SVP [CPCM,eps=78.3553]` のように角括弧が付く）で確認できます。

## 構造ごとの電荷・スピン

断片集合は電荷が揃わないことが多く（脂質の頭部なら choline が +1、
リン酸ジエステルが −1、アルキル鎖が 0）、バッチ全体で 1 つの `charge` では
足りません。**入力 xyz のコメント行に `charge=-1` や `spin=2` と書くと、
その構造だけインスタンス既定値を上書きします。**

```
13
@@E_diMe-Pho@@ C2H6O4P | charge=-1
P   1.234567  ...
```

`spin` は不対電子数（2S）で、コンストラクタの引数と同じ意味です。
実際に使われた値は結果の `charge` / `spin` に入ります。

ライセンス詳細・互換性の考察は [licenses_third_party.md](./licenses_third_party.md) を参照してください。

## 基本的な使い方

### xyz ファイルの最適化

```python
from abmptools.geomopt import QMOptimizerPySCF

opt = QMOptimizerPySCF()           # デフォルト: B3LYP/def2-SVP/D3(BJ)
result = opt.optimize("water.xyz", "water_opt.xyz")

print(result["energy"])            # eV
print(result["energy_hartree"])    # Ha
print(result["converged"])         # bool
print(result["steps"])             # 最適化ステップ数
print(result["out_xyz"])           # 出力ファイルの絶対パス
```

### PDB ファイルの最適化

```python
# PDB 入力でも出力は xyz 形式
result = opt.optimize("molecule.pdb", "molecule_opt.xyz")
```

> **注意**: PDB ファイルは ATOM/HETATM レコードの要素列（cols 77-78）が必要です。
> 空欄の場合は原子名から推定しますが、推定できない場合は `ValueError` を発生させます。

### パラメータのカスタマイズ

```python
opt = QMOptimizerPySCF(
    functional="PBE0",         # 汎関数（PySCF が受け付ける文字列）
    basis="def2-TZVP",         # 基底関数
    dispersion="d3bj",         # "d3bj" / "d3" / "none"
    charge=0,                  # 分子電荷
    spin=0,                    # 不対電子数（2S）
    max_steps=200,             # 最大最適化ステップ
    solver="geometric",        # "geometric" or "berny"
    verbose=3,                 # PySCF 冗長度（0-9）
)
result = opt.optimize("in.xyz", "out.xyz")
```

## 戻り値のキー

| キー | 型 | 内容 |
|---|---|---|
| `energy` | float | 最終エネルギー（eV） |
| `energy_hartree` | float | 最終エネルギー（Ha） |
| `steps` | int | 最適化ステップ数 |
| `converged` | bool | 収束したか |
| `out_xyz` | str | 出力 xyz ファイルの絶対パス |

## サンプルファイル

`sample/qmopt/` に小分子サンプルを用意しています。

```bash
cd abmptools
python -c "
from abmptools.geomopt import QMOptimizerPySCF
opt = QMOptimizerPySCF(basis='sto-3g', dispersion='none', verbose=0)
r = opt.optimize('sample/qmopt/water.xyz', '/tmp/water_opt.xyz')
print(r)
"
```

## xyz ファイルフォーマット

```
3
water molecule
O   0.000000   0.000000   0.119748
H   0.000000   0.756950  -0.478993
H   0.000000  -0.756950  -0.478993
```

- 1行目: 原子数
- 2行目: コメント（任意）
- 3行目以降: `<元素記号> <x> <y> <z>` (Å)

## 既存 optimizer との比較

| 項目 | `MacePdbOptimizer` | `OpenFFOpenMMMinimizer` | `QMOptimizerPySCF` |
|---|---|---|---|
| 入力 | PDB | PDB | **xyz / PDB** |
| 出力 | PDB | PDB | **xyz** |
| エネルギー単位 | eV | kJ/mol | **eV + Ha** |
| ポテンシャル | ML (MACE) | 古典力場 | **DFT (B3LYP等)** |
| 主な用途 | 大系の高速最適化 | 力場適用可能系 | 精度重視の小系 |
| デフォルト精度 | ML 精度 | FF 精度 | **DFT 精度 + D3(BJ)** |

## 拡張可能な設計

- **汎関数**: `functional="PBE0"`, `"M06-2X"` 等 PySCF 対応のものすべて
- **基底関数**: `basis="6-31G*"`, `"def2-TZVP"`, `"cc-pVDZ"` 等
- **電荷・スピン**: `charge=-1`, `spin=2` でイオン・ラジカルにも対応
- **収束条件**: `conv_params={"convergence_grms": 1e-4}` でカスタマイズ可能
- **開殻系**: `spin > 0` の場合は自動的に UKS を使用
