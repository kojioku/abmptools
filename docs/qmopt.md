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
