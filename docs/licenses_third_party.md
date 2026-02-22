# サードパーティライセンス一覧

ABMPTools および FCEWS が依存するサードパーティライブラリのライセンス情報です。

**確認日**: 2026-02-22
**確認方法**: `pip show`、`importlib.metadata`、各プロジェクトの GitHub リポジトリ

---

## ABMPTools 本体

| プロジェクト | ライセンス | 所有者 |
|---|---|---|
| **abmptools** | BSD 2-Clause | Koji Okuwaki |
| **fcews** | *(LICENSE ファイルなし — 要追加)* | — |

---

## コア依存（全モード共通）

| パッケージ | バージョン（確認時） | ライセンス | 商用利用 | 参照 |
|---|---|---|---|---|
| **numpy** | 1.26.4 | BSD 3-Clause | ✅ 自由 | [github](https://github.com/numpy/numpy/blob/main/LICENSE.txt) |
| **pandas** | — | BSD 3-Clause | ✅ 自由 | [github](https://github.com/pandas-dev/pandas/blob/main/LICENSE) |

---

## abmptools.geomopt — MACE バックエンド（`pdbopt --backend mace`）

| パッケージ | バージョン（確認時） | ライセンス | 商用利用 | 備考 |
|---|---|---|---|---|
| **ase** | — | LGPL-2.1+ | ✅ 使用・リンク自由 | ASE 自体の改変配布は LGPL 準拠が必要 |
| **mace-torch** | — | MIT | ✅ 自由 | |
| **torch** (PyTorch) | — | BSD-style | ✅ 自由 | Facebook/Meta 製 modified BSD |

---

## abmptools.geomopt — OpenFF バックエンド（`pdbopt --backend openff`）

| パッケージ | バージョン（確認時） | ライセンス | 商用利用 | 備考 |
|---|---|---|---|---|
| **openmm** | — | MIT | ✅ 自由 | |
| **openff-toolkit** | — | MIT | ✅ 自由 | |
| **rdkit** | — | BSD 3-Clause | ✅ 自由 | |

---

## abmptools.geomopt — QM バックエンド（`qmopt` / `QMOptimizerPySCF`）

| パッケージ | バージョン（確認時） | ライセンス | 商用利用 | 備考 |
|---|---|---|---|---|
| **pyscf** | 2.12.1 | Apache 2.0 | ✅ 自由 | 特許条項あり（Apache 特有）。競合他社への特許攻撃に対する反撃条項 |
| **geometric** | 1.1 | BSD 3-Clause | ✅ 自由 | `pip show` では License フィールド空。GitHub で BSD 3-Clause 確認 |
| **pyberny** | 0.6.3 | MPL-2.0 | ✅ 自由 | ファイル単位コピーレフト。pyberny ファイル自体を改変・配布する場合はソース公開が必要。import して使うだけなら BSD 2-Clause と共存可 |
| **simple-dftd3** | — | MIT | ✅ 自由 | Python バインディングおよび C ライブラリ本体ともに MIT |
| **dftd3** (代替) | — | LGPL-3.0 | ✅ 使用・リンク自由 | Grimme グループ製オリジナル。ライブラリ本体の改変配布は LGPL 準拠が必要 |
| **networkx** | 3.4.2 | BSD 3-Clause | ✅ 自由 | geometric の依存として自動インストール |
| **six** | 1.17.0 | MIT | ✅ 自由 | geometric の依存として自動インストール |

---

## ライセンス種別の概要

| ライセンス | コピーレフト | 改変の公開義務 | 商用利用 | 特記事項 |
|---|---|---|---|---|
| **MIT** | なし | なし | ✅ | 最も制約が少ない |
| **BSD 2-Clause** | なし | なし | ✅ | 本プロジェクトのライセンス |
| **BSD 3-Clause** | なし | なし | ✅ | 名称使用制限が追加 |
| **Apache 2.0** | なし | なし | ✅ | 特許条項あり（ユーザー保護） |
| **LGPL-2.1+** | ファイル単位 | ライブラリ自体の改変のみ | ✅ | リンク・import は自由 |
| **LGPL-3.0** | ファイル単位 | ライブラリ自体の改変のみ | ✅ | リンク・import は自由 |
| **MPL-2.0** | ファイル単位 | 同一ファイルの改変のみ | ✅ | pyberny が該当 |

---

## 結論：ABMPTools (BSD 2-Clause) との互換性

**研究・内部利用の範囲では全て問題ありません。**

配布・公開時の注意点：

1. **pyberny (MPL-2.0)**: pyberny のソースファイル自体を改変して配布する場合は
   改変ファイルのソースを公開する必要があります。
   import して使うだけなら BSD 2-Clause のプロジェクトに組み込んでも問題ありません。
   また、`solver="berny"` はオプションであり、デフォルトは `geometric`（BSD 3-Clause）です。

2. **ase, dftd3 (LGPL)**: ライブラリ本体を改変して配布する場合は LGPL 準拠が必要。
   import して使うだけなら制約はありません。

3. **著作権表示**: MIT / BSD / Apache ライセンスのライブラリを同梱・配布する場合は
   各ライブラリの著作権表示とライセンス文の同梱が必要です。

---

## 論文投稿時の引用（ライセンスとは別途必要）

`qmopt` 機能で成果を発表する際は以下の論文を引用してください。

| ソフトウェア | 引用文献 |
|---|---|
| **PySCF** | Q. Sun et al., *WIREs Comput. Mol. Sci.* **2018**, *8*, e1340; *J. Chem. Phys.* **2020**, *153*, 024109 |
| **geomeTRIC** | L.-P. Wang, C. C. Song, *J. Chem. Phys.* **2016**, *144*, 214108 |
| **B3LYP 汎関数** | A. D. Becke, *J. Chem. Phys.* **1993**, *98*, 5648; C. Lee, W. Yang, R. G. Parr, *Phys. Rev. B* **1988**, *37*, 785 |
| **D3(BJ) 分散補正** | S. Grimme et al., *J. Chem. Phys.* **2010**, *132*, 154104; S. Grimme et al., *J. Comput. Chem.* **2011**, *32*, 1456 |
| **def2-SVP 基底関数** | F. Weigend, R. Ahlrichs, *Phys. Chem. Chem. Phys.* **2005**, *7*, 3297 |
