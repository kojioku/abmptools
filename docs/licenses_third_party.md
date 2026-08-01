# サードパーティライセンス一覧

ABMPTools が依存するサードパーティライブラリのライセンス情報です。
(本ファイル冒頭は以前 "FCEWS" を含んでいましたが、abmptools の文脈に
合わせて修正しました — FCEWS 本体は別 repo (`fcews-workspace/fcews`) で
別途 LICENSE 管理されます。)

**確認日**: 2026-02-22 (geomopt 系)、2026-05-03 (membrane / amorphous 系)
**確認方法**: `pip show`、`importlib.metadata`、各プロジェクトの GitHub リポジトリ

---

## ABMPTools 本体

| プロジェクト | ライセンス | 所有者 |
|---|---|---|
| **abmptools** | **Apache-2.0** (v1.23.0+; ≤ v1.22.0 は MIT) | Koji Okuwaki |
| **fcews** | *(LICENSE ファイルなし — 要追加)* | — |

abmptools は v1.23.0 以降 **Apache License, Version 2.0** で配布される (LICENSE / NOTICE / CITATION.cff の 3 点セット)。Apache-2.0 移行の理由は Citation 強制 (NOTICE-file attribution) と特許 retaliation による defensive な姿勢の確立。詳細は `CHANGELOG.md` `[Unreleased]` → "License migration" 節を参照。

依存ライブラリの license との互換性 (Q&A での確認結果):

- Apache-2.0 ⇔ **GPL-2.0 only**: ❌ incompatible (FSF 公式)、ただし abmptools は `insane` (GPL-2.0) を **subprocess only** で呼び出すため mere aggregation = 実害なし
- Apache-2.0 ⇔ GPL-2.0-or-later: ✅ compatible (GPL-3.0 経路で吸収)
- Apache-2.0 ⇔ GPL-3.0 / LGPL-3.0+: ✅ compatible
- Apache-2.0 ⇔ LGPL-2.1+ (ase): ✅ compatible (import 利用、binary 同梱なし)
- Apache-2.0 ⇔ MIT / BSD: ✅ compatible

---

## コア依存（全モード共通）

| パッケージ | バージョン（確認時） | ライセンス | 商用利用 | 参照 |
|---|---|---|---|---|
| **numpy** | 1.26.4 | BSD 3-Clause | ✅ 自由 | [github](https://github.com/numpy/numpy/blob/main/LICENSE.txt) |
| **pandas** | — | BSD 3-Clause | ✅ 自由 | [github](https://github.com/pandas-dev/pandas/blob/main/LICENSE) |

---

## abmptools.amorphous — multi-component amorphous builder

| パッケージ | バージョン（確認時） | ライセンス | 商用利用 | 備考 |
|---|---|---|---|---|
| **openff-toolkit** | — | MIT | ✅ 自由 | |
| **openff-interchange** | — | MIT | ✅ 自由 | |
| **openmm** | — | MIT | ✅ 自由 | |
| **rdkit** | — | BSD 3-Clause | ✅ 自由 | |
| **packmol** (binary, Martínez et al.) | — | フリーソフトウェア (公式配布) | ✅ 自由 | 学術 + 商用 OK、引用必須 (J. Comput. Chem. 2009) |
| **ambertools** (AM1-BCC, GAFF, …) | ≥ 22 | LGPL / 各モジュール混在 | ✅ 自由 (商用含む) | 配布物全体は free incl. commercial。AmberTools 公式ライセンス参照 |
| **openff-nagl** (任意、ML 電荷の代替) | — | MIT | ✅ 自由 | AmberTools 不要時のオプション |
| **gromacs** (post-MD scripts) | 2020+ | LGPL-2.1 | ✅ 自由 (リンク・実行) | source 改変配布時は LGPL 準拠 |

---

## abmptools.membrane — peptide-bilayer Umbrella Sampling builder

商用利用 OK な権利のみで構成されており、CGenFF Web server / CHARMM-GUI 等の
有償ライセンス対象には依存しません ([`membrane.md`](membrane.md) のライセンス
ルール参照)。

### MD エンジン・配布物

| パッケージ | バージョン（確認時） | ライセンス | 商用利用 | 備考 |
|---|---|---|---|---|
| **gromacs** | 2021+ (or `nompi_cuda*` for GPU offload) | LGPL-2.1 | ✅ 自由 | grompp / mdrun / wham / trjconv |
| **ambertools** (`tleap`, `packmol-memgen`, `antechamber`, `parmchk2`) | ≥ 22 | 各モジュール混在 (公式配布全体は free incl. commercial) | ✅ 自由 | bilayer 構築、力場割り当て、新規小分子 GAFF2 パラメータ化 |
| **parmed** | 4.3+ | LGPL-2.1+ | ✅ 自由 | AMBER prmtop/inpcrd → GROMACS top/gro 変換 |
| **packmol-memgen** (Schott-Verdugo, AmberTools 同梱) | 2023.2.24 | AmberTools 一部として配布 | ✅ 自由 | bilayer + 水 + イオン + ペプチド初期配置 |

### 力場 (parameter values)

| 力場 | ライセンス | 商用利用 | 備考 |
|---|---|---|---|
| **AMBER ff19SB** (タンパク質) | AmberTools 同梱、free incl. commercial | ✅ 自由 | Tian 2020 |
| **AMBER Lipid21** (脂質) | 同上 | ✅ 自由 | Skjevik 2012 |
| **AMBER GAFF2** (任意小分子) | 同上 | ✅ 自由 | Antechamber 経由 |
| **TIP3P / Joung-Cheatham ions** | 同上 | ✅ 自由 | 標準水・対イオン |
| **CHARMM36 / CHARMM36m** (parameter values) | MacKerell 研公式: "free of charge to academic and industrial researchers" | ✅ 自由 | パラメータ値そのものは free。下記 CGenFF / CHARMM-GUI とは別 |
| **CHARMM36 GROMACS port** (Klauda lab) | 上記 CHARMM36 のライセンスを継承 | ✅ 自由 | `MembraneConfig.charmm_ff_dir` で参照、本パッケージは未同梱 |

### 商用利用 NG の経路 (本パッケージは依存しない)

| 経路 | ライセンス | 備考 |
|---|---|---|
| ❌ **CGenFF Web server** (`cgenff.umaryland.edu`) | 商用は Silcsbio 経由のサブスク必須 | 新規小分子のパラメータ化は AMBER GAFF2 経路で代替 |
| ❌ **CHARMM-GUI** 自動生成 | 商用は別契約必要 | bilayer 構築は packmol-memgen で代替 |

### 任意 / 将来的に検討

| パッケージ | ライセンス | 用途 |
|---|---|---|
| **PyMBAR** | MIT | WHAM の代替・検証用 MBAR (現在 stub、将来本実装予定) |

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
| **MIT** | なし | なし | ✅ | 最も制約が少ない (旧 abmptools ≤ v1.22.0) |
| **BSD 2-Clause** | なし | なし | ✅ | 名称使用制限なし |
| **BSD 3-Clause** | なし | なし | ✅ | 名称使用制限が追加 |
| **Apache 2.0** | なし | なし | ✅ | 特許条項 (retaliation 含む) + NOTICE 保持義務。**本プロジェクトのライセンス (v1.23.0+)** |
| **LGPL-2.1+** | ファイル単位 | ライブラリ自体の改変のみ | ✅ | リンク・import は自由 |
| **LGPL-3.0** | ファイル単位 | ライブラリ自体の改変のみ | ✅ | リンク・import は自由 |
| **MPL-2.0** | ファイル単位 | 同一ファイルの改変のみ | ✅ | pyberny が該当 |

---

## 結論：ABMPTools (Apache-2.0、v1.23.0+) との互換性

**研究・内部利用の範囲では全て問題ありません。**

配布・公開時の注意点：

1. **pyberny (MPL-2.0)**: pyberny のソースファイル自体を改変して配布する場合は
   改変ファイルのソースを公開する必要があります。
   import して使うだけなら Apache-2.0 のプロジェクトに組み込んでも問題ありません。
   また、`solver="berny"` はオプションであり、デフォルトは `geometric`（BSD 3-Clause）です。

2. **ase, dftd3 (LGPL)**: ライブラリ本体を改変して配布する場合は LGPL 準拠が必要。
   import して使うだけなら制約はありません。

3. **著作権表示と NOTICE**: Apache-2.0 の §4(d) 規定により、abmptools を再配布
   する場合は LICENSE と NOTICE をすべて保持する必要があります。MIT / BSD /
   Apache ライセンスの依存ライブラリを同梱・配布する場合は、各ライブラリの
   著作権表示とライセンス文の同梱も必要です (`docs/licenses_third_party.md`
   = 本ファイル を維持しておけば実質的に充足できます)。

4. **特許 retaliation (Apache-2.0 §3)**: abmptools またはその Contribution
   に対して特許訴訟を起こした場合、相手側に与えていた特許ライセンスは自動失効
   します。defensive な構造として機能。

5. **insane (GPL-2.0 only) との incompatibility**: Apache-2.0 ⇔ GPL-2.0 only
   は厳密には incompatible (§3 patent retaliation / §4(d) NOTICE 保持義務 が
   GPL-2.0 §6 "no further restrictions" に違反) だが、insane は subprocess
   only で呼び出すため mere aggregation = 実害なし。GPL-2.0-or-later の場合は
   GPL-3.0 経路で吸収されるため問題なし。

---

## 論文投稿時の引用（ライセンスとは別途必要）

### `qmopt` 機能

| ソフトウェア | 引用文献 |
|---|---|
| **PySCF** | Q. Sun et al., *WIREs Comput. Mol. Sci.* **2018**, *8*, e1340; *J. Chem. Phys.* **2020**, *153*, 024109 |
| **geomeTRIC** | L.-P. Wang, C. C. Song, *J. Chem. Phys.* **2016**, *144*, 214108 |
| **B3LYP 汎関数** | A. D. Becke, *J. Chem. Phys.* **1993**, *98*, 5648; C. Lee, W. Yang, R. G. Parr, *Phys. Rev. B* **1988**, *37*, 785 |
| **D3(BJ) 分散補正** | S. Grimme et al., *J. Chem. Phys.* **2010**, *132*, 154104; S. Grimme et al., *J. Comput. Chem.* **2011**, *32*, 1456 |
| **def2-SVP 基底関数** | F. Weigend, R. Ahlrichs, *Phys. Chem. Chem. Phys.* **2005**, *7*, 3297 |

### `membrane` 機能

| ソフトウェア / 力場 | 引用文献 |
|---|---|
| **GROMACS** | M. J. Abraham et al., *SoftwareX* **2015**, *1-2*, 19 |
| **GROMACS GPU offload** | S. Páll et al., *J. Chem. Phys.* **2020**, *153*, 134110 |
| **AmberTools** | D. A. Case et al., *J. Comput. Chem.* **2005**, *26*, 1668 |
| **ff19SB (タンパク質)** | C. Tian et al., *J. Chem. Theory Comput.* **2020**, *16*, 528 |
| **Lipid21 (脂質)** | A. A. Skjevik et al., *J. Phys. Chem. B* **2012**, *116*, 11124 |
| **TIP3P-JC イオン** | I. S. Joung, T. E. Cheatham III, *J. Phys. Chem. B* **2008**, *112*, 9020 |
| **packmol** | L. Martínez et al., *J. Comput. Chem.* **2009**, *30*, 2157 |
| **packmol-memgen** | S. Schott-Verdugo, H. Gohlke, *J. Chem. Inf. Model.* **2019**, *59*, 2522 |
| **WHAM (gmx wham)** | J. S. Hub et al., *J. Chem. Theory Comput.* **2010**, *6*, 3713 |
| **CHARMM36 (使用時)** | R. B. Best et al., *J. Chem. Theory Comput.* **2012**, *8*, 3257; J. B. Klauda et al., *J. Phys. Chem. B* **2010**, *114*, 7830; R. W. Pastor et al., *Chem. Phys. Lett.* **2011**, *517*, 24 |

## abmptools.formulation — AA peptide-formulation builder (1.30.0+)

Hossain et al. 2023 (*Nanoscale* 15, 19180) を **商用利用可能な OSS のみで再現**する builder + analysis stack。

### 採用力場 / ツール

| 役割 | ツール / 力場 | ライセンス |
|---|---|---|
| ペプチド (natural AA) | AMBER ff14SB | パブリックドメイン |
| 小分子 (enhancer / bile salt) | GAFF2 + AM1-BCC | パブリックドメイン |
| 小分子パラメータ生成 | acpype | GPL-3.0 (subprocess のみ = mere aggregation) |
| AmberTools (tleap / antechamber / parmchk2 / sqm) | パブリックドメイン (商用可) |
| 水 | TIP3P (`leaprc.water.tip3p`) | パブリックドメイン |
| イオン | Joung-Cheatham Na/Cl | パブリックドメイン |
| パッキング | packmol | 自由配布 (`amorphous` 経由で reuse) |
| MD | GROMACS | LGPL-2.1 (subprocess のみ) |

### MDAnalysis (GPL-2.0) の lazy-import 方針

Analysis モジュール (`abmptools/formulation/analysis/aggregate_transition.py`,
`contact_map.py`) は MDAnalysis を **関数本体内で lazy import** する。

- **abmptools 配布 wheel には MDAnalysis を同梱しない**
- `pip install abmptools[formulation-analysis]` で user が opt-in 導入
- ImportError 時は install hint メッセージを emit
- これにより abmptools 本体は **Apache-2.0 維持** (FSF の GPL FAQ § "mere
  aggregation" 解釈に沿う; 既存 `abmptools/amorphous/trajectory_ingest.py`
  L115-150 の precedent)

### 不採用 (academic license のため)

- **CHARMM36 + CGenFF**: Hossain 2023 の原典力場だが、CHARMM-GUI / CGenFF
  Web service が academic license のため `formulation` では使わない。
  論文の絶対 PMF 値とは ~5-15 kJ/mol drift する可能性 (定性傾向は保つ)。

### 引用文献

| ソフトウェア | 引用 |
|---|---|
| **acpype** | A. W. Sousa da Silva, W. F. Vranken, *BMC Res. Notes* **2012**, *5*, 367 |
| **GAFF2 / AM1-BCC** | J. Wang et al., *J. Comput. Chem.* **2004**, *25*, 1157; A. Jakalian et al., *J. Comput. Chem.* **2002**, *23*, 1623 |
| **ff14SB** | J. A. Maier et al., *J. Chem. Theory Comput.* **2015**, *11*, 3696 |
| **MDAnalysis** | R. J. Gowers et al., *Proc. 15th SciPy Conf.*, **2016**, 98-105 |
| **NetworkX** | A. Hagberg, P. Swart, D. Schult, *Proc. 7th SciPy Conf.*, **2008**, 11-15 |
