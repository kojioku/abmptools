# Peptide Builders — abmptools の 3 系統選択ガイド

abmptools には現在「**残基配列 (sequence) からペプチドを自動構築する**」機能が
3 系統あります。本ドキュメントは各系統の特徴を比較し、用途に応じてどれを使うべきかの
選択ガイドを示します。

## 一覧

| 系統 | サブパッケージ | バージョン | 解像度 | 系の範囲 | バックエンド |
|---|---|---|---|---|---|
| **A. AA peptide-membrane US** | [`abmptools.membrane`](membrane.md) | 1.17.0+ | 全原子 (AA) | peptide + 脂質膜 + 水 + 塩 | `tleap` (AMBER) / `pdb2gmx` (CHARMM36) + `packmol-memgen` |
| **B. CG peptide** | [`abmptools.cg.peptide`](cg_peptide.md) | 1.18.0+ | Martini 3 (CG) | peptide のみ + 水 + 塩 | `martinize2` + `gmx` (+ `tleap` 推奨) |
| **C. CG peptide-membrane US** | [`abmptools.cg.membrane`](cg_membrane.md) | 1.19.0+ | Martini 3 (CG) | peptide + 脂質膜 + 水 + 塩 | `cg.peptide` sub-call + `insane` (bilayer) |

## 入力フィールドの比較

| 項目 | A. AA membrane | B. CG peptide | C. CG membrane |
|---|---|---|---|
| `sequence` | 1-letter ("AAAAA" 等) | 1-letter ("KGG" 等) | 1-letter ("KGG" 等) |
| **N/C cap** | ACE / NME 指定可 (`cap_n` / `cap_c`) | 非対応 | 非対応 |
| **non-standard residue** (HIP/CYM/HSD/HSE 等) | backend 依存で利用可 | reject (validation で error) | reject |
| **複数ペプチド混合** | `n_copies` で複数本同種 | `peptides: List[PeptideSpec]` で異種混合 OK | 1 種のみ (v1.19) |
| **脂質膜** | POPC ほか lipid table から指定 (混合 OK) | なし | POPC / DOPC / DPPC ほか insane 対応 lipid (v1.19 は単一種) |
| **PMF (umbrella sampling)** | あり (13-31 windows × 1-5 ns AA) | なし | あり (13-31 windows × 1-5 ns CG、AA の 30-100× 速い) |
| **入力スクリプト** | `MembraneConfig` JSON / YAML | `PeptideBuildConfig` JSON / YAML | `MembraneCGBuildConfig` JSON / YAML |

## 用途別の選択ガイド

### 「ペプチドを単独で水中に入れて MD したい」

→ **B. `abmptools.cg.peptide`** (CG、軽量) または独自に AA を組む。
abmptools は AA peptide-only ビルダーは現状提供せず (membrane なし AA は
ユーザーが直接 `tleap` を使うのが推奨)。

```bash
python -m abmptools.cg.peptide example > kgg.json
python -m abmptools.cg.peptide build --config kgg.json --ff-dir ./ff -o ./out
bash out/run/run.sh
```

### 「peptide が脂質膜を透過する PMF を計算したい (AA、FF 比較や原子解像度の議論)」

→ **A. `abmptools.membrane`** (CHARMM36 / Lipid21 backend、wall ~3 時間 GPU)。
`peptide.sequence: "AAAAA"` + `cap_n: "ACE"` + `cap_c: "NME"` でビルド開始。

```bash
# AMBER backend
python -m abmptools.membrane example > polyala_amber.json
# CHARMM36 backend (config に "backend": "charmm36" と "charmm_ff_dir")
python -m abmptools.membrane example > polyala_charmm.json
```

### 「peptide が脂質膜を透過する PMF を、もっと安く / 長くサンプリングしたい (CG)」

→ **C. `abmptools.cg.membrane`** (Martini 3、wall 5-45 分 CPU)。
sequence + lipid 種を指定すれば 1 コマンドで build → MD → WHAM まで。AA の
30-100× 速いので long peptide や parameter sweep (k, spacing, sampling time) に向く。

```bash
python -m abmptools.cg.membrane example > kgg.json
python -m abmptools.cg.membrane build --config kgg.json --ff-dir ./ff -o ./out
bash out/run.sh
```

## 解像度別の比較

| 観点 | A. AA membrane | B/C. CG (Martini 3) |
|---|---|---|
| 計算コスト (KGG-POPC、13-31 window × 1-5 ns) | wall ~3 時間 GPU | wall **5-45 分** CPU |
| 原子レベルの相互作用 | あり (FF 比較・構造詳細議論可能) | なし (4:1 mapping) |
| 大規模 peptide (>20 残基) | コスト爆発 | 実用的 |
| 文献値との直接比較 | AMBER ff19SB / CHARMM36 文献多数 | Martini 3 文献は急増中 |
| FF 選択肢 | AMBER ff19SB / Lipid21 / TIP3P or CHARMM36 / TIP3P-modified | Martini 3001 のみ (current) |
| **典型的な用途** | 機構研究、FF 比較、原子解像度の議論 | スクリーニング、長 peptide、parameter sweep、定性的傾向 |

## 共通する設計

3 系統は abmptools 流儀で設計が揃えられています:

- `@dataclass` ベースの config (JSON / YAML 往復)
- `python -m abmptools.<sub>` で `example` / `validate` / `build` の CLI 提供
- `output_dir/run.sh` を生成 → `bash run.sh` で全 MD stage 自動実行
- `forcefield_check.py` で外部依存性を事前検証
- subprocess only 戦略で外部ツールのソース改変 / 同梱なし (license 問題回避)

## ライセンス上の留意

| 系統 | バンドル | サブプロセス依存 | 商用利用 |
|---|---|---|---|
| A. AA membrane | abmptools 本体 (MIT) | `packmol-memgen` (AmberTools, LGPL) / `tleap` (LGPL) / `pdb2gmx`+`gmx` (LGPL) / CHARMM36 ff (academic ok / industrial 別途確認) | ✅ (CHARMM-GUI を使わない設計、CGenFF も非依存) |
| B. CG peptide | abmptools 本体 (MIT) | `martinize2` (Apache-2.0) / `gmx` (LGPL) / `tleap` (LGPL、推奨任意) | ✅ Martini 3 ITP は cgmartini.nl 配布物 (license 未明記、ユーザー自取得) |
| C. CG membrane | abmptools 本体 (MIT) | 上記 B + `insane` (**GPL-2.0**) | ✅ insane は subprocess only (mere aggregation、GPL に感染しない、FSF 公式見解) |

## 各系統の最小サンプル config

### A. AA membrane (poly-Ala 5、CHARMM36)

```json
{
  "backend": "charmm36",
  "lipids": [{"resname": "POPC", "n_per_leaflet": 32}],
  "peptide": {
    "name": "aa5", "sequence": "AAAAA",
    "cap_n": "ACE", "cap_c": "NME", "n_copies": 1,
    "initial_z_offset_nm": 3.0
  },
  "ions": {"cation": "Na+", "anion": "Cl-",
           "salt_concentration_M": 0.15, "neutralize": true},
  "box_xy_nm": 6.0,
  "umbrella": {
    "z_min_nm": -1.5, "z_max_nm": 1.5,
    "window_spacing_nm": 0.25,
    "force_constant_kj_mol_nm2": 1000.0,
    "window_nsteps": 500000, "window_dt_ps": 0.002
  },
  "charmm_ff_dir": "/path/to/charmm36-feb2026_cgenff-5.0.ff"
}
```

完全なサンプル: [`sample/membrane/charmm_phaseD/input/config_phaseD.json`](../sample/membrane/charmm_phaseD/input/config_phaseD.json)

### B. CG peptide (KGG ×5、water box)

```json
{
  "peptides": [{"name": "kgg", "sequence": "KGG", "count": 5}],
  "box_size_nm": 8.0,
  "solvent_enabled": true,
  "neutralize": true,
  "nacl_molar": 0.15,
  "temperature": 310.0,
  "martini_itp_dir": "./ff",
  "output_dir": "./out"
}
```

### C. CG membrane (KGG ×1、POPC bilayer)

```json
{
  "lipids": [{"resname": "POPC", "n_per_leaflet": 64}],
  "peptide": {"name": "kgg", "sequence": "KGG",
              "initial_z_offset_nm": 3.0},
  "insane_d_nm": 8.0, "box_z_nm": 14.0,
  "nacl_molar": 0.15,
  "umbrella": {
    "z_min_nm": -1.5, "z_max_nm": 1.5,
    "window_spacing_nm": 0.10,
    "force_constant_kj_mol_nm2": 500.0,
    "window_nsteps": 250000
  },
  "martini_itp_dir": "./ff"
}
```

完全なサンプル: [`sample/cg_membrane/kgg_popc_production.json`](../sample/cg_membrane/kgg_popc_production.json)

## 関連リソース

- [`docs/membrane.md`](membrane.md) — A 系統 (AA membrane) subsystem reference
- [`docs/cg_peptide.md`](cg_peptide.md) — B 系統 (CG peptide) subsystem reference
- [`docs/cg_membrane.md`](cg_membrane.md) — C 系統 (CG membrane) subsystem reference
- [`docs/tutorial_membrane_us.md`](tutorial_membrane_us.md) — A 系統 step-by-step
- [`docs/tutorial_cg_membrane_us.md`](tutorial_cg_membrane_us.md) — C 系統 step-by-step
- [`docs/overview.md`](overview.md) — abmptools 全体の subpackage 一覧

## バージョン履歴

| version | 日付 | 変更 |
|---|---|---|
| (本ドキュメント新設) | 2026-05-06 | 3 系統横断の選択ガイドとして新規作成 |
