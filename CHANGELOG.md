# Changelog

## [Unreleased]

## [2.2.0] - 2026-06-24

### Fixed — `abmptools.cg.dpd` build-udf が cognac でロードできない不具合

- **重大**: 旧 `udf_writer.py` (positional plain-text) の出力は **UDFManager で
  ロード不可** (`Pair_Interaction` を `Molecular_Attributes` 配下に置く /
  `Interaction_Site_Type[]` 未定義で cognac112 スキーマ非準拠)。実機で
  `RuntimeError: ... near "0.0"` を確認。ユニットテスト 66 件は文字列生成のみ検証
  でロードを見ておらず見逃していた
- 新規 `cg/dpd/udf_writer_udfm.py` (`write_dpd_udf_udfm`): **UDFManager の
  named-path put** で組み立て、cognac112 スキーマ準拠を保証。`Pair_Interaction`
  を正しく `Interactions[]` 配下へ、`Molecular_Attributes.Interaction_Site_Type[]`
  を定義。`CGDpdBuilder.build_udf` / CLI `build-udf` をこちらに切替
- **設計変更**: build-udf は **UDFManager(OCTA) 必須**に (旧「UDFManager 非依存
  plain text」設計を廃止)。DPD UDF は元々 OCTA で使うもので、生成できても
  ロード不可では無意味なため
- テスト刷新: `tests/cg_dpd/test_udf_writer.py` を **loadability ベース**に書換
  (UDFManager でロード → `Interactions.Pair_Interaction[]` size / `DPD.a` 等を検証)。
  `conftest.requires_cognac` で OCTA 非導入環境では skip。cg_dpd 全体 67 passed
  (OCTA 有り環境)

### Added — `abmptools.cg.dpd assign-aij` (既存 DPD UDF への aij 割り当て)

- 新規 `cg/dpd/aij_assign.py`: 既存 Cognac DPD UDF の
  `Interactions.Pair_Interaction[].DPD.a` を aij.dat で割り当て直す
  (build-udf の新規生成と別。既存 UDF をパッチ):
  - `Site1_Name` / `Site2_Name` で粒子名ペアを読み、aij.dat と **順不同照合**
  - aij モード=値直入れ、chi モード=`a = aii + χ/0.306` (= aii + 3.268·χ、
    `AijMatrix.to_a_values()`)
  - `UDFManager` round-trip (`abmptools.udfcharge` と同方式、`float()` cast で
    numpy silent-0 回避)。`Potential_Type=="DPD"` のみ対象 (default)
  - `build_a_lookup` / `match_aij_to_pairs` は UDFManager 非依存で単体テスト可能
- CLI `python -m abmptools.cg.dpd assign-aij --udf X.udf --aij aij.dat [--aii 25.0]
  [--output Y.udf]`
- テスト `tests/cg_dpd/test_aij_assign.py` 7 件 (chi→a 変換 / 順不同照合 / 未照合
  検出 / モック UDFManager で I/O 配線)。cg_dpd 全体 66 passed

### Changed (docs) — cg_segmenter DPD 経路の表記を実態に整合 + 公開スコープ明記

- `cg_segmenter` の DPD 入力生成 (`dpdgen_exporter`) の下流が **`abmptools.cg.dpd`**
  (自前実装、外部 dpdgen / UDFManager 非依存) であることを docstring / CLI 出力 /
  Jupyter UI / docs / README に統一。誤解を生む「外部 `makeudf_dpd.py` でそのまま
  処理できる」表記を `python -m abmptools.cg.dpd build-udf ...` に修正
- `docs/cg_segmenter.md` に **公開版スコープ (open-core)** 節を追加: 公開版は汎用
  segmentation + FMO/DPD reference 入力生成まで。系別最適パラメータ・命名 bead-typing・
  production チューニングは範囲外 (= 別途 private 拡張)。DPD 既定値 (1.661 / 2.502 /
  stiffness) は cholesterol 由来の汎用既定値である旨を明記

### Added — `abmptools.udfcharge` (new sub-package)

- 単分子 UDF (電荷あり) の per-atom partial charge を抽出し、 バルク系 UDF の
  **同名分子すべて**へ転写する OCTA/COGNAC UDF 電荷割り当てモジュール
  (例: 量子化学 / FMO で求めた電荷を MD バルク系へ反映)
- `read_molecule_charges(udf, mol_name=/mol_index=)` → `MoleculeChargeTemplate`
  (mol_name / n_atoms / charges[e] / atom_type_names / net_charge)。
  `assign_charges_to_bulk(bulk, template, out)` → `AssignResult`
- 電荷規約は gro2udf / udf2gro と共通: `electrostatic_Site[].Type_Name="POINT_CHARGE"`、
  `.ES_Element = charge[e] × 18.224159264`、 `.atom[0] = atom index`
  (`UDFManager.put` の numpy silent-zero を避けるため `float()` cast)
- 割り当ては **atom index 対応**。 atom 数 (必須) と `Atom_Type_Name` 列
  (`verify_atom_types`) を検証してから書き込み、 不一致は strict で例外 / 非 strict で skip。
  出力は別ファイル (入力 bulk 無改変)
- 更新するのは `electrostatic_Site` のみ。 `Structure.Position` (座標) / `Unit_Cell` /
  `bond`・`angle`・`torsion` (結合トポロジー) は **無改変** (座標・トポロジー保持を回帰テストで固定)
- CLI `python -m abmptools.udfcharge --template mol.udf --bulk bulk.udf --out out.udf`
  (`--mol-name` / `--mol-index` / `--no-verify-types` / `--non-strict`)
- 単体テスト `tests/udfcharge/test_udfcharge.py` 11 件 (default_template.udf から
  methanol 系をプログラム生成、 座標保持 + トポロジー保持テスト含む)。 サンプル
  `sample/udfcharge/` は **topology + 電荷 + 実座標 + bond/angle/dihedral + Unit_Cell**
  を持つ OCTA viewer で分子の形ごと開ける UDF (methanol 単分子 + 2×2×2 格子 bulk の
  end-to-end)、 docs `docs/udfcharge.md` + `docs/tutorial_udfcharge.md`

### Added — `abmptools.formulation` Phase 2: Windows native OpenFF route

- `FormulationBuildConfig.force_field_route` で力場経路を 2 つから選択:
  - `"amber"` (default) — tleap + acpype + parmed。**Linux / macOS のみ** (AmberTools)
  - `"openff"` — **全 OS (Windows native)**。PDBFixer + OpenFF `Topology.from_pdb`
    + ff14SB SMIRNOFF で tleap / acpype を一切使わず peptide formulation を build
- OpenFF route の実装 (`topology_openff.py` / `peptide_atomistic_openff.py` /
  `builder._build_openff`):
  - **PDBFixer** で water 除去 + 欠損補完 + explicit H 付加 (OpenFF 必須)
  - **`Topology.from_pdb`** (≠ `Molecule.from_polymer_pdb`) で multi-chain protein を
    1 分子認識 + **disulfide 自動検出** (insulin 2G4M の S-S を手動宣言なしで処理)
  - protein = ff14SB library charges、small molecule = precomputed (gasteiger/am1bcc)、
    ion = ff tip3p charges を `protein_flags` で振り分け
  - **O(N²) 回避**: full mixture の `Interchange.from_smirnoff` は nonbonded exception
    生成が O(N²) で爆発する (insulin×6 で CPU 2h+) ため、**単一コピーを parametrize
    (~4 s) → `[molecules]` の count を実値に書き換え**て複製。`.gro` は packmol mixture
    から `gmx editconf`、溶媒和は `gmx solvate` + `gmx genion` (Joung-Cheatham)
  - sequence 入力は PeptideBuilder + Biopython で 3D extended chain 生成 (natural L-AA)
- amber route の multi-chain / disulfide 対応 (`topology._patch_mixture_pdb_with_ter_records`
  + `builder._disulfide_tleap_bonds`): packmol が drop する TER record を chain ID 遷移 +
  resnum reset 検出で再挿入、disulfide Cys を CYS→CYX rename + global resid map で
  `bond sys.<i>.SG sys.<j>.SG` 注入
- `ndx.py`: >99999 atom 系 (insulin 156k) で gro の atom-serial 列が 100000 で wrap し
  nvt grompp が "Invalid atom number 0" になる bug を 1-based 連番に修正 + 回帰テスト
- 新規 extras `[formulation-openff]` (openff-toolkit / openff-interchange /
  openff-amber-ff-ports / openmm / pdbfixer / PeptideBuilder / biopython)
- 新規 sample: `octreotide_{l,d}_aggregation_10mM/`、`hexarelin_l_aggregation_10mM/`、
  `insulin_aggregation_10mM/` (論文 Hossain 2023 準拠 10 mM = peptide × 6)
- 実機検証: octreotide L/D 10mM (81k/90k atoms、~10h GPU)、insulin 10mM
  (156k atoms、15.2h / 143 ns/day on RTX 4070 Ti) を 100 ns 完走。insulin は
  tetramer→hexamer の irreversible aggregation、TYR14(A)/PHE22(B-Phe1) が凝集 hotspot
- docs: `platform_support.md` (OS 別対応表 + OpenFF route 実装点)、`formulation.md` /
  `dependencies.md` に OpenFF route 節

### Added — `abmptools.formulation` Phase 2-D: Windows native CI

- `.github/workflows/windows-native.yml` (abmptools 初の CI workflow)。`windows-latest`
  runner で OpenFF route と `abmptools.trajectory` の **Windows native 動作を実機検証**:
  - `windows-pure-python` — pip のみで trajectory + formulation models/ndx/mdp unit test
    (pathlib / subprocess shell=False の OS 依存箇所を検証)
  - `windows-openff-smoke` — micromamba (conda-forge) で OpenFF stack →
    `tests/integration/test_openff_windows_smoke.py` (sequence → PeptideBuilder →
    PDBFixer → `Topology.from_pdb` → ff14SB Interchange の parametrization smoke、@slow)
  - trigger は `workflow_dispatch` (手動) + formulation / trajectory 変更の PR。
    packmol / gmx は外部ツールなので CI では呼ばず AmberTools 非依存の Python 層だけ検証
  - **windows-latest runner で実機検証 green (2026-06-14)**: pure-python 52 unit +
    openff-smoke 2 slow が PASS。 検証で判明 — **Windows conda では `openff-toolkit-base`
    必須** (メタパッケージ `openff-toolkit` は AmberTools に hard-depend、 Windows conda
    build が無く solve 不能)。 base は RDKit backend で ambertools 非依存、 ff14SB
    library charges のみの本 route には十分。 docs/platform_support.md + pyproject の
    `[formulation-openff]` に注記反映
  - **小分子電荷を NAGL で Windows native 化 (2026-06-16 green)**: protein は library
    charges で電荷計算不要だが、 小分子 (caprate 等) の AM1-BCC は `sqm` (AmberTools、
    Windows 無し) が要る。 → **`openff-nagl` (ML graph neural net、 pure-Python、 sqm 不要、
    total charge 保存) で代替**。 `molecule_prep` は `NAGLToolkitWrapper` を明示渡しに
    robust 化 (global registry 非自動登録対策)、 smoke に caprate anion → NAGL → Sage
    Interchange を追加 (3 slow PASS)、 CI env + `[formulation-openff]` extra に
    `openff-nagl` + `openff-nagl-models` 追加

### Added — `abmptools.trajectory` (new sub-package)

- cross-platform GROMACS trajectory post-processor (`thin_nojump` / `nojump` / `thin` /
  `wrap_pbc` / `energy`)。subprocess shell=False + pathlib + shutil.which で **Windows
  native** 動作。旧 bash script (wrap_pbc.sh 等) を置換。`python -m abmptools.trajectory`
  CLI + `--gmx` で gmx version 明示 (GROMACS 2026 tpr は古い gmx で読めない gotcha 対応)

### Added — `cg_segmenter` から FCEWS `segment_data.dat` 生成

- 新規 module `abmptools/fragmenter/cg_segmenter/fcews_export.py`。cg_segmenter の
  segment 分割 (ring 検出 + chain MW walk、cholesterol 等) を流用し、**FCEWS** が
  `abmptools.abinit_io.config_read` 経由で読む `segment_data.dat` (FMO フラグメント
  定義、`mode='FMO'`) + monomer `.xyz` を生成する:
  - `build_fcews_segment_data()` — segment partition → FCEWS form entry。atom 番号は
    monomer 内 local 1-origin、H atom は親 heavy atom の fragment に割当
  - `write_fcews_segment_data()` / `write_molecule_xyz()` / `export_fcews()`
  - `CGSegmenter.export_fcews()` メソッド + CLI `python -m
    abmptools.fragmenter.cg_segmenter fcews`
- **connect の atom 順序は `[BDA, BAA]`** (ABINIT-MP の AJF fragment 接続情報の並び)。
  **`connect_num` は BAA 数** = 各 cut bond の BAA を含む fragment に加算
  (LOGManager の `fbaas` と同じ。`connect` の第2要素 BAA が `connect_num>0` の
  fragment 内)。FCEWS 手書き test-data (nafion) に一致。BDA/BAA の役割は
  `auto_split.decide_bda_baa_for_manual_cut` を流用
- FMO フラグメントは atom を共有できないため `allow_atom_sharing=False` (partition)
  必須。共有検出時は `ValueError`。CLI `fcews` は自動で強制
- 単体テスト `tests/cg_segmenter/test_fcews_export.py` (7 件)。`config_read`
  in-process round-trip / connect `[BDA,BAA]` 順序 / cholesterol (C27H46O, 74 atoms)
  の完全 partition + round-trip / 共有検出エラー / xyz↔seg_info 整合 / CLI。
  cg_segmenter 全体で 36 passed
- 第一弾は `mode='FMO'` のみ (rigid cluster `solv` / oligomer `term`・`repeat` は対象外)
- **注**: `abmptools.fragmenter` 本体 (log2config 経路、pdb2fmo 用) は変更なし

### Added (docs)

- `docs/hbond.md` §7 と `sample/hbond/imc_amorphous/README.md` に距離 / 角度分布
  (`distance_dist`) の出力例画像 (IMC amorphous, `imc_hbond_distplot_distance_*`)
  を埋め込み。サンプル参照出力 (CSV / PNG) 自体は別コミットで追加済み。

## [2.1.0] - 2026-06-08

### Added — `abmptools.hbond` 距離 / 角度分布 (`distance_dist.py`)

- 検出済み H-bond の `d(D...A)` / `∠(D-H...A)` を全 record 横断で集計し、
  3 視点のプロット + 統計 CSV を 1 ランで生成:
  - A: `<prefix>_distance_hist.png` — 全 H-bond の `d_DA` 1-D ヒストグラム
    (mean / peak 縦線注釈付き)
  - B: `<prefix>_distance_by_class.png` — クラス別重ね描き (step + 半透明
    fill)。imc mode は `COOH-COOH (dual)` / `COOH-COOH (chain/single)` /
    `COOH-amide` の 3 群、generic mode では `(donor_type, acceptor_type)`
    pair 別。dual 識別は `FunctionalGroupClassification.carboxyl_roles[]
    .dual_partners` を参照
  - C: `<prefix>_distance_angle_2d.png` — `(d_DA, ∠(D-H...A))` 2-D heatmap
    (default 0.05 Å × 5°)
  - `<prefix>_distance_stats.csv` — クラス別 `n / mean / median / std /
    peak / p25 / p75` 表
  - `<prefix>_distance_hist.csv` — long-form (`label, bin_center_DA, count`)
    で再描画用
- 新規 module `abmptools/hbond/distance_dist.py`。`__init__.py` に export 追加
- `AnalyzerConfig` に `do_distance_plots=True` / `distance_d_min=2.0` /
  `distance_d_max=3.6` / `distance_bin_width=0.05` / `angle_bin_width=5.0`
  フィールド追加。CLI に対応する `--no-distance-plots` / `--distance-d-min`
  / `--distance-d-max` / `--distance-bin-width` / `--angle-bin-width` option
- 単体テスト `tests/hbond/test_distance_dist.py` (8 件) 追加。dual ペア
  識別、空入力時の `None` 返却、generic mode pair-type 分け、CSV
  round-trip を回帰固定
- IMC sample (`SI/IMC_result450.0_out_rec900.bdf`、1 record) 実機検証:
  全 81 H-bonds, mean d_DA=2.82 Å, peak 2.72-2.78 Å。`COOH-amide` (N=50)
  が最も短い側に立ち上がる (peak 2.72 Å)、`COOH-COOH (dual)` (N=10) は
  peak 2.82 Å — Yuan 2015 IMC NMR の cyclic dimer / chain / COOH-amide
  振り分けと整合
- `docs/hbond.md` 出力一覧 + CLI option 一覧 + 「7. 距離 / 角度分布」セク
  ション追加

### Fixed — `addsolvfrag` が `AutoFrag='ON'` を出力していた不具合

- `addsolvfrag` は「溶質を手動フラグメント分割した雛形 AJF に、スナップショット
  ごとの溶媒フラグメントを追加する」ツールだが、既定で `aobj.autofrag = True` と
  なっていたため、出力 AJF が `AutoFrag='ON'` かつ `&FRAGMENT` ブロック空という、
  本来の目的と矛盾した内容になっていた（`-ma/--manual` を明示しないと正しい出力に
  ならなかった）。
- 既定を `aobj.autofrag = False` に修正。`-ma` 無しでも出力は `AutoFrag='OFF'` ＋
  雛形フラグメント表に溶媒（HOH/WAT/NA 等）フラグメントを連結した完全な
  `&FRAGMENT` テーブルとなる。`6lu7-covhip` サンプルで NF=311（溶質）→ 1793
  （+1479 水 +3 Na）を確認。
- `-ma/--manual` は既定挙動と同一になったため後方互換目的で残置（no-op）。
- 回帰テスト参照 (`tests/regression/reference/main/addsolvfrag_{covneu,6lu7-covhip}`)
  はバグ挙動を golden master として固定していたため、修正後の正しい出力で再生成。

### Added — `abmptools.trajectory` (new sub-package)

- Cross-platform Python wrapper around `gmx trjconv` / `gmx energy` for
  trajectory post-process (Linux / macOS / Windows、 旧 bash script の置換)
- Public API: `thin_and_nojump`, `nojump`, `thin`, `wrap_pbc`, `gmx_energy`,
  `run_trjconv` (low-level)
- CLI: `python -m abmptools.trajectory {thin_nojump,nojump,thin,wrap_pbc,energy}`
- `aggregation` 系の基本セット (100 ns prod.xtc → `-pbc nojump -skip 10` で
  100 frame、 1 ns stride、 ~300 MB) を Python 1 行で生成
- `subprocess` を `shell=False` + stdin で呼ぶので Windows でも動作、
  `pathlib.Path` で path 区切り差を吸収
- 18 unit tests (`tests/test_trajectory_postprocess.py`)、 gmx subprocess を
  monkeypatch して引数組立 / output 命名 / center group の 2 行 stdin / ndx flag
  / energy term 列挙 / FileNotFoundError / GmxError を検証

### Added (`abmptools.amorphous` cluster center + posres、 branch `fix/cluster-cut-minimum-image-wrap`)

- `BuildConfig.cluster_pdb_path` (str) — pre-built rigid cluster PDB (e.g.
  water trimer の H-bond triangle)。 packmol input に `fixed <center> 0 0 0`
  constraint で box 中央に rigid block 配置、 first component の `n_mol` を
  cluster 分減算 (UDF route の `cluster_file` 等価)。
- `BuildConfig.frozen_atom_indices` (List[int]) — 1-based global GROMACS ndx
  atom indices。 set すると:
  - `system.ndx` に `[ FrozenAtoms ]` index group 追加
  - first moleculetype に `[ position_restraints ]` (`#ifdef POSRES_TRIMER`
    ガード) 追加、 homo pair で件数が cluster mol 数より多い場合は
    `<name>_TRIMER` (cluster 分) + `<name>` (残り) に moletype split。 split
    時 GROMACS の "Only one moltype with [settles] allowed" 制約回避のため
    TRIMER 側 `[ settles ]` を等価な `[ constraints ]` (LINCS-based) に変換
  - 02_nvt_highT 以降の全 mdp に `define = -DPOSRES_TRIMER` (EM は skip、
    packmol overlap relax での numerical 不安定回避)
- `BuildConfig.posres_force_constant` (float, default 10000 kJ/mol/nm²) —
  harmonic restraint の force constant。
- `udf_io.getcontactstructure` の `cutmode='around'` mode — solute (center
  cluster) atoms から `contact_criteria` Å 以内の atom を持つ mol を neighbor
  とする (hydration shell 切り出し)。 既存 legacy mode (`cutmode='contact'`、
  COM 距離 < 2 (r_i + r_j)) と切替可。
- `udf_io.getcontactstructure` の `contact_criteria` default 4.0 Å (1st
  hydration shell)。
- 中央 cluster (mixflag/clusterflag mode) の周辺 27-cell image-neighbor を
  minimum-image wrap する fix — 修正前は中央 cluster の周りに「+L/-L 方向に
  飛び離れた phantom cluster」が現れ FMO 非収束の主因だった。

### Fixed (`fix/cluster-cut-minimum-image-wrap` branch)

- `trajectory_ingest._settles_per_moltype` — TIP3P/TIP4P 水で `[ bonds ]`
  が空 (settles で剛体拘束) の場合、 MDAnalysis Universe に O-H bond が登録
  されず `make_whole` が境界跨ぎ water を unwrap できない問題。 `[ settles ]`
  block から OW を head に O-H bond を synthesize。
- `molcalc.Exportardpos` を obabel 経由から direct Python PDB writer に置換
  — obabel の XYZ→PDB 自動 bond perception が境界跨ぎ broken water の H を
  orphan (`ATOM/UNK` 行) として誤判定する問題を解消。 nameAtom (atom_type) +
  molnames (Mol_Name) を信頼して PDB strict column layout で write。
- `udf_io.moveintocell` の `np.float32` silent-zero bug を回避するための
  cell float cast cherry-pick (`51cd7b3`、 旧 `fix/trajectory-ingest-float-cast`
  branch)。

### Changed

- `sample/formulation/_postprocess/trajectory_thin_nojump.sh` を DEPRECATED 化
  (機能は `abmptools.trajectory` に移行)
- `sample/formulation/octreotide_{l,d}_aggregation_100ns/README.md` に
  Post-process 節を追加 (Python CLI + API の呼び出し例)
- **`amorphous.mdp_protocol.write_wrap_script`** が生成する script を
  `wrap_pbc.sh` (bash) → **`wrap_pbc.py` (Python)** に変更。
  生成 file は `from abmptools.trajectory import wrap_pbc` を呼び、
  Linux / macOS / Windows のどこでも `python wrap_pbc.py` で実行可能。
- **`amorphous.mdp_protocol.write_udf_export_script`** が生成する script を
  `gen_for_udf.sh` (bash) → **`gen_for_udf.py` (Python)** に変更。
  生成 file は `abmptools.trajectory.gmx_energy` + `nojump` を呼ぶ。
- 全 `sample/amorphous/*/README.md` + `run_sample.sh` + `docs/*.md` の
  `bash wrap_pbc.sh` / `bash gen_for_udf.sh` 表記を
  `python wrap_pbc.py` / `python gen_for_udf.py` に統一。
- `getcontactstructure` の `freezegrps` → `position_restraints` (harmonic)
  に switch — freezegrps は LINCS/SETTLE + Domain Decomposition で
  `determinant = -inf` で abort するため。 ~0.04 nm の wiggle のみで cluster
  geometry を維持する (branch `fix/cluster-cut-minimum-image-wrap`)。
- `mdp_protocol.write_all_mdp` の `freeze_group` 引数を `define_posres`
  に rename (同 branch)。


### Added — `formulation.analysis` workflow (Hossain 2023 Fig 1b/1c/2/4 1-コマンド再現)

- `run_analysis(traj, top, out_dir, n_peptides, enhancer_resnames, bile_salt_resnames)` で
  集合体形成 (Fig 1b/1c/Fig 2) + per-residue contacts (Fig 4) + plot を一括実行
- 新規低レベル API: `per_frame_clusters_heavy_atom` (PBC-aware heavy-atom min distance)、
  `compute_per_residue_contacts` (cap 除外 + per-peptide 平均化)
- `aggregate_transition.compute_aggregate_transitions` を atom-index 分割 + heavy-atom
  cutoff + PBC 補正に refactor (GROMACS の resid 1-N リセット慣習に対応)
- CLI: `python -m abmptools.formulation analyze --traj ... --top ... --out ... --n-peptides 6 --enhancer-resnames CPN,CPC --bile-salt-resnames TCH`
- plots: `plot_aggregate_timeseries` / `plot_max_size_distribution` / `plot_per_residue_contacts`
- 出力: `aggregate/{cluster_states,aggregate_size_timeseries}.csv` + `aggregate_summary.json`、
  `contacts/per_residue_contacts.{npy,csv,json}`、 `plots/*.png`
- L 体 10mM 100 ns で smoke 完了: max size mode=6 (49%)、 pct_agg mean=92.8%、
  Phe3/Trp4 が enhancer 主要 contact (10.0/9.8) で論文と整合

### Added (docs)

- **`docs/platform_support.md`** — OS 別 (Linux / macOS / Windows native /
  WSL2) で各 sub-package + 依存 module がどこまで動くかの早見表 + 用途別
  推奨 setup シナリオ + Phase 1 (`formulation` Windows route) 計画。
- `docs/dependencies.md` と `docs/formulation.md` の冒頭に
  `platform_support.md` への誘導 1 行を追加。

### Added (sample)

- `sample/formulation/octreotide_l_aggregation_10mM/` — Hossain 2023 主要系
  (10 mM peptide、 Fig 1-2 出典) 準拠で `n_copies = 6` の設定。 既存
  `octreotide_l_aggregation_100ns/` (peptide × 2 = 3.3 mM) は論文 Table 1 の
  どの系とも一致しない hybrid 構成だったことが判明し、 訂正用として追加。

### Added — `formulation` OpenFF route (Phase 1 完了)

- `force_field_route: "openff"` config 切替で **AmberTools 依存を完全回避**:
  - peptide: `Molecule.from_polymer_pdb` (要 `pdb_path` 入力)
  - small molecule: SMILES → OpenFF Sage 2.x SMIRNOFF
  - charge: peptide=Gasteiger / small mol=AM1-BCC default、 NAGL 切替可
  - typing: amorphous の `create_interchange` を再利用 (Apache-2.0 互換 OpenFF
    + Interchange 経路)
  - GROMACS export: `Interchange.to_gro()` + `.to_top()`
- 新規 module 3 本 (本実装):
  - `formulation/peptide_atomistic_openff.py` (~135 行)
  - `formulation/small_molecule_openff.py` (~105 行)
  - `formulation/topology_openff.py` (~140 行)
- `builder.py` に `_build_openff` ルート (~110 行) 追加、 既存 amber route と
  独立並行、 stage 1/2/4 のみ差替、 stage 3/5/6/7 (packmol / ndx / mdp /
  run_script) は共有
- 実機検証: kggggg × 2 + caprate × 16 + TC × 1 (box 6 nm、 water なし) で
  end-to-end build PASS (707 atoms)、 grompp em PASS
- Amber route 28 unit tests 回帰なし

### TODO (Phase 2)

- **water solvate + ions balance を OpenFF route 内に追加**: 現状 Phase 1 の
  build 出力は dry-mixed system (water なし)。 OpenMM `Modeller.addSolvent`
  経由で TIP3P 充填 + Joung-Cheatham Na/Cl + 0.15 M NaCl を Interchange に
  combine する経路を追加予定。
- **sequence からの peptide build** (PDBFixer 経由): 現状 OpenFF route は
  `PeptideSpec.pdb_path` 必須。 sequence 単独入力からの 3D build 経路を追加。
- **ff14SB SMIRNOFF for natural L-AA peptide**: 現状 whole-peptide Sage 経路
  のみ。 `openff-amber-ff-ports` で peptide 部分だけ ff14SB SMIRNOFF を適用
  する FF stack 経路を追加。
- **Windows native 環境 (実機 Win10/11) での実走確認 + CI 化**。

## [2.0.0] - 2026-05-28

メジャーリリース。v1.15.4 (2026-04-19) 以降に develop で積み上げた 176 commits 分の機能・修正をまとめて公開。

### Breaking changes

- **gro2udf の UDF atom field semantics 変更** (commit `6b32e92`):
  - `Set_of_Molecules.molecule[].atom[].Atom_Name` = **element symbol** (`C` / `H` / `O` 等、旧版は GROMACS atom name `ca1` / `ha0` 等)
  - `Set_of_Molecules.molecule[].atom[].Atom_Type_Name` = system.top atomtype 列、OpenFF SMIRNOFF の `MOL0_<N>` は `<element><N>` (例: `C4` / `H19`) に rewrite
  - `Set_of_Molecules.molecule[].atom[].Atom_ID` = per-molecule **local 0-indexed** (旧 global counter)
  - 1.x との UDF 出力差は test_regression で確認可能
- **abmptools の license 表記更新**: v1.23 以降は Apache-2.0 + NOTICE + CITATION.cff (≤ v1.22.0 は MIT 固定)

### Added (`abmptools.gro2udf` Time/Pressure/Density/Temperature 同期 — v1.x.x 候補)

- `MdpParams` に新規 accessor 追加: `dt`, `nsteps`, `nstxout_compressed`
  (nstxout-compressed / nstxout_compressed / nstxout / nstvout fallback chain),
  `nstenergy`
- `TopModel` に新 field 追加: `dt_ps`, `nsteps`, `nstxout_compressed`,
  `nstenergy` (defaults: 0.001 / 0 / 0 / 0)
- `TopAdapter` で mdp 値を上記 field に格納
- `top_exporter._set_default_condition()` で
  `Simulation_Conditions.Dynamics_Conditions.Time` の以下を mdp から書込み:
    - `delta_T` ([ps] -> [tau] 変換)
    - `Total_Steps`
    - `Output_Interval_Steps` (nstxout-compressed 優先、fallback nstenergy)
- xvg → UDF Statistics_Data mapping を拡張 (旧 `_XVG_TO_UDF_ENERGY` を
  `_XVG_TO_UDF_STATS` に rename + path/unit tuple 形式に):
  - Total Energy / Hamiltonian は既存 mapping にあったが、複数 LJ / Coulomb
    の合算で潰されていた → "Total Energy" -> Energy.Instantaneous.Total
  - 新規追加: Temperature -> Temperature.Instantaneous ([K])
  - 新規追加: Pressure -> Pressure.Instantaneous ([bar])
  - 新規追加: Density -> Density.Instantaneous ([kg/m^3])
  - 新規追加: Volume -> Volume.Instantaneous ([nm^3])
- `_aggregate_energy_per_frame` を `_aggregate_statistics_per_frame` に rename
- `_append_structure` の energy_values 形式を `{(path, unit): value}` に変更
  (旧 `{field: value}`)、Statistics_Data.{Energy, Temperature, Pressure,
  Density, Volume} 全部に対応
- tests/test_mdp_parser.py に 8 件追加 (dt/nsteps/nstxout_compressed/
  nstenergy の default + parsed + fallback chain)、16/16 PASS
- 実機検証 (ketoprofen amorphous):
    delta_T=0.0205 [tau] (= 0.001 ps)
    Total_Steps=500000
    Output_Interval_Steps=5000
    各 frame で Bond/Total/Temperature/Pressure/Density/Volume が record に embed

### Added (`abmptools.gro2udf` multi-frame trajectory + xvg energy — v1.x.x 候補)

- 新規 `abmptools/gro2udf/trajectory_io.py`:
  - `frames_from_multi_gro(path)` — `gmx trjconv -pbc nojump -o output.gro`
    形式の multi-frame `.gro` を pure-Python で parse (stdlib only)
  - `read_xvg(path)` — gmx energy 出力の xvg を `(times, {legend: values})`
    に parse
- `TopExporter.export()` / `.export_model()` に
  `trajectory_path` / `energy_path` / `energy_times` / `energy_series`
  パラメータを追加。CLI で `--trajectory <.gro/.xtc>` / `--energy <.xvg>`
  を指定すると全 frame の Structure record + Statistics_Data.Energy が
  1 UDF に embed される
- xvg → UDF energy field map (`_XVG_TO_UDF_ENERGY` dict):
  Bond / Angle / Proper Dih.+Improper Dih.→Torsion / LJ-14+LJ(SR)+
  Disper.corr.→Nonbonding / Coulomb*→Electrostatic / Potential / Kinetic /
  Total
- frame time と xvg time grid は nearest-neighbour で照合 (501 row vs
  101 frame など denser xvg からも正しく取れる)
- schema 上に energy field が無い cognac 古版では silently skip (try/except)
- docs/gro2udf.md にトラブルシューティングの前 section
  「Multi-frame trajectory + energy を 1 UDF に embed」を新規追加
- tests: `_capture` helper signature を `energy_values=None` 受け入れに更新
  + 既存 6/6 PASS
- 実機検証 (ketoprofen amorphous 50 mol × 33 atom × 101 frame):
    totalRecord=101、各 frame で Cell.a + Bond + Angle 値が正常に embedding

### Added (`abmptools.gro2udf` OCTA8.4 対応 + 単位 fallback — v1.x.x 候補)

- **`abmptools/gro2udf/default_template_cognac101.udf` (新規 bundle)** —
  OCTA8.4 / OCTA8.4 (cognac10.1 までしか同梱されていない環境)
  向けの schema 互換 minimal template。`\include{"cognac101.udf"}` +
  `Unit_Parameter:{"","",1.0,4.184,0.1}` を含み、cognac10.1 でも
  `[nm]` / `[ps]` / `[kJ/mol]` の unit alias が解決される。cognac11.2 default
  template と **出力データは完全同一** (Cell.a, Unit_Parameter.Length 等)。
- `top_exporter._put_with_unit_fallback()` helper を追加 — unit 引数つきの
  put が失敗した場合に unit なしで retry する fallback。`Unit_Parameter` が
  欠けた template でも動くようにする保険。
- `top_exporter._rewrite_cognac_include()` helper + `--cognac-version` CLI
  option を追加。
- CLI: `--cognac-version 101` / `102` 指定時に bundled の cognac10.x template
  を自動選択 (`--template` を渡さなくても OK)。
- `pyproject.toml` の `package-data` を `gro2udf/default_template*.udf` に
  glob 化、両 template を install パッケージに同梱。
- `UDFExportError` の hint メッセージを cognac10.1 環境向けの具体ガイドに
  刷新 (推奨 `--cognac-version 101`、代替 GOURMET Save As + `--template`)。
- docs/gro2udf.md にトラブルシューティング section
  「OCTA8.4 / OCTA8.4 で `file not found:cognac112.udf` エラー」
  を新規追加、`--cognac-version 101` を推奨対処として記載。
- bundled cognac10.1 template の元 (`A20B40A20_in.udf`) はユーザー提供
  (粗視化 sample、`Unit_Parameter` を追加して unit alias を有効化済み)。

### Added (`abmptools.gro2udf` 診断付きエラー — v1.x.x 候補)

- `top_exporter.UDFExportError` 例外クラス + `_section()` context manager を
  追加。`TopExporter.export_model()` の各 stage (`template-copy`,
  `UDFManager-open`, `erase-existing-records`, `Set_of_Molecules`,
  `Structure[record=N]`, `default_condition`, `Molecular_Attributes`,
  `Interactions`) を context manager で wrap して、UDFManager の cryptic な
  例外 (RuntimeError 等) を **どの section / どの template / どの output / 元の
  underlying error / 対処 hint** を含む `UDFExportError` に再 raise する。
- 主な動機: OCTA84 で gro2udf がエラーになった時、UDFManager がどの field
  で落ちたか分からず原因特定が困難。section 名で範囲を絞り、hint で OCTA
  version 差異を案内 (bundled template の利用 / 当該 OCTA から template を
  再生成等)。
- `import UDFManager` 失敗時は `UDFExportError` で OCTA PATH 設定 hint を表示。
- 存在しない template path 渡し時は明示的に `UDFExportError` を pre-flight
  check で投げる。
- docs/gro2udf.md にトラブルシューティング section 追加。
- tests/test_top_exporter_frames_override.py に 2 件追加 (missing template /
  section failure context)、5/5 PASS。

### Added (`abmptools.amorphous` OCTA viewer export — v1.x.x 候補)

- `mdp_protocol.write_udf_export_script(output_dir, ndx, stage, n_energy_terms)`
  新規 API。OpenFF amorphous protocol 後に OCTA viewer (GOURMET) で読み込み可能な
  energy.xvg + nojump gro を生成する `gen_for_udf.sh` を書き出す。
- 内容:
  - `seq <N> | gmx energy -f <stage>.edr -o <stage>_energy.xvg`
    (default N=50、gmx は存在しない term 番号を silently skip するので
    50 で標準 energy term を網羅)
  - `echo 0 | gmx trjconv -f <stage>.{trr,xtc} -s <stage>.tpr -pbc nojump
    -o <stage>_nojump.gro` (`wrap_pbc.sh` の `-pbc mol -ur compact` と
    違い、分子を box 内に wrap せず PBC を跨いで連続的に追跡。
    OCTA viewer (GOURMET) での軌跡再生に適する)
- `FormulationBuilder.build()` から `wrap_pbc.sh` の隣に自動生成。
  返り値 dict に `udf_script` key を追加。
- 既存 sample (`pva_amorphous/`、`ketoprofen*/`) の `md/` 配下にも script
  配置済み。
- ドキュメント: `docs/amorphous.md` (workflow + tree 図 + 仕様説明)
- テスト: `tests/test_mdp_protocol.py::test_write_udf_export_script_*`
  4 件追加 (default / custom stage + n_energy_terms / omits ndx / executable bit)、
  17/17 PASS。
- 実機検証 (ketoprofen amorphous 50 mol × 33 atom):
  `05_npt_final_energy.xvg` (320 KB) + `05_npt_final_nojump.gro` (7.5 MB,
  101 frames × 1650 atoms) 生成 OK。

### Added — `abmptools.formulation` (v1.30.0 候補)

- AA mixed-solution peptide-formulation builder modeled after Hossain et al.
  2023 (*Nanoscale* 15, 19180-19195) — peptide drug + permeation enhancer +
  intestinal bile salt in a cubic water box. **Commercial-permissive force
  fields only**: AMBER ff14SB + GAFF2 + TIP3P + Joung-Cheatham ions
  (CHARMM36 + CGenFF は学術ライセンスのため不採用).
- `FormulationBuilder.build()` 7-stage orchestrator:
  ff_staging → peptide_atomistic (tleap) → small_mol_parameterize
  (acpype GAFF2/AM1-BCC, lazy RDKit SMILES → 3D) → packmol multi-component
  → solvate_ions (tleap solvatebox + parmed) → index (Python resname table)
  → mdp_render → run_script.
- `python -m abmptools.formulation {example, validate, build, analyze, release_us}`
  argparse CLI; `example` emits a smoke JSON, `release_us` writes pull +
  N-window US MDPs for peptide-from-aggregate PMF.
- Analysis stack (opt-in `[formulation-analysis]` extra): aggregate
  transition matrix + cluster timeseries (MDAnalysis + networkx, lazy
  import), per-residue contact map, `gmx dssp` / `gmx sasa` / `gmx hbond`
  wrappers, matplotlib heatmap + timeseries plots.
- Sample configs in `sample/formulation/`: `kggggg_smoke/` (6 nm box,
  ~8k atoms, 200 ps smoke; runs in ~10 min on 4-core CPU) and
  `insulin_smoke/` (10 nm box, 99k atoms, 1 ns; PDB 2G4M downloaded
  on demand from RCSB, not redistributed).
- License posture: acpype (GPL-3.0) and MDAnalysis (GPL-2.0) are mere
  aggregation — acpype via subprocess only, MDAnalysis via lazy import
  inside `analysis/` modules with `pip install abmptools[formulation-analysis]`
  opt-in install (precedent: `amorphous/trajectory_ingest.py`). The
  abmptools wheel itself remains Apache-2.0 with zero GPL code shipped.
- Tests: 51 (models + mdp + ndx) + 19 (packer + small_molecule + peptide_atomistic)
  + 9 (topology) + 5 (builder) + 11 (analysis) + 5 (umbrella_release)
  = 100 unit tests, all mock-based (no gmx/tleap/acpype required). Slow
  integration smoke is gated by `@pytest.mark.slow + skipif(tool missing)`.

### Added (`abmptools.hbond` element + bond-graph fallback — v1.28.0 候補)

- **`fallback_tag_by_element()` 新規追加** (``func_tags.py``): atom_type が
  ``None`` または per-atom unique な値 (OpenFF SMIRNOFF の ``MOL0_X`` 等) の
  atom について、element + bond graph で機能タグを自動付与:
  - O atom: H と bond → ``hydroxyl_O``、H なし → ``carbonyl_O``
  - H atom: O と bond → ``hydroxyl_H``、N と bond → ``amide_H``
  - N atom: C と bond → ``amide_N`` (amide/amine 区別は後段)
  - C atom: ``carbonyl_O`` (= O without H) と bond → ``carbonyl_C``
- ``AnalyzerConfig.use_element_fallback`` フィールド (default True) +
  CLI ``--no-element-fallback`` (strict mode への切替) 追加
- ``functional_groups._tag_atoms_of_mol`` が None tag のみ fallback で補完
  (FF mapping が hit した atom は touch せず、mapping 優先)
- 効果: **OpenFF Sage 経由の UDF を antechamber 経由の GAFF type patch なし
  で直接解析できる**。PVA amorphous sample の手順から antechamber step が
  消える (build_bdf.py 実行 + hbond CLI 直叩きで完結)
- 実機検証: PVA 10-mer × 30 で element fallback ON / antechamber patch なし
  と OFF / antechamber patch ありで完全同一の結果 (rec=0 で 188 H-bonds、
  ratio_donor_busy=62%、ratio_acceptor_busy=60%)
- IMC は GAFF2 type が既に正しく付いているので fallback は no-op (既存
  baseline dual=10 / chain=41 / single=38 / free=36 を維持、57 → 64 tests
  PASS)
- ``tests/hbond/test_element_fallback.py`` 7 件追加:
  alcohol_OH / carboxyl pattern / amide NH / preserves_existing /
  detect_carboxyls_unknown / detect_hydroxyls_unknown / disabled_no_groups

### Added (`abmptools.hbond` generic mode — v1.28.0 候補)

- **`classify_mode={imc, generic}` を新規追加** (`AnalyzerConfig.classify_mode`、
  CLI `--classify-mode`、Jupyter UI dropdown)。default = imc で既存挙動維持
- **`imc` mode** (既存): COOH 中心 4-species (dual/chain/single/free)
- **`generic` mode** (新規): donor-type × acceptor-type の pair 統計 + atom
  単位の role tag (Donor/Acceptor/Both/Candidate)。PVA / peptide / アルコール /
  混合系等、COOH を持たない任意系で動作
- 新規 module `abmptools/hbond/pair_type_stats.py`:
  `PairTypeStat` + `GenericPairClassification` dataclass、`classify_generic()` /
  `summarize_pair_stats()`
- colorizer に generic 版 3 関数: `write_hbond_attributes_generic` /
  `colorize_udf_action_generic` / `write_show_python_script_generic`
- `DEFAULT_GENERIC_COLORS`: Donor=red / Acceptor=cyan / Both=magenta /
  Candidate=faint gray
- 新規出力: `<prefix>_pair_stats.csv` (generic mode のみ)、
  pairs.csv の `kind` 列は generic では `<donor_type>-><acceptor_type>` 形式、
  Attributes 値は `Donor` / `Acceptor` / `Both` (Candidate は skip)
- tests: `test_pair_type_stats.py` 5 件追加 (no_hbonds / single_pair / both_role /
  multiple_pair / unique_dedup)、計 57/57 PASS
- IMC 系での generic mode 動作確認 (`--donor-groups carboxyl
  --acceptor-groups carboxyl_O,amide_O`): carboxyl→amide_O=50, carboxyl→carboxyl_O=31
  が pair_stats に出る (= imc mode の hb_cc=31 + hb_ca=50 と一致)

### Changed (`abmptools.hbond` 4-species classifier + NMR 比較 plot — v1.27.0 候補)

- **per-COOH 分類を 4 species に拡張**: dual / **chain** / single / free。
  Yuan et al. (2015) *Mol. Pharm.* 12, 4518 (DOI 10.1021/acs.molpharmaceut.5b00705)
  の amorphous IMC 13C SSNMR deconvolution (cyclic dimer ~179 / disordered chain
  ~176 / COOH-amide ~172 / free ~170 ppm) に対応する分類軸。
  - **chain** = cyclic dimer ではない COOH-COOH 一方向 H-bond の参加者 (donor or
    acceptor 側)。論文の "disordered chains having various lengths" + chain end
    + ring-larger-than-dimer を一括して捕捉
  - 優先度: dual > chain > single > free (per-COOH も mol 代表 role も同じ)
  - `CarboxylRole.chain_partners` / `FunctionalGroupClassification.n_carboxyls_chain`
    + `ratio_carboxyl_chain` + `MolRole.n_carboxyls_chain` を新規追加
  - summary.csv に `n_carb_chain` / `ratio_carb_chain` / `n_chain_mols` カラム追加
  - colorizer: Mol_Name リネーム経路に `_CHAIN` (Magenta)、action/script 経路に
    `chain: [0.85, 0, 0.85, 1]` (magenta) を追加、write_hbond_attributes の
    value_map にも `chain: Chain` を追加
- **IMC ベースライン値の再変更**: 旧 (3-species) single=49 / free=66 → 新 (4-species)
  chain=41 / single=38 / free=36 (合計 125 = 旧 single+free 115 のうち chain に
  41 が分流)。dual=10、amide accept=49 / free=76 は不変。
- **`sample/hbond/imc_amorphous/plot_nmr_comparison.py` 新規追加**: 3-row
  plot (Yuan Figure 5 image / Yuan Table 1 deconv bars / MD per-COOH bars) を
  同一 ppm 軸で並列描画。生成 PNG `output/imc_hbond_nmr_comparison.png` を
  同梱。比較で大きいギャップ (NMR dual 58.5% vs MD 8.0% 等) が視覚化される

### Changed (`abmptools.hbond` per-functional-group 統計 — v1.27.0 候補)

- **分類モデルを per-functional-group に変更**: v1.26 まで「分子単位 1 役割」
  だった分類を「官能基単位 (COOH ごと / amide ごと)」に置き換え。1 分子内に
  複数の COOH や amide がある場合に役割が混在するケースに正しく対応する。
  - 新 dataclass: `CarboxylRole`(role=dual/single/free) + `AmideRole`(role=accept/free)
  - 旧 `MolRole` には per-mol の `n_carboxyls_dual/single/free` + `n_amides_accept/free`
    フィールドを追加。`role` 属性は **分子代表 role** (色付け用、優先度
    dual > single > free)
  - `ClassificationResult` は `FunctionalGroupClassification` の alias
    (backward compat、外部 import パス維持)
- **summary.csv 拡張**: per-functional-group カラム (`n_carb_dual`,
  `n_carb_single`, `n_carb_free`, `n_amides`, `n_amide_accept`,
  `n_amide_free`, `ratio_carb_*`, `ratio_amide_*`) を追加。
  従来の `n_dual_mols/n_single_mols/n_free_mols` (mol-level representative)
  も末尾に残す。
- **`<prefix>_classification.csv` 新規追加**: 全 carboxyl / amide ごとの
  (record, group_type, mol_index, group_index, role, partner_count, partners)
  テーブル。NMR の COOH-C / amide-C 信号分離と直接対応する物理量。
- **`<prefix>.bdf` (OCTA viewer プリ描画用コピー) 新規出力**: Mol_Name 維持 (元の
  `molecular` 等) の単純コピー。`<prefix>_colored.bdf` (Mol_Name リネーム済) は
  OCTA viewer のプリ描画で空表示になる問題への対処。CLI option `--no-copy-uncolored`
  でスキップ可。
- **IMC ベースライン値の変更**: 旧 single=73 (COOH→amide H-bond の両当事者を
  カウント) → 新 single=49 (COOH 状態が single の COOH のみカウント、amide
  acceptor 側 mol は COOH free なら free に入る)。`test_imc_regression.py`
  の許容値を更新。物理的には新値が NMR の COOH 信号と直接対応する。
- **verbose log の表現変更**: `COOH dual/single/free=10/49/66 (8%/39%/53%),
  amide accept/free=49/76 (39%/61%)` のように官能基単位の比率を表示。
- **gourmet 可視化手順の docs 改訂**: Python パネルで
  `show.all("line","mol","molname",...)` に書換が必要なこと、OCTA viewer プリ描画
  には `<prefix>.bdf` を使うことを明記。
- **Python action 経路の追加** (`colorize_mode="action"`): GOURMET
  `Draw_Attributes` schema が Mol_Name 維持での per-functional-group
  色付けに対応していないため、`<prefix>_show.act` (`autorun: showHbond()`)
  + `<prefix>_action.bdf` (Mol_Name 維持コピー + Action ヘッダパッチ) を
  併出する経路を追加。各 carboxyl atoms (c/o/oh/ho) と amide atoms (c/o/n)
  を role に応じた色 (dual=red, single=blue, free=gray, accept=cyan) で
  sphere overlay 描画。1 分子内に複数官能基が異なる役割で参加するケースも
  正しく可視化される。CLI `--colorize-mode {molname,action,both}`、default は
  backward compat の `molname`。`colorize_udf_action()` API 新規 export。
- **`<prefix>.bdf` の Attributes[] に hbond タグを append** (OCTA viewer Attribute
  フィルタ対応): 各 functional-group atom (carboxyl c/o/oh/ho、amide c/o/n)
  の `Set_of_Molecules.molecule[].atom[].Attributes[]` 末尾に
  `Name='hbond' Value='Dual'/'Single'/'Free'/'Accept'` を append。既存
  Attributes (`Name='1' Value='molecular:<id>'` 等の OCTA viewer 内部用) は
  touch せず idempotent (再実行で重複 entry 作らない)。OCTA viewer で
  `<prefix>.bdf` を開いた後、Attribute フィルタで `hbond=Dual` 等の atom
  のみ可視化できる (色分けではなくカテゴリフィルタ)。CLI
  `--no-write-attributes` / `--attributes-name NAME` option、default `hbond`。
  `write_hbond_attributes()` API 新規 export、`AnalyzerConfig` に
  `do_write_attributes=True` / `attributes_name='hbond'` フィールド追加。
- **`<prefix>_show.py` 併出** (OCTA viewer (GOURMET) 対応): OCTA viewer (GOURMET) は
  `<prefix>_show.act` の autorun action 形式で落ちることがあるため、同じ描画
  ロジックを autorun ラッパー無しの **プレーン Python script** として
  `<prefix>_show.py` に出力。OCTA viewer で `<prefix>.bdf` (Mol_Name 維持 copy)
  を開いた後、Python パネルから `Load…` → `Run` で同じ役割色 overlay が描画
  される。`write_show_python_script()` API 新規 export。`colorize_mode in
  {action, both}` で .act / .py 両方を併出。

### Added (`abmptools.cg.dpd` サブパッケージ — v1.26.0 候補)

- **DPD 系入力ファイルビルダー**: `cg_segmenter` (R0) で生成した CG segment + fcews
  `aij.dat` (Python 辞書 script) から Cognac DPD 入力 UDF / OCTA viewer dpm を生成する
  `abmptools.cg.dpd` サブパッケージを追加。 dpdgen (Koji Okuwaki 本人作品) のロジックを
  参考に **subprocess も import もせず abmptools 内で自前実装**。
- **3 ルート構成**:
  - R0 (既存): `cg_segmenter dpdgen` で `{name}_monomer` + `{name}_calc_sett` 生成
  - R1 (新): `monomer + calc_sett + aij.dat` → Cognac DPD 入力 UDF (`*_uin.udf`)、
    plain text writer、 UDFManager 非依存 (abmptoolsenv 等 non-OCTA 環境でも動く)
  - R2 (新): `monomer + aij.dat` → OCTA viewer dpm + `monomer-lib/<seg>/Virtual.mom` +
    `#Message.txt`、 **B 案 (user template + abmptools patch)** で権利配慮
- **権利配慮の設計**:
  - **R1**: 冒頭 1 行 `\include{"cognac112.udf"}` で class 定義を OCTA viewer install dir 経由で
    resolve、 abmptools は OCTA / OCTA viewer spec を一切持たない
  - **R2**: user が OCTA viewer で空 dpm template を 1 回作成 → abmptools が `\begin{data}` 内の
    5 ブロック (SegmentModel / SegmentPairModel / PolymerModel / DpdInput / FcewsParam)
    のみを `patch_dpm` で brace-aware に差し替え。 class 定義 (商用 OCTA viewer spec) は
    template のまま温存。
  - **Virtual.mom**: 全 segment dir に user 提供を copy するのみ
- **モジュール構成** (`abmptools/cg/dpd/`, 9 files, ~1435 行):
  - `models.py`: `AijMatrix` / `MonomerSpec` / `CalcSett` / `DpdSystemSpec` データクラス
  - `aij_io.py`: fcews aij.dat read/write + chi → a 自動変換 (Groot-Warren: `a = aii + chi/0.306`)
  - `monomer_io.py` / `calc_sett_io.py`: cg_segmenter 出力との互換 reader
  - `dpm_writer.py`: R2 = B 案 patch (5 block) + `propagate_virtual_mom` + `write_message_txt`
  - `udf_writer.py`: R1 = plain text Cognac DPD UDF writer
  - `orchestrator.py`: `CGDpdBuilder` クラス
  - `__main__.py`: CLI (`build-udf` / `build-dpm`)
- **Cognac PDF → md 変換**: `man/octa/Cognac_jpn_old.pdf` (2.1 MB, 253 ページ) を
  `pdfminer.six` で `Cognac_jpn_old.md` (504 KB) に変換、 ABMPTools repo 外の `man/octa/`
  に保存 (TOC + 254 ページ本文)。
- **fcews aij.dat 統一形式 (確定)**: `aij = [['seg_i','seg_j', value], ...]` または
  `chi = [...]` の Python list of [str, str, float]、 `exec()` で読み込み。
- **cholesterol 検証 (R1)**: UDF 10.5 KB / 615 行、 `\include cognac112.udf` 1 行、
  4 top-level sections、 Bond_Potential 4 / Angle_Potential 3 / Pair_Interaction DPD 15、
  brace+bracket balance OK。
- **cholesterol 検証 (R2)**: dpm 18.6 KB、 def section 温存、 SegmentModel 5 (P0..P4)、
  DpdBond 4 / DpdAngle 3、 `monomer-lib/{P0..P4}/Virtual.mom` 5 個 + `#Message.txt`。
- **多 monomer 混合系対応** (`from_multi_files`、 CLI ``--multi-monomer JSON``):
  cholesterol + water、 lipid + protein 等を 1 builder で扱う。 各 monomer に独立した
  particle 名指定。 segment 名が aij.dat に未定義なら ``logger.warning`` で通知し、
  default ``aii=25.0`` に fallback。 single (`--monomer`) と排他、 どちらか必須。
- **Particle 名マッピング docs**: cg_segmenter 汎用ラベル (P0..Pn-1) と fcews aij の
  segment 名対応を `docs/cg_dpd.md` で命名規約 + 不一致時症状込みで明文化。
- **NOTICE 追記**: DPDgen (本人作品) のロジック移植を「Re-implementation attribution」
  section で明記、 ABMPTools リポジトリ内に DPDgen ソースを含まないことを宣言。
- **Jupyter UI (`open_panel`、 G1)** (`notebook_ui.py`): ipywidgets で
  interactive build panel。 Summary (monomer 一覧 + aij segments + validate
  結果) / R1 UDF build / R2 DPM build (template path 入力付き) / Re-verify
  ボタンを 1 画面に集約。 `from abmptools.cg.dpd import open_panel; open_panel(builder)`。
- **aij skeleton 生成** (`create_empty_aij(segments, aii, mode, off_diagonal)`):
  fcews の自動 aij.dat が使えない場合の hand-craft 出発点。 全 N*(N+1)/2 ペアを
  default 値で生成 (mode='a' なら同種 aii / 異種 off_diagonal、 mode='chi' なら
  同種 0 / 異種 off_diagonal)。
- **Monomer hand-craft** (`build_monomer(name, particle_names, bond12, angle13, ...)`):
  cg_segmenter なしで MonomerSpec を構築する helper。 1-粒子 solvent、 simple
  dimer/trimer 等。 bond/angle potential 値は default 共通 (0.86/50.0 / 余角 0=180°/5.0)。
- **a → chi 逆変換 (G2)** (`AijMatrix.to_chi_values`): Groot-Warren 逆向き
  ``chi = (a - aii) * 0.306``。 ``write_aij`` に ``out_mode='auto'/'a'/'chi'`` 追加で
  出力形式を強制変換可能 (a モード保存の aij.dat を chi で再保存等)。 同種ペア (a==aii) は
  chi=0 になり、 a ↔ chi の round-trip は完全可逆。
- **整合性検証 (`validate` / `verify` CLI)**: 実機 `build-udf` / `build-dpm` 前の
  dry-run。 `CGDpdBuilder.validate()` で warning list を返却 + CLI
  ``python -m abmptools.cg.dpd verify`` (整合 OK = return 0、 missing/extra
  segment あり = return 1)。 監査項目: ① monomer particle 名が aij.dat にあるか
  (missing) ② aij.dat に monomer 未使用 segment が無いか (extra)。
- **テスト**: 59 件 (`tests/cg_dpd/`、 ~1.9s で全 PASS):
  - models + I/O round-trip 12 件 (chi → a 変換、 missing file エラー、 symmetric lookup 含む)
  - R1 UDF writer 3 件 (cholesterol e2e、 chi mode auto-convert、 custom include)
  - R2 DPM writer 7 件 (propagate_virtual_mom、 patch_dpm 5 block、 invalid template/field)
  - 多 monomer 5 件 (`from_multi_files` e2e、 missing key エラー、 aij mismatch warning、
    CLI ``--multi-monomer`` JSON、 ``--monomer``/``--multi-monomer`` 排他)
  - orchestrator + CLI 8 件 (`build-udf` / `build-dpm` subprocess test 含む)
  - a→chi 逆変換 + validate / verify 11 件 (to_chi_values round-trip、 out_mode 強制変換、
    validate ok/missing/extra、 verify CLI return 0/1)
  - notebook_ui 3 件 (import + signature + panel 構築 with ipywidgets)
  - helper (create_empty_aij + build_monomer) 10 件 (default / chi mode / write-roundtrip /
    1-粒子 solvent / dimer / trimer / custom params / CGDpdBuilder 連携)
- **サンプル**: `sample/cg_dpd/cholesterol/` に cholesterol → R1 UDF の生成例 (再現コマンド付き)。

### Added (`abmptools.hbond` 拡張 — v1.26.0 候補)

- **FF 抽象化** (`func_tags.py`): GAFF2/OPLS-AA/CHARMM36/OpenFF の 4 force field
  に対応。atom type → 機能タグ (`carbonyl_C`/`hydroxyl_O`/`amide_N` 等) の
  マッピング辞書 + 自動 FF 検出 (`detect_force_field`) + ユーザ拡張
  (`add_mapping`)。`functional_groups.py` を tag ベースに refactor。
- **任意官能基対選択**: donor `{carboxyl, amide_donor, amine_donor, hydroxyl}` ×
  acceptor `{carboxyl_O, amide_O, hydroxyl_O, ether_O}` から CLI/Python API/
  Jupyter UI で自由に組み合わせ可能。デフォルトは v1.25.0 互換 (COOH→{COOH O,
  amide O})。
- **Secondary amide donor 対応**: `AmideGroup.tert` フラグで tertiary/secondary
  を区別、`detect_amine_donors()` で N-H を donor 集合に取得可能。
  peptide 主鎖 H-bond network の解析が可能になった。
- **Lifetime + autocorrelation** (`lifetime.py`):
  - Continuous lifetime: 連続存在区間の strict 集計
  - Intermittent lifetime: `gap_tolerance` で許容 gap 指定
  - Luzar-Chandler 自己相関 `C(t) = <h(0)h(t)>/<h(0)>` (unbiased estimator)
  - τ_HB (= ∫C(t)dt) 台形則積分
  - Multi-record CLI: `--gap-tolerance N --dt FLOAT --autocorr-max-lag N`
  - 追加出力: `<prefix>_lifetime.csv` + `<prefix>_autocorr.{csv,png}`
- **Jupyter UI 拡張**: donor/acceptor 官能基チェックボックス, lifetime 設定 box,
  FF 自動検出表示, sec amide N-H donor 数表示。
- **テスト 22 件追加** (合計 42、全 PASS):
  - `test_func_tags.py` (9): 各 FF mapping + auto-detect
  - `test_lifetime.py` (8): continuous/intermittent/autocorr/τ_HB
  - `test_amine_donor.py` (5): synthetic N-methylacetamide (GAFF2/CHARMM36) + IMC 否定確認

### Added (`abmptools.hbond` サブパッケージ — v1.25.0 候補)

- **非晶質 MD トラジェクトリ用 H-bond 解析**: COGNAC UDF/BDF を入力に、
  カルボキシル基 (COOH) 同士の dual H-bond (環状二量体) と COOH→アミド C=O の
  single H-bond を幾何条件で検出・分類し、gourmet で 3 色可視化できる UDF
  を出力する新規サブパッケージ。
  - **官能基自動検出** (`functional_groups.py`): GAFF2 atomtype (`c`/`oh`/`ho`/`o`/`n`)
    + bond graph で carboxyl / amide / hydroxyl を検出。Tertiary amide (N-H なし)
    も `tert=True` でマーキング (インドメタシン対応)。
  - **Luzar-Chandler 幾何判定** (`hbond_detector.py`): `d(D-A) ≤ 3.5 Å` かつ
    `∠(D-H-A) ≥ 120°` を default、`strict` モード (`d(H-A) ≤ 2.5 Å` かつ
    `∠ ≥ 150°`) と `custom` モード (任意閾値) も選択可能。直交 cubic box の
    最短像法で PBC 対応。
  - **dual/single/free 分類** (`classifier.py`): 各分子に 3 役割を割り当て
    (優先度 dual > single > free)。Dual は両方向 COOH↔COOH H-bond が成立する
    分子ペアのみ。
  - **gourmet 色付け** (`colorizer.py`): `Set_of_Molecules.molecule[i].Mol_Name`
    を `IMC_DUAL` / `IMC_SINGLE` / `IMC_FREE` にリネームし、
    `Draw_Attributes.Molecule[]` に Red / Blue / Gray の named color を書き込む。
    **GOURMET の Draw_Attributes color は select 型 (9 色名のみ) で RGBA tuple
    は受け付けない** ことを実機検証で確認。`transparency` は 1.0 = 不透明
    (直観に反する) も確認済。
  - **3 経路インターフェース**:
    - CLI: `python -m abmptools.hbond <bdf> --criteria luzar-chandler -o prefix`
    - Python API: `from abmptools.hbond import Analyzer, AnalyzerConfig`
    - Jupyter UI: `open_panel(bdf_path)` (RDKit 2D 構造図 + ipywidgets コンパネ
      + matplotlib count plot のインライン表示)
  - **出力**: per-record summary CSV + H-bond pair CSV + colored BDF + count PNG
- **IMC amorphous サンプル**: `sample/hbond/imc_amorphous/` に CLI/notebook
  ワンライナーと期待値 (T=450 K, 125 IMC: dual=10, single=73, free=42)
- **テスト**: 20 unit + integration テスト
  - 8 角度パターン + PBC wrap geometry の検出器単体テスト
  - IMC count regression (±2-5 mol tolerance)
  - colorize round-trip 検証
- **依存**: `extras_require['hbond'] = ["matplotlib>=3.5"]` 追加。コア機能は
  `numpy` + `UDFManager` (OCTA 同梱) のみ。Jupyter UI は `[jupyter]` +
  `[fragmenter]` (rdkit) を別途併用。
- **ドキュメント**: `docs/hbond.md`、`README.md` (このセクション)

### Added (`abmptools.fragmenter.cg_segmenter` DPDgen export — v1.24.0 候補)

- **DPDgen 入力生成**: CG segments から [DPDgen](https://github.com/kojioku/dpdgen)
  (Koji Okuwaki 作の DPD UDF 生成ツール、OCTA COGNAC エコシステム) 用の
  `{name}_monomer` + `{name}_calc_sett` 2 ファイルを生成する `dpdgen_exporter.py`
  サブモジュールを追加。
  - **bond ポテンシャル (距離制約)**:
    - **bond12** (path 1): 共有 atom (fused boundary) または直接 atom bond で検出。
      - 両 segment が ring* kind → distance **0.60** / stiff 200 (cholesterol 流儀)
      - その他 → distance **0.86** / stiff 50 (DPD 均等配置 default)
      - **boundary atom filter**: 複数 segment に shared assigned された atom は
        has_direct check から除外し、 cholesterol B 環内部 bond が A-C / B-D の
        誤検出を起こすのを防ぐ
    - **bond13_150** (新): fused ring chain の BFS path 2 (= 1-skip ring) を検出。
      distance **1.661** (= 2·0.86·sin(75°)、 余弦定理 150° 仮定) / stiff 200。
      cholesterol で A-C, B-D。
    - **bond14_150** (新): fused ring chain の BFS path 3 (= 2-skip ring) を検出。
      distance **2.502** (= 4 環 chain extension, 150°/165° 想定) / stiff 200。
      cholesterol で A-D。
  - **angle ポテンシャル (角度制約)** — `bond13_180 / bond13_120` (bond ポテンシャル)
    を **置き換え**:
    - `bond12` graph path length 2 のペア `(a, b, c)` を抽出、 b を中央 segment に。
    - `angle13 = [[a, b, c], ...]` + `angle13data = [[a, b, c, eq_余角, stiffness], ...]`
      の DPDgen format (cognac 流) で出力。
    - **平衡角は cognac 流の余角 (180° - θ) で指定**:
      - 両端 ring + 中央 ring → eq **30** (= 150° ring bend, cholesterol-like)
      - cis 二重結合 周辺 (RDKit `BondType.DOUBLE`) → eq **60** (= 120°)
      - その他 (chain / mixed) → eq **0** (= 180° 直線想定)
    - stiffness 一律 **5.0** (DPD 用やや弱め)。
  - **オプション (コメントアウト template)**: `bond13_180h` (dist 1.72) /
    `bond13_120h` (dist 1.49) は default 生成せず、 距離制約も併用したい場合に
    ユーザーが monomer file で手動 uncomment する設計。
  - **aij.dat 自体は生成しない** (calc_sett に file path のみ書き込み、χ パラメータは
    `fcews-manybody` 等で別途計算済の aij.dat を配置する)
  - `CGSegmenter.export_dpdgen(...)` メソッド + CLI `dpdgen` subcommand +
    Jupyter UI `[Export DPDgen]` ボタンの 3 経路から呼び出し可能
  - tests: 11 件 (`test_dpdgen.py`、 path hierarchy / angle eq 30/60/0 / cis 検証
    / bond13 オプション コメントアウト検証 含む)、計 29/29 PASS in 0.46s
  - **検証 (cholesterol)**:
    - bond12 = 4 pair (A-B, B-C, C-D, A-tail)
    - bond13_150 = 2 pair (A-C, B-D), distance 1.661 / stiff 200
    - bond14_150 = 1 pair (A-D), distance 2.502 / stiff 200
    - angle13 = 3 entry (B-A-tail eq=0, A-B-C eq=30, B-C-D eq=30), stiff 5.0
    - bond13_180 / bond13_120 = コメントアウト template のみ

### Added (`abmptools.fragmenter.cg_segmenter` Jupyter UI 拡張 — v1.24.0 候補)

- **Jupyter UI (`open_panel`)**: fragmenter UI と同じ操作感で CG segment を
  interactive 編集できる ipywidgets パネル。
  - `notebook_ui.py` (新規 ~290 行): SVG output / Segments list / Move atom /
    Re-segment / Export ボタン
  - `exporter.render_segments_svg` (新規): 各 segment を別色 (palette 10 色) で
    highlight、shared atom を黒縁取り (`bold_shared=True`)、atom 番号 + shared `*`
    注記 (`show_atom_numbers=True`)
  - `CGSegmenter.{move_atom, toggle_cap, delete_segment, re_segment, _recompute_caps}`
    の edit API を追加。すべて cap atom を自動再計算する。
  - **move_atom**: exclusive default + `shared=True` で fused ring 化
  - **toggle_cap**: cap_index 単位で H ↔ CH3 切替、位置も新結合長で再投影
  - **delete_segment**: 任意 segment を完全削除 (atoms はどの seg にも属さない状態)
  - **re_segment**: `target_mw` / 3 flags 変更で全 segment 再計算 (上書き)
  - tests: 7 件追加 (`test_edit.py`)、計 18/18 PASS in 0.42s

### Added (`abmptools.fragmenter.cg_segmenter` 新サブモジュール — v1.24.0 候補)

- **CG (粗視化) セグメント構築ツール**: FMO 用の `fragmenter` (BDA/BAA で擬似分割) に対し、
  本サブモジュールは **物理的に分子を分割し、cap atom (H or CH3) を付与** する。
  CG MD の前処理として、コレステロール等の環構造を含む系で「環ごとの segment +
  fused ring atom 共有 + chain MW walk 切断 + cap 自動配置」を 1 コマンドで行う。
  - `ring_detector.py`: RDKit `GetRingInfo().AtomRings()` で SSSR ring 抽出、fused
    の共有 atom を両 Segment に含める (`allow_atom_sharing=True`)、ring atom の 1
    heavy 置換基 (-OH/-NH2/-F 等) を吸収 (`absorb_single_substituent=True`)
  - `chain_splitter.py`: ring 外 chain を connected component に分け、各 component を
    `auto_split._atom_total_mw` (流用) ベースの target_mw walk で切断
  - `cap_attach.py`: 境界 atom が C/halogen → H cap、N/O/S/P → CH3 cap (central C +
    tetrahedral 3 H)。位置は元 bond 方向ベクトル × 結合長で配置。
    Hetero に CH3 cap を選ぶのは「不要な水素結合スポット (擬似 H-bond) を MD で
    出さない」ため (user-specified rule)。
  - `exporter.py`: per-segment **PDB + XYZ** + summary JSON (`shared_atom_pairs`
    を含む)。1 atom が複数 segment に属することを許容する。
  - `orchestrator.py`: `CGSegmenter` クラスで pipeline 統合
    (`from_pdb()` → `segment()` → `export()`)。
  - `__main__.py`: `python -m abmptools.fragmenter.cg_segmenter {build,example}` CLI。
  - `tests/cg_segmenter/test_basic.py`: 11 tests (dataclass roundtrip /
    propane / octane / methyl acetate / benzene / naphthalene / cyclohexanol /
    cyclohexyl-octyl / 出力ファイル / CLI) PASS in 0.33s。
  - `docs/cg_segmenter.md` (新規)、`README.md` / `docs/overview.md` 更新。
- `sample/fragmenter/eude_block55.pdb`: Eudragit E mimic block 共重合体 (PMMA×5 +
  DMAEMA×5、heavy=90、MW=1288.7)。側鎖ありビニル系ポリマー (ester + tertiary amine)
  のテスト用、cg_segmenter で 7 chain segments に分割される。

### Added (`abmptools.membrane` sample)

- **`sample/membrane/amber_phaseD/`** — AMBER backend
  (ff19SB + Lipid21 + GAFF2 + TIP3P) 用 umbrella sampling 入力サンプル。
  既存 `sample/membrane/charmm_phaseD/` (CHARMM36 backend) と完全に同じ
  系構成・プロトコル (poly-Ala 5-mer + POPC 32×2 + 0.15 M NaCl、
  13 windows × 1 ns、`window_spacing_nm=0.25`) で力場のみ差し替えた、
  Phase D = L9 verification の AMBER ベースライン側 reference。
  - `input/config_phaseD.json` — `MembraneConfig` JSON
    (`backend="amber"`, `charmm_ff_dir=""`)
  - `README.md` — 概要 / CHARMM 版との差分表 / 結果サマリ
    (PMF +86.7 kJ/mol、CHARMM +97.9 と Δ-11.3 kJ/mol、典型 FF 差) /
    5-stage パイプライン (`packmol-memgen` → `tleap` → `parmed` →
    GROMACS → `gmx wham`) / 実行例 / 重量級結果保管先
  - 商用利用クリーン: `packmol-memgen` + `tleap` + `parmed` 経路で
    CGenFF / CHARMM-GUI 非依存

### Added (`abmptools.fragmenter` 拡張)

- **BDA/BAA 役割の決定ロジック**: 各切断 bond で `(BDA, BAA)` を decide し、
  `CutSite.bda_atom_idx` / `baa_atom_idx` に格納。`segment_data.dat` 出力は
  ABINIT-MP の log2config 形式に完全準拠 (BDA holder fragment のみ
  `connect_num=1`、BAA 側は 0、`connect=[(BDA, BAA)]`)。決定ルール:
  - **C-C 単結合**: walk 起点側 (graph diameter path[0] 含む側) を BDA
    (アミノ酸 N→C 方向のアナロジー、deterministic)
  - **C-X 単結合** (`include_c_heteroatom=True` 時): **C 側を BDA** (ABINIT-MP
    慣習、ユーザー指定通り)
  - **peptide N→C**: 本サブパッケージ対象外 (`abmptools.log2config` 経路で対応)
- **SVG カラースキーム強化**: `_render_svg` で BDA atom (青) / BAA atom
  (オレンジ) / cut bond (赤) を色分けハイライト。CLI ヘッドレス経路と
  Jupyter UI 両方に効く。`AllChem.Compute2DCoords` で 2D 座標を再計算する
  default 化も同時導入 (PE / PP のような長鎖が zigzag 構造式として綺麗に表示)。
- **`include_c_heteroatom` config option**: C-X (X=N/O/S/P/F/Cl/Br/I) 単結合
  切断を opt-in で許可。CLI flag は `--include-c-heteroatom`。
  `exclude_heteroneighbor` フィルタは C-C bond にのみ適用される。
- **Jupyter UI 拡張** (`open_panel`):
  - **Show bond numbers** checkbox: SVG に heavy_mol の bond_idx を `bondNote`
    で表示 (Add cut で参照)
  - **Per-cut Remove button**: enable/disable とは別に CutSite を完全削除
  - **Add cut**: bond_idx (heavy_mol) を入力欄 + "Add cut" ボタンで新規 cut
    を追加。BDA/BAA は新 helper `auto_split.decide_bda_baa_for_manual_cut`
    で自動 decide (C-X = C 側 / C-C = atom_idx 若い側)
  - **Re-suggest**: Target MW を変更 + ボタンで `suggest_cuts()` 再実行
    (全 cut を上書き)
- **SVG カラースキーム refinement** (P10 → 後追い更新): solid blue/orange
  だと BDA/BAA の区別が瞬時にしづらいというユーザー指摘を受け、化学教科書
  慣習に近い見た目に変更。**BDA は青の点線円 (塗りつぶしなし)**、
  **BAA はオレンジ塗り**、両端に `atomNote='BDA'`/`'BAA'` 文字を併記。
  実装: BDA を `highlightAtoms` から外し、`drawer.GetDrawCoords()` で取得
  した座標に SVG `<circle ... fill="none" stroke-dasharray="3,3"/>` を
  post-process で挿入。

### Fixed (`abmptools.fragmenter`)

- **`_atom_total_mw` の H カウント抜け修正**: PDB 由来の Mol は H を
  explicit atom として持つため、`atom.GetTotalNumHs(includeNeighbors=False)` は
  implicit H = 0 を返す → CH3 の MW が ~12 (本来 ~15) と過少評価され、
  累積 walk が早く target_mw に達して cut 位置が中央より 1 bond 左にずれていた。
  implicit + explicit H neighbors の両方を加算するよう修正、結果として PP N=10
  の cut が主鎖中央 (atom 13-15) に正しく配置される。
- **PP N=10 sample SMILES 修正**: 末尾に余分な `C` が混入していた
  (`"CC(C)" + "CC(C)" * 9 + "C"` = 主鎖 21 C, 31 heavy)。正しい
  `"CC(C)" * 10` に直し、主鎖 20 C / 30 heavy / 92 atoms / MW 422.8 に。

### Changed

- **License migration: MIT → Apache-2.0** (適用範囲: v1.23.0 以降の新規
  リリース)。abmptools の citation 基盤を強化するため、Apache-2.0 の
  NOTICE-file attribution 機構 (§4(d)) と特許 retaliation (§3) を採用。
  - `LICENSE` を MIT (24 行) から **Apache License, Version 2.0** (200 行
    appendix 付き、Copyright 2022-2026 Koji Okuwaki) に差し替え
  - 新規 `NOTICE` を追加 (attribution + 引用依頼 + 第三者ソフトウェアへの
    pointer)。Apache-2.0 §4(d) により再配布時の verbatim 保持が必要
  - 新規 `CITATION.cff` を追加 (GitHub の "Cite this repository" ボタン
    対応、Zenodo DOI は最初の release tag で発行・追記予定)
  - `pyproject.toml` の `license = "MIT"` を `license = "Apache-2.0"` に
    変更 (PEP 639 SPDX expression)
  - `README.md` に **License** / **How to cite** セクションを追加
  - `docs/licenses_third_party.md` 本体 license 行を `BSD 2-Clause` (誤記)
    から `Apache-2.0 (v1.23.0+; ≤ v1.22.0 は MIT)` に修正、互換性表 +
    結論セクションを Apache-2.0 ベースに更新 (NOTICE 保持義務 / 特許
    retaliation / GPL-2.0 only との incompatibility 注記を追加)
  - `docs/{cg_membrane,cg_peptide,grest,mmgbsa,fragmenter,peptide_builders}.md`
    + `README.md` の "abmptools 本体 MIT" 系の記述を Apache-2.0 ベースに
    一括更新
- **理由**: 単独開発者 (Koji Okuwaki) が abmptools を学術論文で引用して
  もらうための仕組み化。MIT は permissive すぎて attribution が薄く、
  Apache-2.0 の NOTICE 機構 + CITATION.cff + Zenodo DOI の 3 点セットで
  citation を強制する。
- **依存ライブラリとの互換性 (確認済み)**:
  - Apache-2.0 ⇔ GPL-2.0 only (insane): ❌ incompatible だが subprocess
    only なので mere aggregation = 実害なし
  - Apache-2.0 ⇔ GPL-2.0-or-later / GPL-3.0 / LGPL-2.1+ / LGPL-3.0+: ✅
    全て互換 (subprocess / import いずれでも問題なし)
  - Apache-2.0 ⇔ MIT / BSD: ✅ 互換 (依存ライブラリとして自由に取り込み可)
- **Per-file copyright header** (Apache-2.0 appendix の `[optional]` 表記)
  は **追加しない**: single developer + LICENSE/NOTICE/CITATION.cff の 3
  点セットで法的に十分。crystal 等の新規モジュールでも当面 header は無し。
- **Contributor 同意**: 単独開発者のため外部 contributor の同意取得は不要。
- 旧版 (≤ v1.22.0) は MIT のまま PyPI に固定済 (`pip install abmptools==1.22.0`
  などで取得可能)。Apache-2.0 化は v1.23.0 以降の新規取得分のみに適用。

### Added (Phase A — `abmptools.crystal` skeleton)

新サブパッケージ `abmptools.crystal` の **Phase A skeleton** を追加。
有機物結晶 CIF から FMO 計算入力 (ajf) を生成するワークフローを
`abmptools.{readcif,pdb2fmo,ajf2config,pdbmodify,getifiepieda}` の flat
配置から段階的にサブパッケージ化していく多段リリースの 1 段目。

- **`abmptools.crystal` skeleton** (cg.peptide / genesis.mmgbsa と同構造)
  - `abmptools/crystal/__init__.py` — docstring + `legacy` namespace re-export
  - `abmptools/crystal/__main__.py` — `python -m abmptools.crystal` entry point
  - `abmptools/crystal/_subprocess.py` — `cg.peptide._subprocess` の薄い再 export
    + crystal logger 用 `setup_logging`
  - `abmptools/crystal/legacy/__init__.py` — flat 配置の 5 module
    (`readcif`, `pdb2fmo`, `ajf2config`, `pdbmodify`, `getifiepieda`) を
    `abmptools.crystal.legacy.<name>` で参照可能にする薄い再 export
    (Phase C で本格的なリファクタを行うまで挙動は完全保持)
  - `abmptools/crystal/forcefield_check.py` — `ase` / `pyyaml` / `abinitmp_*`
    / `mkinp_openver1rev20.py` の presence check (Phase A はすべて optional)
  - `abmptools/crystal/cli.py` — `example` / `validate` の 2 subcommand
    (Phase C で `expand` / `fragment` / `jobs` / `pipeline` / `postproc` /
    `nearest` を追加)
- **`pyproject.toml`**
  - `[project.optional-dependencies] crystal = ["ase>=3.22", "pyyaml>=6.0"]`
    extras を追加 (Phase A は使わないが Phase C 以降必須)
  - `[project.scripts] abmp-crystal = "abmptools.crystal.cli:main"`
- **`abmptools/__init__.py`** — try/except で `from . import crystal` を追加
  (`amorphous` / `fragmenter` と同パターン、未 install でも import 失敗しない)

### Notes

- 既存 CLI (`python -m abmptools.readcif` 等) と既存 Python API は
  完全に維持。`abmptools.crystal.legacy.readcif` 経由でも同じ module
  を取得可能 (re-export)。
- 既存テスト 1302 件は全て pass (Phase A の skeleton 追加で挙動変更ゼロ)。

### Added (Phase B — regression fixture & smoke sample)

- **`tests/regression/reference/main/crystal_csp7/`** (~700 KB)
  - `R00001` / `R00002` / `R00004` — 3 構造の入力 CIF + 出力
    (`for_abmp/*layer5Zp1-around_ar6.0.{ajf,pdb}`) を fixture 化
  - 共有ドライバ (`input_param`, `segment_data.dat`, `UNK.ajf`) は
    fixture root に同梱
  - 採取は v1.22.0 の `python -m abmptools.readcif -an 32 -l 5` →
    XYZ を pdb 横に copy →
    `python -m abmptools.pdb2fmo -p input_param -xyz` 経由
    (**直接座標 AJF モード**: `&XYZ` block にフル精度浮動小数点で
    座標を埋め込み、PDB 経由の `%8.3f` 切捨を回避)。Phase C 以降の
    リファクタはこの fixture との byte-equivalence を維持すること
  - `R00001layer5Zp1` の `&XYZ` block (832 行 + ヘッダ) は
    `abmptools-sample/.../csp7_ciftest/cifout/layer5/pdb/for_abmp/` の
    既存出力と byte-equivalent を確認済み (namelist 部分のみ
    abmptools 1.22.0 で `&COUPLING/&CIS/&CISGRD/&GF2` 追加・
    `CPFVER=10` skip 化の正常進化差分あり)
- **`tests/test_crystal_regression.py`** — 4 test (3 構造 × pipeline +
  legacy namespace identity)、いずれも `@pytest.mark.slow` 付き。
  `_compare_output_dir` で for_abmp/ 配下を ajf+pdb 両方 byte-compare
- **`sample/crystal/csp7_smoke/`** — 1 構造 (R00001) 分の最小 smoke
  pipeline (`run.sh` で readcif → pdb2fmo を 5 秒で完走、出力は
  fixture と byte-equivalent)。README で Phase C 以降への移行手順も
  併記

### Added (Phase C-1〜C-4 — v1.24.0 候補の前段)

Phase C は `abmptools.crystal` 本体実装。`abmptoolsenv` に
`ase>=3.22` (LGPL-2.1+) を追加し、サブパッケージに **8 dataclass +
ASE バックエンド + HPC ジョブテンプレ + 距離検索ユーティリティ** を
配置 (CrystalOrchestrator と CLI 7-subcommand 拡張は次フェーズ)。
54 tests 全 pass、既存テストの regression ゼロ。

- **`abmptools/crystal/models.py`** — 7 leaf dataclass
  (`CIFInputSpec` / `CIFEngineConfig` / `FragmentTemplate` / `FMOMethod` /
  `HPCJobSpec` / `PostProcessSpec`) + top-level `CrystalBuildConfig`、
  `to_yaml/from_yaml/to_json/from_json` round-trip 対応。`is_xyz=True`
  を `FMOMethod` のデフォルトに (Phase B での精度確認結果を反映)。
- **`abmptools/crystal/cif_engine_legacy.py`** — `run_legacy(cif, layer,
  atoms_in_mol, odir, cwd)` で `abmptools.readcif.run_legacy_cif_pipeline`
  を Python API として呼び出すアダプタ。
- **`abmptools/crystal/cif_engine_ase.py`** — ASE バックエンド:
  `read_cif_to_atoms` (内部で対称展開済み Atoms 取得) /
  `expand_supercell` (`Atoms.repeat((N,N,N))` で `layer³` cells、
  legacy と原子数一致) / `detect_molecules` (`ase.neighborlist` +
  union-find で connected components 分割、`atoms_in_mol` post-validate) /
  `run_ase` orchestrator。csp7 R00001/R00002/R00004 (P21/n / Z=8 / P-1)
  の 3 空間群で動作確認済み。
- **`abmptools/crystal/job_templates.py`** — PJM/SLURM/PBS/local の 4
  scheduler を `string.Template` で抽象化。`render_jobscript(spec, ajf)` +
  `write_batch_runner(specs, output_dir, submit_command)` でシェル
  シェル脱出 (`$$VAR`) も整理。`spec.template_override` で任意
  テンプレ差し替え可能。
- **`abmptools/crystal/atom_distance.py`** —
  `find_nearest_atoms(pdb, center_res_seq, n_neighbors)` で
  PDB centroid からの近接原子検索 (`tips/pdbtips/readatomdistpdb.py`
  を最小 API として復元)、独立 PDB パーサで重い依存なし。
- **`abmptools/readcif.py`** — `__main__` ブロック (line 355-994) を
  `_build_parser()` + `run_legacy_cif_pipeline(args)` の 2 関数に抽出。
  挙動完全保持 (Phase B regression 4 test で担保)。
- **`pyproject.toml`** — `[crystal]` extras に `ase>=3.22, pyyaml>=6.0`
  を実体追加 (Phase A は空セクションだった)。
- **新規 unit test (54 個、`pytest tests/test_crystal_*.py`)**:
  - 25 × `test_crystal_models` — 8 dataclass の round-trip + validation
  - 2 × `test_crystal_cif_legacy` — readcif の API 経由呼び出し
  - 9 × `test_crystal_cif_ase` — read/repeat/detect/run_ase
  - 8 × `test_crystal_job_templates` — PJM/SLURM/PBS/local + override + batch
  - 6 × `test_crystal_atom_distance` — synthetic PDB nearest-atom
  - (4 × `test_crystal_regression` は Phase B で配置済み)

### Added (Phase C-5/C-6 — v1.24.0 候補完成)

Phase C 完成。**`abmp-crystal pipeline --config crystal.yaml` 1
コマンドで CIF → for_abmp/*.{ajf,pdb} + HPC ジョブスクリプト** が
生成される。csp7 R00001 で Phase B fixture と完全 byte-equivalent
(ajf/pdb cmp exit 0)。

- **`abmptools/crystal/builder.py`** — `CrystalOrchestrator` クラス、
  5 stage の orchestration:
  1. `expand_cif()` — `cif_engine_legacy.run_legacy` 経由で
     CIF → cifout/layer<L>/{pdb,xyz}/
  2. `generate_fmo()` — XYZ + drivers を pdb dir に staging、
     `abmptools.pdb2fmo.run_pdb2fmo` を **in-process 呼び出し**
     (param dict を直接渡し、interpreter 再起動を避ける)。CLI 互換は
     `pdb2fmo.main()` が同 API を呼ぶため完全保持
  3. `write_jobs()` — `job_templates.render_jobscript` で per-AJF
     jobscript 生成、`write_batch_runner` で submitter
  4. `run_abinit()` — Phase D 実装予定 (NotImplementedError stub)
  5. `postprocess()` — `postproc.run_postprocess` (Phase C-6 で実装)
  - drivers (`input_param`, `segment_data.dat`, `UNK.ajf`) は
    `_emit_input_param` / `_emit_segment_data` で config から自動生成
    (UNK.ajf テンプレのみ `FragmentTemplate.template_ajf` で path 指定
    必須)
- **`abmptools/crystal/postproc.py`** — IFIE/PIEDA + 最近接原子統合:
  `run_getifiepieda` で `subprocess` 経由 `getifiepieda` 起動 →
  `csv/` 取得、`annotate_with_nearest_atoms` で nearest-atom 列を
  追記、`run_postprocess` でこの 2 段を統合
- **`abmptools/crystal/cli.py`** — Phase A の 2 cmd から **8 cmd** に拡張:
  `example` / `validate` / `expand` (Stage 1) / `fragment` (1+2) /
  `jobs` (1+2+3) / `pipeline` (全部 + `--run-local` Phase D) /
  `postproc` (Stage 5 standalone) / `nearest` (距離検索 standalone)。
  YAML/JSON config 両対応 (拡張子で dispatch)
- **`abmptools/crystal/__init__.py`** — `CrystalOrchestrator` と
  全 7 dataclass を eager export
- **`sample/crystal/csp7_smoke/crystal.yaml`** — Phase C+ の推奨
  YAML config 例。1 コマンドで `run.sh` (Phase A/B legacy) と
  byte-equivalent な出力 + PJM jobscript + runbatch.sh
- **新規 unit test (10 個追加、累計 64)**:
  - 5 × `test_crystal_builder` — `_emit_input_param` / `_emit_segment_data`
    + R00001 byte-equivalence end-to-end + Phase D/C-6 stub raises
  - 5 × `test_crystal_cli` — example YAML round-trip + validate +
    `--help` (8 subcommand 全列挙) + nearest standalone +
    YAML-driven pipeline byte-equivalence

### Notes (Phase C 全体)

- v1.22.0 → v1.24.0 候補で `abmptools.crystal` サブパッケージ完成
  (Phase A/B = v1.23.0、Phase C = v1.24.0)
- 64 crystal tests 全 pass、既存 test の regression ゼロ
  (1541 pass、12 fail / 2 error は HEAD でも同じ pre-existing)
- **`pdb2fmo.py` には in-process API `run_pdb2fmo` を追加**
  (`main()` から `run_pdb2fmo` + `_resolve_oname` + `_load_param_file`
  に分解、CLI 動作は完全保持)。同パッケージ内のサブモジュール
  (`abmptools.crystal.builder`) は subprocess なしの直接呼び出しに
  切替済み (smoke 実測: csp7 R00001 1 構造で 1.7s、subprocess 版より
  高速)
- `readcif.py` の唯一の変更点は `__main__` ブロックの関数化
  (`_build_parser` + `run_legacy_cif_pipeline` に抽出、CLI 動作は完全保持)
- ase バックエンドは Phase C-3 で **read pipeline のみ実装**
  (read_cif_to_atoms / expand_supercell / detect_molecules /
  run_ase)、PDB emit + AJF 生成は Phase D 候題 (`engine='ase'` で
  pipeline を呼ぶと NotImplementedError)

### Added (Phase D — v1.25.0 候補完成)

Phase D で `abmp-crystal pipeline --run-local` の実機 abinitmp 実行と
ase engine 経由の計算入力生成を実装。実機 smoke (abinitmp v2r8 + csp7
R00001 layer2 + HF/STO-3G) で end-to-end 動作確認済み。

- **`CrystalOrchestrator.run_abinit()`** (D-1) — `--run-local` 実装。
  abinitmp バイナリの 3 段階解決 (`abinit_dir/binary_name` → PATH の
  `binary_name` → PATH の `abinitmp`)、`subprocess.run([abinit],
  stdin=ajf, stdout=log, stderr=err)` で per-AJF 実行、`fail_fast` で
  非ゼロ exit 時の挙動制御。Phase D 以前の `NotImplementedError` stub
  を本実装に差し替え
- **`cif_engine_ase.write_pdb_for_abmp` / `write_xyz_for_abmp` /
  `run_ase_pipeline`** (D-2) — ASE backend 経由で legacy 互換の
  PDB + XYZ を出力。`builder.expand_cif()` の `engine='ase'` 分岐を
  実装、`pdb2fmo -xyz` 経由で valid な FMO 入力 (`Natom > 0` +
  matching `&XYZ` block) まで貫通。csp7 R00001 で smoke 確認済み
  (legacy: Natom=832 / ase: Natom=704、supercell origin と分子順序
  差は設計上想定済みで byte-equivalence は legacy のみ保証)
- **`docs/crystal.md`** (D-3) — 10 節 API/設計ノート (scope / config
  schema / legacy vs ASE / cutmode / `is_xyz` / HPC templates / CLI /
  output layout / failure modes / license)
- **`docs/tutorial_crystal_fmo.md`** (D-3) — 8 節 step-by-step
  tutorial (Environment / Smoke / csp7 reproduction / ASE engine /
  HPC / `--run-local` / postprocessing / Failure modes)。Phase D-3
  追補で「Establishing a numeric reference for your crystal」
  Section 8 (7 サブ節 8-1 〜 8-7) を追加し、計 9 節に拡張
- **`tests/integration/run_crystal_smoke.sh`** (D-3) — 実機 abinitmp
  smoke。`abinitmp` 不在時は SKIP、layer=2 + HF/STO-3G で 1 構造を
  pipeline + run-local 実行
- **9 個の新規 unit test (累計 73)**:
  - 5 × `test_crystal_builder` 追加 (`_resolve_abinit_binary` 3 パス /
    fake binary 実行 / 非ゼロ exit propagation)
  - 4 × `test_crystal_cif_ase` 追加 (`write_pdb_for_abmp` 列レイアウト /
    `write_xyz_for_abmp` フル精度 / `run_ase_pipeline` PDB+XYZ 出力 /
    orchestrator `engine='ase'` end-to-end)

### Notes (Phase D 完了)

- v1.25.0 候補完成。`abmp-crystal pipeline --config crystal.yaml
  --run-local` 1 コマンドで CIF → for_abmp/*.{ajf,pdb} →
  abinitmp 実行 → log 取得まで貫通
- 73 crystal tests 全 pass、実機 smoke PASS (abinitmp v2r8 + HF/STO-3G)
- ASE backend は valid な FMO 入力を生成 (legacy と byte-equivalent
  ではないが Natom と `&XYZ` block の整合性は確保)
- README に `### Crystal-FMO Pipeline (crystal)` セクション追加は
  別 commit で扱う (working tree に既存の README 変更があるため)

### Added (Phase D-3 — numeric reference & verification doc)

Phase D 完了後、abinitmp v2r8 で **csp7 R00001 layer3 HF/6-31G の
本格 FMO 計算**を実機実行 (1 core で 9h 14min 完走) し、Total energy
+ 24 monomer energies + 95 ESP-AOC IFIE 値を numeric reference として
凍結。

- **`tests/regression/reference/main/crystal_csp7/R00001/expected_layer3_hf_631g.json`**
  (21 KB) — abinitmp v2r8 reference snapshot:
  - FMO2-HF / 6-31G、24 fragments / 768 atoms / 276 dimer pairs
    (95 ESP-AOC + 181 ESP-PTC)
  - Total energy = -34450.7633976498 hartree
  - 24 monomer HF energies (~-1435.43 hartree each, ±0.005 hartree
    結晶環境による分散)
  - 95 ESP-AOC IFIE pairs (kcal/mol range -9.47 ~ +1.65)
  - tolerances: total/monomer ±1e-5 hartree, IFIE ±1e-3 kcal/mol
- **`tests/regression/reference/main/crystal_csp7/R00001/excerpt_layer3_hf_631g.log`**
  (10 KB) — full log (322 KB) からの必要セクション抜粋。test fixture
  として extract script の roundtrip 検証に使う
- **`tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py`**
  — abinitmp log から JSON reference を生成する CLI ツール。
  再現用に repo に同梱
- **`tests/test_crystal_numeric_regression.py`** (3 test、累計 63 collected /
  62 passed + 1 gated):
  - `test_extract_script_roundtrip` — excerpt log → extract → JSON
    bit-perfect equality (~0.05 s、常時実行)
  - `test_expected_json_shape` — JSON 必須フィールド検証
  - `test_live_layer3_hf_631g_against_reference` —
    `@pytest.mark.slow` + `ENABLE_FMO_LIVE_REGRESSION=1` gate、
    実 abinitmp で 9h 計算 → 凍結 reference と tolerance 比較
    (default skip)
- **`docs/crystal_verification.md`** — 7 節 verification record:
  test matrix (8 ファイル × カバー機能) / 実行コマンド / 検証済み
  環境 (Python 3.11.11 + ase 3.28.0 + abinitmp v2r8) / numeric
  reference summary / failure mode 逆引き表 / 既知の限界
  (HF/6-31G 1 構造のみ、MP2 未取得、abinitmp v2r8 OMP 並列無効) /
  reference 再現 step-by-step

### Notes (Phase D-3)

- 計算実行時の所要時間 (2026-05-07 → 2026-05-08): 9h 14min 33s
  (33253 sec, exit 0) on WSL2 + 1 core。OMP_NUM_THREADS=4 環境
  変数は abinitmp v2r8 build に効かず、MPI 並列前提と判明
- 完全 log (322 KB) は OneDrive `abmptools-dump/crystal/csp7_r00001_l3_hf_631g/`
  に退避 (memory `feedback_artifact_placement` 準拠)、repo には
  excerpt のみ commit
- `getifiepieda` 経由 CSV 化は今回 skip (log の `## HF-IFIE` block を
  Python regex で直接抽出)。CSV 経路は将来追加予定 (※ 2026-05-09 で
  下記 D-3 revised により実装、Python regex 経路は廃止)

### Changed (Phase D-3 revised — 2026-05-09: getifiepieda 経路に統合)

ユーザー方針 (「同じ機能の別モジュールは作らない」) に沿って、
`abmptools.getifiepieda` を numeric reference 抽出に再利用する形に
切替。HF log 対応のため `anlfmo.py` 5 箇所を修正。

- **`abmptools/anlfmo.py` 5 箇所修正**:
  1. `getlogorpdbfrag` (line 648-660) — `ReadGeom = ` 空欄 (`is_xyz=True`
     経路) で `items[2]` IndexError → sibling PDB basename へ fallback
  2. `readmultiifie` — `self.logMethod = self.getlogmethod(...)` +
     method 別 `self.icolumn` 設定を追加 (multi mode に欠落していた、
     `readsingleifie` は line 2608-2620 で同等処理あり)
  3. `getfiltifpifd` — HF log で MP2 関連列 (MP2-IFIE/PR-TYPE1/GRIMME/
     JUNG/HILL) が無い場合 0.0 で埋め、`KeyError` 回避
  4. `getmomenedf` — HF log の monomer 行は 2 列 (`Frag.`, `HF`)、MP2 は
     3 列。method 別に column 設定 + MP2 を 0 で pad
  5. `getdimenedf` — HF log dimer 行は 3 列、MP2 は 4 列。同様に
     method 別

- **`extract_layer3_hf_631g.py`** — Python regex 抽出を廃止し、
  `subprocess.run([... abmptools.getifiepieda ... --multi 1 -dimeres
  -imd -zp 5 -t <id> <id> 1 -nof90 -i '["pfx","sfx"]'])` の wrapper
  に書き換え。出力 CSV を `--out-dir` へ rename copy
- **新 reference CSVs (canonical getifiepieda outputs)**:
  - `tests/regression/reference/main/crystal_csp7/R00001/expected_layer3_hf_631g_ifiesum.csv`
    (250 bytes) — target frag 1 の sum 1 行 (HF-IFIE/ES/EX/CT-mix/
    `MonomerEnergy(1)` = -1435.433574 hartree 等)
  - `tests/regression/reference/main/crystal_csp7/R00001/expected_layer3_hf_631g_ifiedt.csv`
    (1797 bytes) — dimer-es=False の 16 pair 詳細、`UNK<i>(<i>)`
    resname 注釈付き
- **`tests/test_crystal_numeric_regression.py`** — excerpt log
  roundtrip を廃止。新 test:
  - `test_expected_csvs_have_committed_shape` (always run、~0.01 s):
    committed CSVs の shape sanity (列存在 / MP2 zero-padding /
    monomer energy range / dimer-es=False filter)
  - `test_live_layer3_hf_631g_against_reference`
    (`@pytest.mark.slow` + `ENABLE_FMO_LIVE_REGRESSION=1` gate):
    pipeline + run-local + extract → byte-compare with committed CSVs
- **`docs/tutorial_crystal_fmo.md` Section 8-4/8-5/8-6/9** —
  getifiepieda 経路 + in-crystal monomer energy (= MonomerEnergy 列、
  ifiesum.csv 由来) の default 推奨、isolated monomer 用別 ajf
  (cutmode='around' criteria=0.0) の生成手順、Failure modes に
  anlfmo.py 5 箇所修正点に対応する 4 行追加 (各 KeyError /
  ValueError → 修正済の旨)

### Notes (Phase D-3 revised)

- legacy `expected_layer3_hf_631g.json` (21 KB) と
  `excerpt_layer3_hf_631g.log` (10 KB) は archival のため残置。
  test 経路は読まない (削除しない方針 = memory `feedback_no_deletions`)
- production (Fugaku csp7 1500 構造、MP2/6-31Gdag) は影響なし —
  anlfmo.py 修正は HF log 用の防御を追加しただけ、MP2 path は不変
- 1 構造 + HF 6-31G という小規模で getifiepieda 経路が成立
  したことで、今後別 system (R00002/R00004 / 他 space group) も
  同じ extract script を流用可能 (zero-padding 桁違いだけ調整)
- crystal tests: `pytest tests/test_crystal_*.py` で 61 passed +
  2 skipped (numeric live + 別 1 件)、regression なし

### Changed (Phase D-3 revised, separation — 2026-05-09)

公開リポジトリ (abmptools) には csp7 構造データを含めない方針に
基づき、入力構造を **abmptools-sample (非公開)** へ分離。整形手順 +
結果数値 CSV のみ abmptools 残置。

- **abmptools-sample に移動した入力構造**:
  - `<sample>/sample/csp7_ciftest/crystal_reference/R00001_layer3_hf_631g/`
    — Phase D-3 numeric reference 用 (cif + UNK.ajf + crystal.yaml +
    run_local.sh + excerpt log + JSON snapshot)
  - `<sample>/sample/csp7_ciftest/crystal_reference/phase_b_layer5/`
    — Phase B byte-equivalence 用 (UNK.ajf + input_param +
    segment_data.dat + R00001/R00002/R00004 ごとの cif + 期待
    `*-around_ar6.0.{ajf,pdb}`)
  - 移動先 README: `<sample>/.../crystal_reference/README.md` に
    test 起動方法 (env var ABMPTOOLS_SAMPLE_DIR) を明記

- **abmptools-sample 公開分離の test 改修**:
  - `tests/test_crystal_regression.py` (Phase B): fixture path を
    `${ABMPTOOLS_SAMPLE_DIR:-~/repos/abmptools-sample}/sample/csp7_ciftest/
    crystal_reference/phase_b_layer5/` に変更。fixture 不在時 skip
  - `tests/test_crystal_numeric_regression.py::test_live_*`: 入力 cif と
    UNK.ajf の取得元を abmptools-sample 経由に変更。env var 不在時 skip
  - `tests/integration/run_crystal_smoke.sh`: 同上、cif/UNK.ajf
    取得元を abmptools-sample に変更
  - `sample/crystal/csp7_smoke/run.sh`: cif/UNK.ajf を abmptools-sample
    から auto-stage (ローカル不在時のみ)
  - `sample/crystal/csp7_smoke/README.md`: 書き換え、 cif/UNK.ajf は
    abmptools-sample 経由で stage する旨を明記

- **削除候補 (user 自身が `git rm`、memory `feedback_no_deletions` 準拠)**:
  abmptools 公開リポジトリから取り除くべき csp7 構造ファイル一覧は
  本セッションのログ末尾を参照。anlfmo.py / extract script /
  数値 CSV (`expected_*_ifiesum.csv` / `_ifiedt.csv`) は残す

### Changed (Phase D-4 — 2026-05-09: --run-local の MPI/OMP 並列対応)

`abmp-crystal pipeline --run-local` を mpirun 経由起動に書き換え。
Phase D-3 の reference 計算で `OMP_NUM_THREADS=4` が効かず 9h かかった
原因は **`abinitmp` が MPI flat 版**だったため (user 共有、memory
`reference_abinitmp_parallelism` に保存)。

- **`HPCJobSpec` に `mpi_launcher` フィールド追加** (default `"mpirun"`)。
  空文字列で direct invocation (mpirun 不在環境向け、`proc_per_node=1`
  + flat binary 限定で valid)
- **`CrystalOrchestrator._build_run_command()` 新設**: バイナリ命名規則で
  flat / 混成を判定 → cmd と env を構築:
  - `abinitmp` (suffix なし) = MPI flat: `mpirun -np <proc_per_node>
    abinitmp` 経由 (single rank + `mpi_launcher=""` の場合だけ direct)
  - `abinitmp_omp` (suffix `_omp` or `_omp_` 含む) = MPI/OMP 混成:
    `mpirun -np <proc_per_node> abinitmp_omp` + `OMP_NUM_THREADS=
    <omp_threads>` env override
  - 多 rank or 混成で `mpi_launcher=""` の場合 `ValueError` で明示拒否
- **`run_abinit()` を新 cmd builder 経由に書き換え**、`subprocess.run` の
  `env` 引数で OMP override を反映
- **5 個追加 unit test** (`tests/test_crystal_builder.py`):
  flat single-rank direct / flat multi-rank mpirun / `_omp` hybrid /
  custom launcher (`srun`) / 多 rank での空 launcher 拒否
- **既存 `test_run_abinit_with_fake_binary` / `test_run_abinit_propagates_failure`**
  は `proc_per_node=1, mpi_launcher=""` を明示して direct invocation
  path を維持
- **`tests/integration/run_crystal_smoke.sh`** の crystal.yaml に
  `mpi_launcher: ""` を追加 (smoke 用途で MPI 不在環境にも対応)
- **`docs/tutorial_crystal_fmo.md`**:
  - Section 8-2 wall time 表に **バイナリ別 launcher** 表を追加、
    `OMP_NUM_THREADS` が flat では無視される旨を明記
  - Section 8-3 の crystal.yaml テンプレに `mpi_launcher: ""` 追加
  - Section 9 (Failure modes) に MPI runtime 不在 / 空 launcher 多 rank /
    omp_threads 効かない の 3 行追加 (各々 fix 提示)

検証: `pytest tests/test_crystal_*.py` -> 62 passed, 6 skipped (5 個増、
fail なし、production csp7 1500 構造 path には影響なし)

### Added (Phase D-5 — 2026-05-10: ASE PBC unwrap + 公開分子 MP2 reference)

ASE backend の supercell 拡張で **境界またぎ分子** が放置され
ABINIT-MP の Monomer SCC が振動 (urea で ±1420 hartree、glycine で
±0.04 hartree、unwrap 前) する root issue が判明。さらに公開可能な
公開文献値分子の MP2/6-31G(d) reference を追加。

- **`abmptools/crystal/cif_engine_ase.py: unwrap_molecules()`** 新規:
  BFS along bond graph + minimum-image translation で `Atoms.repeat`
  後に分子が unit-cell 境界をまたいでいる atom を contiguous に再配置。
  `run_ase()` で `detect_molecules()` の後に自動呼び出し。
  `__all__` にも export 追加
- **新規 unit test (`tests/test_crystal_cif_ase.py`)**:
  - `test_unwrap_molecules_keeps_molecules_intact`: 5 Å unit cell 境界を
    またぐ N2 分子を `Atoms.repeat((2,2,2))` 後 unwrap → 各分子の atom
    間距離 < 1.5 Å (= 結合長) を確認
  - `test_run_ase_unwraps_molecules_for_csp7`: csp7 R00001 layer 2 で
    各 32-atom 分子の最大 intramol 距離 < 15 Å を確認
- **`tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py`**:
  `args.log.resolve()` → `args.log.absolute()` に変更し、
  symlink を dereference しない (新 sample で zero-padded structure id
  経由 getifiepieda 起動するための前提)
- **公開可能な分子 sample 3 種** を `sample/crystal/<molecule>/`
  に追加 (公開文献値ベース、ASE `crystal()` で構築):
  - `urea/` (Worsham, Levy & Peterson 1957, P-421m, Z=2):
    cif + UNK.ajf + crystal_layer3.yaml + monomer.yaml + run.stdout +
    `reference/expected_layer3_mp2_631gd_{ifiesum,ifiedt}.csv` +
    `reference/expected_isolated_monomer_mp2_631gd.txt` + README.md
    - 14 fragment / 112 atoms (layer3 cutmode around 6 Å)
    - 2.5 分で完走 (mpirun -np 4)
    - in-crystal vs isolated monomer の MP2 差 +5.2 kcal/mol
      (polarization 由来、structural deformation ではない)
  - `glycine/` (Marsh 1958, P21/n, Z=4): 同 layout、layer3 21 frag
  - `benzene/` (Cox/Cruickshank/Smith 1958, Pbca, Z=4): 同 layout
  - (`naphthalene/` は当初 hand-craft cif で座標重複だったが、Phase
    D-5 続編 (2026-05-10、後続 commit) で **COD 2311088** (Hoser &
    Madsen 2017、P21/c、293 K X-ray、CC0) を取得し追加。HF-IFIE
    -1.5 kcal/mol、MP2-IFIE -26.3 kcal/mol で benzene の 2× の dispersion
    contribution を示す aromatic stacking 例)
- **`docs/tutorial_crystal_fmo.md` Section 9 (Failure modes)** に
  ASE Monomer SCC 振動 → unwrap fix 1 行追加

検証 (`pytest tests/test_crystal_cif_ase.py -v`):
abmptoolsenv 経由で 1 unwrap test pass + 既存 13 test 維持。
公開分子 reference は実 abinitmp v2r8 で完走確認済み。

## [1.22.0] - 2026-05-06

新サブパッケージ `abmptools.genesis.mmgbsa` を追加。GENESIS atdyn の
`[ENERGY] implicit_solvent=GBSA` を使った protein-ligand 単フレーム
MM/GBSA ΔG_bind 計算を end-to-end に組める **GENESIS 系統 2 番目** の
モジュール (1.20.0 grest との並走)。

### Added

- **`abmptools.genesis.mmgbsa` — GENESIS MM/GBSA single-point ΔG_bind**
  - `MMGBSAOrchestrator.run()` の 4 stage:
    1. **PDB 分割** (Biopython): receptor (指定残基除外) + ligand (指定残基のみ)
    2. **パラメタライズ** (acpype + tleap): GAFF/GAFF2 + AM1-BCC ligand →
       complex / ligand / receptor の 3 系 prmtop+inpcrd
    3. **GBSA 単フレーム実行** (mpirun atdyn): NOBC + implicit GBSA で
       1-step minimize による single-point energy 抽出
    4. **解析** (matplotlib): `[STEP4] Compute Single Point Energy` パース
       → `ΔG_bind = (E+S)_complex − (E+S)_ligand − (E+S)_receptor` →
       CSV + 棒グラフ
  - データクラス (6 個): `TargetSpec` / `ForceFieldSet` /
    `LigandParameterization` / `EnergyProtocol` /
    `MinimizationProtocol` / `MMGBSABuildConfig`
    (`@dataclass` + `to_json/from_json`、5 段の nested JSON 往復)
  - 入力モード 2 系統 (mode "C"):
    - **JSON config** (推奨, 再現性): `targets: [...]` で N target を明示
    - **folder mode** (POC 互換): `-i input_dir -r ligand_resno [-c chain]`
      で `*.pdb` を一括処理 (`pipeline` サブコマンド)
  - CLI: `python -m abmptools.genesis.mmgbsa
    {example,validate,divide,parameterize,run,analyze,pipeline}`
    (7 sub-command、Step ごと再開可能 + 一気通貫 + folder shortcut)
  - 力場 default: AMBER **ff14SB** + DNA.OL15 + RNA.OL3 + TIP3P +
    GAFF/GAFF2 (POC 通り、grest の ff19SB とは意図的に別、MM/GBSA で
    広く使われる安定 default)
  - 既存 `abmptools.cg.peptide.forcefield_check.ToolStatus` を import
    経由で再利用、`_subprocess` も同流儀
  - 4 target POC 再現用 sample (`example_config.json`) + smoke
    (`smoke_config.json`) + README を `sample/gbsa/` 直下に配置
  - 116 unit tests (`test_mmgbsa_{models,pdb_splitter,system_builder,
    inp_writer,analysis,forcefield_check,builder_mocked,cli}.py`) +
    1 slow integration smoke (`test_mmgbsa_integration.py`)
  - **3772L gold value** (POC 実 log fixture から実測):
    `ΔG_bind = -38.2700 kcal/mol` を test で regression 化
    (`tests/fixtures/gbsa_logs/3772L-rename/{complex,ligand,receptor}.log`)。
    GENESIS doc `05_Energy.rst:564` の `U = U_FF + ΔG_solv` 規約に従い、
    `ΔG_bind = E_complex - E_ligand - E_receptor` (ENERGY 列差)。
  - 解析関数: `compute_dg_bind` (合計差) + `compute_dg_components`
    (POC 風の `{dg_mm, dg_solv, dg_bind}` 分解、3 値の和は dg_bind と一致)

### Changed

- `sample/gbsa/gbsa-input/{1_devidepdb,2_setupmd,3_rungbsa,4_analyse}.py`
  を `sample/gbsa/legacy_scripts/` に `git mv` で移動 (POC 保存)、新規
  `sample/gbsa/example_config.json` / `smoke_config.json` /
  `README.md` / `legacy_scripts/README.md` を `sample/gbsa/` 直下に配置
  (mode "C")。POC scripts は研究の再現性のため温存、新規利用者は
  abmptools.genesis.mmgbsa CLI 経由を推奨。
- POC `4_analyse.py` の ΔG_bind 計算式
  `(egas+S)_c - (egas+S)_l - (egas+S)_r` は v1.22.0 の
  `E_c - E_l - E_r` と **代数的に同等** (egas + S = ENERGY、GENESIS
  doc 05_Energy.rst:564)。POC 値と v1.22.0 値は bit-exact 一致する。

### External dependencies (new, optional `[mmgbsa]`)

- `biopython >= 1.80` (Biopython License + BSD-3-Clause): PDB splitter (Stage 1)
- `matplotlib >= 3.5` (PSF/BSD): ΔG_bind 棒グラフ (Stage 4、grest と shared)
- GENESIS atdyn >= 2.1 (**LGPL-3.0-or-later**): subprocess 経由のみ、
  未同梱、`https://github.com/genesis-release-r-ccs/genesis` から build
- AmberTools `tleap` (free academic+commercial): conda 推奨
- **acpype** (>=2022.7.21, **GPL-3.0**): subprocess 経由のみ、未同梱
  (mere aggregation per GPL FAQ; abmptools 独立)
- mpirun (OpenMPI / MPICH): `mpirun -np 1 atdyn`

### Deferred to v1.22.x / v1.23.x

- per-residue energy decomposition (POC `4_analyse.py` にも非実装)
- MD ensemble averaging (現状は単 frame minimize-only single-point のみ)
- alanine scanning
- 並列 target 処理 (multiprocessing): 現状は逐次のみ
- `forcefield = CHARMM` 経由 (現状 AMBER 想定で kcal/mol)
- explicit solvent (TIP3P box + PBC + PME) 版 MM/PBSA

## [1.21.0] - 2026-05-06

新サブパッケージ `abmptools.fragmenter` を追加。FMO 計算用のフラグメント分割を
PDB から自動提案する。タンパク質 / DNA は対象外 (既存 `log2config` 経路を維持)、
小分子 / 脂質 / ポリマーに特化した **canonical SMILES グループ化 + C-C 切断
MW walk + Jupyter / CLI 双方の UI** を提供する。リリース予定は v1.21.0。

### Added

- **`abmptools.fragmenter` — FMO automatic fragment splitter** (新サブパッケージ)
  - `models.py`: `FragmenterConfig` / `CutSite` / `MoleculeGroup` /
    `FragmentResult` データクラス (JSON 往復可能)
  - `pdb_loader.py`: PDB → RDKit Mol、4 段階 bond perception フォールバック
    (proximity → CONECT-only → sanitize=False → obabel)、連結成分分解
  - `grouping.py`: heavy-atom-only canonical SMILES でのグループ化、
    `RemoveHs` 失敗時のフォールバックあり
  - `auto_split.py`: C-C 切断候補の自動提案
    - graph diameter (heavy-atom-only 2-pass BFS) で主鎖検出
    - 主鎖を walk し、累積 MW + 側鎖 MW を加算
    - 累積 ≥ `target_mw` (default 200 g/mol) のたびに candidate bond で切断
    - フィルタ: 環内除外 / 多重結合除外 / ヘテロ隣接除外
  - `cut_apply.py`: `RWMol` で破壊的に bond 削除 → `GetMolFrags` で fragment 抽出、
    formal charge sum と BAA pairs を計算
  - `expand_to_system.py`: 1 group の cut_sites を全コピーへ展開、
    `log2config` 互換 `segment_data.dat` を出力 (`pdb2fmo` がそのまま読める)
  - `headless_io.py`: review bundle 経路 (C 経路)
    - `export_review_bundle`: `config.json` + `review.json` +
      group_NNN.svg / group_NNN.json (RDKit SVG で cut_sites 赤色ハイライト)
    - `import_edited_review`: 編集後 JSON を読み戻し、
      `resync_member_indices` で PDB との対応を再計算
  - `__main__.py`: CLI (`python -m abmptools.fragmenter {suggest,apply,example}`)
    - `--target-mw`、`--no-exclude-ring/multibond/hetero`、
      `--split-by-resname` のフラグ
  - `notebook_ui.py`: Jupyter UI (A 経路)
    - `AutoFragmenter` クラス (state 管理 orchestrator)
    - `open_panel`: ipywidgets dropdown / SVG / checkbox / export
  - `polymer.py`: γ 経路 (`declare_same_pattern`)
    - 異なる SMILES (PE N=10 / N=11 等) を明示的に同一視
    - master (= 最も cut が多い group) のパターンを atom-path-index 対応で
      短い chain にも転送
  - **依存追加**: `pyproject.toml` に
    `[fragmenter]=rdkit-pypi>=2022.09`、
    `[jupyter]=ipywidgets>=8.0,ipykernel`
  - **テスト**: `tests/fragmenter/test_basic.py` (10 tests) +
    `test_polymer.py` (4 tests)。fixtures は `conftest.py` で RDKit から
    生成 (propane / octane / propane×3+acetone×2 / PE N=10+N=11)。
  - **ドキュメント**: `docs/fragmenter.md` (本サブパッケージリファレンス)、
    `abmptools/fragmenter/README.md` (モジュール構成早見表)、
    `docs/overview.md` の "Where to Start Reading" に項目 11 として追加

## [1.20.0] - 2026-05-06

新サブパッケージ `abmptools.genesis.grest` を追加。GENESIS gREST_SSCR
(generalized Replica-Exchange with Solute Tempering — Solute Side Chain
Repartitioning) を end-to-end に組める最初の **GENESIS 系統** サブパッケージ。
AMBER ff19SB + TIP3P を tleap で勾配し、prmtop+coor を経由して
minimize → equilibrate → grest → remd_convert → 1D distance PMF まで自動化する。

### Added

- **`abmptools.genesis.grest` — GENESIS gREST_SSCR builder + analysis** (新サブパッケージ)
  - JSON 入力から GENESIS 入力一式
    (`system.tleap` / `system.prmtop` / `system.coor` /
    `step1_minimize.inp` / `step2_equilibrate.inp` /
    `step3_grest.inp` / `step5_remd_convert.inp` + `run.sh`)
    を生成する end-to-end builder。
  - `GrestBuilder.build()` で 5 stage を 1 呼び出し:
    1. **`tleap`** で AMBER prmtop + coor + reference PDB を生成
       (ff19SB + TIP3P solvateBox + 自動中和)
    2. **REST 残基確定** (`mode="explicit"` 直接指定 / `mode="around"`
       cpptraj `<:radius` mask 経由)
    3. **温度ラダー生成** (`mode="auto"` 幾何級数 / `mode="manual"` 直接指定)
    4. **`inp_writer`** で 4 種の GENESIS control file を生成
       (POC + tutorial 12.3 syntax)
    5. **`grest_runner`** で `mpirun -np {n_total} spdyn` 用 `run.sh`
       テンプレート出力 (Fugaku PJSUB / SLURM scaffold は comment-only)
  - データクラス (6 個): `RESTSelectionSpec` / `ReplicaTemperatureSpec` /
    `MinimizationStage` / `EquilibrationStage` / `GrestStage` /
    `GrestBuildConfig` (`@dataclass` + `to_json/from_json`、5 段の
    nested JSON 往復)。
  - 解析: `analyze` サブコマンドで以下を実装:
    - `remd_convert` パラメータ sort (最低 T レプリカ抽出)
    - replica transition plot (`.rem` ファイル → matplotlib)
    - acceptance ratio plot (REMD log 集計、burn-in 100 exchange skip)
    - 1D 距離 PMF (`-kT log P(r)`、cpptraj 距離時系列経由)
    - 2D PMF / k-means / MBAR は v1.20.x 以降 (NotImplementedError)
  - CLI: `python -m abmptools.genesis.grest {build,validate,example,analyze}`。
  - **GENESIS LGPL-3.0-or-later、AmberTools (free academic + commercial)**:
    binary は本パッケージ未同梱、`subprocess` で呼ぶのみ
    (mere aggregation per LGPL §5/§6)。abmptools 自体は MIT (1.20.0 時点;
    Apache-2.0 化は別途 release 後に予定)。
  - 既存 `abmptools.cg.peptide._subprocess` を import 経由で再利用、
    コード重複なし。
  - 実行は **ローカル mpirun のみ** サポート。Fugaku (pjsub) / SLURM は
    `run.sh` を真似た comment-only scaffold が同梱されるが、abmptools 側では
    実機検証しない (= ユーザー責務)。
  - 159 unit tests (test_grest_models.py / _replica_temperatures /
    _rest_selection / _inp_writer / _system_builder / _builder_mocked /
    _analysis_mocked / _cli) で 6 dataclass + 5-stage flow + 4 解析関数
    を CI 化 (tleap / spdyn / mpirun 不要)。実機 smoke は
    `tests/test_grest_integration.py` (`@pytest.mark.slow`) で gated。

### External dependencies (new, optional)

- `matplotlib >= 3.5` (PSF License/BSD-compatible): replica transition /
  acceptance / PMF プロット。`pip install abmptools[grest]` で自動 install。
- GENESIS spdyn / atdyn / remd_convert >= 2.1 (**LGPL-3.0-or-later**):
  https://github.com/genesis-release-r-ccs/genesis から build。subprocess
  経由のみ、ソース改変・同梱なし。Apache-2.0 化後も互換 (mere aggregation)。
- AmberTools `tleap` (必須) + `cpptraj` (around-mode 時に必須): conda 経由
  install (`mamba install -c conda-forge ambertools`)。
- mpirun (OpenMPI / MPICH): 並列レプリカ実行。

### Notes

- `param_type` default は `["C", "L"]` (CHARGE + LJ): POC ログ
  (`Setup_Remd_Solute_Tempering>` 出力の `CHARGE=T LJ=T`) と tutorial 12.3
  step3 の SSCR-typical 設定に合わせた。
- 温度ラダー auto mode は v1.20.0 では `geometric` のみ。Patriksson-van der
  Spoel acceptance-ratio formula は v1.20.x で追加予定 (`NotImplementedError`)。
- POC build caveat: icx で `fileio_data_.c` の `ftello64`/`fseeko64` が
  undeclared になる場合は `sed 's/ftello64/ftello/g; s/fseeko64/fseeko/g'`
  で patch。

## [1.19.0] - 2026-05-05

新サブパッケージ `abmptools.cg.membrane` を追加。Martini 3 + insane で
ペプチド-脂質膜系の **PMF (umbrella sampling 経由)** を end-to-end に組める
CG 系統 2 番目のモジュール。`cg.peptide` (1.18.0) を内部 sub-call し、
`abmptools.membrane` (AA umbrella + WHAM) の generic helper を import 経由で
直接再利用する -- コード重複ゼロ。

### Added

- **`abmptools.cg.membrane` — Martini 3 peptide-membrane US system builder** (新サブパッケージ)
  - JSON 入力から GROMACS MD 入力一式 (CG bilayer + peptide CG +
    `topol.top` + `index.ndx` + `em/nvt/npt/pull.mdp` + 13 window MDPs +
    `run.sh`) を生成する end-to-end builder。
  - `MembraneCGBuilder.build()` で 7 stage を 1 呼び出し:
    1. 4 Martini 3 ITPs を output_dir に自動コピー
    2. **`PeptideCGBuilder` sub-call** (`solvent_enabled=False`,
       `mdp_*=False`) で M3 CG ペプチド (`<name>_cg.pdb` + `.itp`) 生成
    3. **`insane`** で peptide を POPC bilayer に埋め込み + W solvent +
       NaCl (`-charge auto` で peptide 電荷を考慮した salt 配置)
    4. **`topology_composer.compose_topology`** で insane の生 topology を
       post-process (4 ITP includes 群への書き換え、`Protein` →
       `molecule_0` 置換、`NA+`/`CL-` 名の正規化)
    5. **`normalize_ion_atom_names_gro`** で .gro の ion atom name 正規化
    6. **`write_ndx_from_gro_cg`** で Bilayer / Peptide / W / NA / CL /
       Non_Bilayer の named groups を Python 直書き (gmx make_ndx の
       stdin 戦法不要)
    7. CG MDP 群 (em / nvt / npt-semiisotropic /
       direction-periodic pull / 13 window static umbrella) + `run.sh`
  - データクラス: `MembraneCGBuildConfig` / `LipidMix` /
    `PeptideMembraneSpec` / `EquilibrationCGProtocol` /
    `PullingCGProtocol` / `UmbrellaCGProtocol`
    (`@dataclass` + `to_json/from_json`、5 段の nested JSON 往復)。
  - CLI: `python -m abmptools.cg.membrane {build,validate,example,make-windows,wham}`
    (argparse、`cg.peptide` / `membrane` 流儀準拠)。
  - **GPL-2.0 / Apache-2.0 / MIT 互換のみ**: insane (GPL-2.0) と
    vermouth (Apache-2.0) を `subprocess` で呼ぶのみ、改変 / 同梱なし。
    abmptools 自体は MIT のまま (subprocess は GPL FAQ の "mere
    aggregation" 扱い、`martinize2` で確立した同流儀)。
  - 既存サブパッケージとの統合 (コード重複ゼロ):
    - Stage 1 で `abmptools.cg.peptide.PeptideCGBuilder` を
      sub-call (`solvent_enabled=False`, `mdp_*=False`) で CG ペプチドのみ
      生成。
    - Pulling helpers (`parse_pullx_xvg`, `find_pbc_center_atom`,
      `extract_window_frames`, `estimate_initial_pull_coord`,
      `_gmx_trjconv_dump`, `_override_mdp_field`,
      `_read_gro_positions`, `_read_ndx_groups`) を
      `abmptools.membrane.pulling` から **import 経由で再利用**。
    - `render_pull_block` を `abmptools.membrane.mdp_us_protocol` から
      duck-typed import (`UmbrellaCGProtocol` と `USProtocol` の関連
      フィールド名が一致)。
    - `gmx wham` 呼び出し (`run_wham`) は `abmptools.membrane.pmf` を
      delegate (`config.gmx_path` / `config.equilibration.temperature_K`
      の duck-typing で動作)。
  - **Default umbrella protocol**:
    - 13 windows (z = -1.5 to +1.5 nm, spacing 0.25 nm)
    - k = 1000 kJ/mol/nm²
    - 1 ns/window (50,000 steps × dt=20 fs)
    - 5 ns pulling (250,000 steps × dt=20 fs, 1 nm/ns)
    - AA membrane と同じ z 軸プロトコルを CG dt=20 fs に scale。
  - **Default MDP chassis** (Martini 3 標準):
    - cutoff = 1.1 nm, reaction-field, ε_r = 15
    - 2-group thermostat (Bilayer / Non_Bilayer, V-rescale, τ_t = 1 ps)
    - semiisotropic c-rescale (τ_p = 12 ps, compressibility = 3e-4)
    - constraints = none (CG だから h-bond constraint 不要)
    - pull stage は **NVT-chassis + direction-periodic** (動的 box 不可)、
      window stage は **NPT-semiisotropic + direction** (動的 box 互換)
  - `forcefield_check`: REQUIRED Martini 3 files = 4 ITP
    (`martini_v3.0.0.itp` / `_solvents_v1` / `_ions_v1` /
    **`_phospholipids_v1.itp`** -- POPC topology を含む M3 phospholipid
    一括ファイル)。 `validate` サブコマンドで存在確認 + cgmartini.nl
    からの取得手順を表示。
  - `MembraneCGBuilder` は build 時に REQUIRED ITP を `output_dir/` に
    自動コピー (`gmx grompp` の bare-name `#include` 解決のため、
    `cg.peptide` 流儀踏襲)。
  - **`topology_composer`** モジュール: insane の出力 topology を
    Martini 3 用に post-process。`get_moleculetype_name_from_itp`,
    `compose_topology`, `normalize_ion_atom_names_gro` の 3 関数。
  - **`system_packer.write_ndx_from_gro_cg`**: gmx make_ndx を呼ばず、
    Python だけで `.gro` の residue 名から groups を分類して `.ndx`
    生成 (`POPC` -> Bilayer, 標準 amino acids -> Peptide, `W` -> W,
    `NA`/`NA+` -> NA, `CL`/`CL-` -> CL, それ以外 -> Non_Bilayer)。
  - `tests/test_cg_membrane_*.py` で **166 件の unit test** を整備
    (insane / gmx / martinize2 不要、CI 通せる)。実機 smoke は
    `tests/test_cg_membrane_integration.py` (`@pytest.mark.slow`、4 件)
    で gated -- abmptoolsenv (insane 1.2.0 + vermouth 0.15 +
    GROMACS 2021.3 + AmberTools tleap) で KGG x1 / POPC 32-per-leaflet /
    13 windows × 100 steps の build → `gmx grompp -f em.mdp` /
    `pull.mdp` / `windows/win_006/window.mdp` の **3 stage 全て**
    PASS まで 19.11s で実機検証済。

### External dependencies (new)

- `insane` (PyPI、**GPL-2.0**): Martini bilayer assembly。
  abmptools は subprocess 経由のみ、ソース改変・同梱なし
  (GPL-2.0 接触なし)。`pip install abmptools[cg]` で自動 install。

### Documentation

- **`docs/cg_membrane.md`** (~330 行): サブシステム reference doc。license
  ルール (insane GPL-2.0 + vermouth Apache-2.0 + Martini 3 ITP 未同梱)、
  7-stage pipeline 詳解、API (6 dataclass)、CLI、6 つの設計判断
  (`cg.peptide` sub-call / AA helper の duck-typing 再利用 / NVT-pull /
  NPT-window split / pbcatom windows-only / topology post-process /
  Python ndx writer)、既知 caveat、デフォルト値の根拠。AA `docs/membrane.md`
  の CG 版。
- **`docs/tutorial_cg_membrane_us.md`** (~510 行): step-by-step
  ops tutorial。env 準備 → smoke (5-6 分) → smoke の限界解説 →
  production (45 分) → smoke vs production 比較 plot →
  別 lipid / 別 peptide / 失敗パターン → 参考リンク。
  AA `docs/tutorial_membrane_us.md` の CG 版。
- **`sample/cg_membrane/kgg_popc_production.json`** (新規): 31 windows ×
  5 ns × k=500 × Δz=0.10 の production reproduction config (~45 分)。
  smoke (`kgg_popc_smoke.json`) と並べて配置。
- **Top-level `README.md` + `docs/overview.md`**: cg_membrane / tutorial
  へのリンクを追加。

### Notes

- abmptools の **MO-AAMD-CGMD マルチスケール基盤** の CG 系統 2 番目の
  モジュール。先行 `abmptools.cg.peptide` (1.18.0) を内部で sub-call
  する形で再利用、コード重複ゼロ。AA 系の `abmptools.membrane`
  (1.17.3、CHARMM36 / AMBER Lipid21 backend) と並走する CG 版位置付け。
- 後続予定: 複数 lipid 種混合 (POPC/POPE/CHOL 等)、
  small-molecule permeant (`cg/smallmol/` Auto-Martini 経由) は
  v1.20+ に deferred。

## [1.18.0] - 2026-05-04

新サブパッケージ `abmptools.cg.peptide` を追加。abmptools の MO-AAMD-CGMD
マルチスケール基盤化に向けた CG (粗視化) 系統 (`abmptools/cg/`) の最初の
モジュール。後続で `cg/polymer/` (polyply 経由) や `cg/smallmol/`
(Auto-Martini 経由) を計画している。

### Added

- **`abmptools.cg.peptide` — Martini 3 peptide CG system builder** (新サブパッケージ)
  - JSON/YAML 入力から GROMACS MD 入力一式を生成する end-to-end builder。
    abmptools の MO-AAMD-CGMD マルチスケール基盤化に向けて新設した
    `abmptools/cg/` namespace の最初のモジュール。
  - `PeptideCGBuilder.build()` で 6 stage を 1 呼び出し:
    atomistic PDB (tleap or extended-backbone fallback) →
    `martinize2 -ff martini3001` で M3 CG mapping →
    `gmx insert-molecules` で peptide 配置 →
    `gmx solvate` で Martini W →
    `gmx grompp` + `gmx genion` で NaCl 中和 + 0.15 M →
    em/nvt/npt/md.mdp + index.ndx + run.sh を生成。
  - データクラス: `PeptideSpec` / `PeptideBuildConfig`
    (`@dataclass` + `to_json/from_json`、JSON 往復、YAML は optional)。
  - CLI: `python -m abmptools.cg.peptide {build,validate,example}`
    (argparse、`amorphous` / `membrane` 流儀準拠、Click 不使用)。
  - **Apache-2.0 互換のみ**: vermouth-martinize (Apache-2.0) を `subprocess`
    経由で利用、改変 / 同梱なし。Martini 3 の `.itp` は本パッケージ未同梱で、
    `validate` サブコマンドで存在確認 + cgmartini.nl からの取得手順を表示。
  - ライセンス未明記の cgmartini 配布物 (martinize-dna.py 等) は除外し、
    Martini 3 ペプチド機能のみを切り出した。
  - `tests/test_cg_peptide_*.py` で 104 件の unit test を整備
    (martinize2/gmx/tleap 不要、CI 通せる)。実機 smoke は
    `tests/test_cg_peptide_integration.py` (`@pytest.mark.slow`) で gated --
    abmptoolsenv (vermouth 0.15 + GROMACS 2021.3 + AmberTools tleap) で
    KGG x1 / 4 nm cubic box の build → `gmx grompp -f em.mdp` PASS まで
    実機検証済。
  - 新 extras: `pip install abmptools[cg]` で vermouth + pyyaml が入る。
  - `forcefield_check`: REQUIRED Martini 3 files は ITP 3 ファイル
    (`martini_v3.0.0.itp` / `_solvents_v1.itp` / `_ions_v1.itp`)。
    `martini_v3.0.0_water.gro` は cgmartini.nl が直接配布していないため
    OPTIONAL 扱い (`solvent_enabled=True` のときのみ要)。
  - `PeptideCGBuilder` は build 時に REQUIRED ITP を `output_dir/` に
    自動コピーし、`topol.top` の bare-name `#include` を `gmx grompp` が
    `cwd` から解決できるようにする (amorphous / membrane 流儀)。
  - **`water_box.make_martini_water_box`**: cgmartini.nl が
    `martini_v3.0.0_water.gro` を直接配布していないため、ff_dir に
    water box が無い場合は `gmx insert-molecules` で W bead を
    Martini 標準密度 (~8.36 W/nm^3、5 nm cubic で約 1045 beads) で詰めて
    自動生成する。builder の `_stage4_solvate` で透過的に呼ばれる。
    ユーザーが Martini 3 tutorial archive 等から自分で water.gro を
    用意した場合はそちらが優先される。
  - **実機検証 (abmptoolsenv で 17.42s)**:
    - `solvent_enabled=False`: packed.gro -> `gmx grompp` PASS
    - `solvent_enabled=True`: water_box auto-gen -> `gmx solvate` ->
      `gmx genion` (NaCl 中和 + 0.15 M) -> system_ions.gro ->
      `gmx grompp` PASS

### External dependencies (new, optional via `[cg]` extra)

- `vermouth` (PyPI、Apache-2.0): `martinize2` CLI 提供。
- `gmx` (GROMACS, conda): `solvate` / `genion` / `grompp` / `make_ndx` /
  `insert-molecules`。
- `tleap` (AmberTools, conda; **推奨**): atomistic PDB 生成。**不在時は
  extended-backbone fallback** だが、芳香族残基 (W/F/Y) で CG sidechain
  bead が NaN になる可能性あり。研究用途では tleap 推奨。

### Notes

- Pydantic + Click ベースの非公開 m3-peptide リポジトリ (Martini 2/3 +
  ssDNA 対応) からの再実装。DNA 関連 (`external/martinize-dna.py` 等の
  ライセンス未明 derivative) は除外し、`abmptools.amorphous` /
  `abmptools.membrane` と同じ流儀 (`@dataclass` + `argparse`) に統一。

## [1.17.3] - 2026-05-04

CHARMM36 backend 実機検証で見つかった 7 件の互換性 bug を修正。Klauda lab
GROMACS port (`charmm36-feb2026_cgenff-5.0.ff/` / `charmm36-jul2022.ff/`) で
peptide + bilayer + water + ion の smoke build を end-to-end で完走できるように
なった。

### Fixed

- **PDB 4-char residue 名の切り詰め回避** (`parameterize_charmm.translate_pdb_amber_to_charmm`):
  packmol-memgen `--charmm` 出力は POPC / TIP3 / SOD 等を col 18-21 の 4-char で
  書き込む (chain ID 列にはみ出す)。3-char `line[17:20]` fallback で読むと
  POPC → POP に切り詰められて pdb2gmx が `Residue 'POP' not found` で fatal
  終了していた。`line[17:21].strip()` に統一。
- **Klauda port 規約に合わせた residue rename テーブル縮小** (`AMBER_TO_CHARMM_RESNAME`):
  `NME → CT3`, `CYM → CYS`, `HIP → HSP` を削除。Klauda port は `[ NME ]`
  `[ CYM ]` `[ HIP ]` を直接定義しており、CHARMM オリジナルの慣例に
  rename すると `Residue type 'CT3' not found` で fatal。
- **ACE atom 名 mapping を Klauda 規約に修正** (`PER_RESIDUE_ATOM_MAP['ACE']`):
  Klauda port は AMBER 寄りの命名 (CH3/HH31/HH32/HH33/C/O) を採用しており、
  CHARMM オリジナルの CAY/HY1/HY2/HY3/CY/OY ではない。AMBER tleap 由来の
  H1/H2/H3 を HH31/HH32/HH33 にマップ、CH3/C/O はそのまま通す。
- **NME atom 名 mapping を新規追加** (`PER_RESIDUE_ATOM_MAP['NME']`):
  AMBER NME (`N H C H1 H2 H3`) → Klauda port NME (`N HN CH3 HH31 HH32 HH33`)。
  C → CH3、H1/H2/H3 → HH31/HH32/HH33。`SKIP_BACKBONE_RENAME` から `'CT3'`
  を削除し、NME は universal H → HN を適用。これにより `atom C in residue
  NME 7 was not found in rtp entry NME with 6 atoms` 解決。
- **terminus 'None' index の hardcode** (`KNOWN_CHARMM_TERMINUS_NONE` /
  `_resolve_terminus_none_indices`): pdb2gmx `-ter` のメニュー順は ff の
  `.n.tdb` / `.c.tdb` 読み込み順で決まり、`charmm36-feb2026` /
  `charmm36-jul2022` 系では `None = (8, 7)` (N-term: 8、C-term: 7)。
  デフォルトの `0/0` (`NH3+/COO-`) では ACE/NME cap で
  `atom N not found in buiding block 1ACE` で fatal。subprocess による
  動的 probe は pdb2gmx 2021.3 が無効 input で **99% CPU 無限 spin**
  する未解決 bug を踏むため使わず、ff ごとに手動登録した dict で hardcode。
- **TIP3 spurious O-H-H angles の post-process** (`_strip_water_spurious_angles`):
  CHARMM port `solvent.rtp` の `[ TIP3 ]` には rigid 制約用の H1-H2 "bond"
  が並び、pdb2gmx は通常 bond と解釈して 1 water あたり 3 angles
  (1 real + 2 spurious) を生成。`ffbonded.itp` は `HT-OT-HT` のみ定義
  なので、各 spurious angle が `No default U-B types` error 1 件を出す
  (3500 waters で 7000+ errors)。`run_pdb2gmx` の post-process で全
  chain itp の `[ angles ]` セクションをパースし、middle atom が `OT`
  以外の 3-tuple を削除。除去件数を log 出力。
- **subprocess 自動 probe の禁止 (regression 防止)**: pdb2gmx 2021.3 が
  無効 input (`999\n999\n` 等) で 99% CPU 無限 spin する bug を踏まないよう、
  terminus index 解決から subprocess probe を完全削除し、`KNOWN_CHARMM_TERMINUS_NONE`
  辞書 lookup のみに simplify。新 ff port 追加時はメニューを手動で
  1 度確認して dict に登録する運用。

### Tests

- **`tests/test_membrane_charmm_translate.py` の更新と新規テスト**:
  - `test_histidine_tautomers`: HIP は港 (port) 規約で verbatim、HID/HIE
    のみ HSD/HSE に rename
  - `test_protomers`: ASH/GLH/LYN → ASPP/GLUP/LSN、CYM は verbatim
  - `test_caps`: ACE/NME ともに residue 名は verbatim、atom のみ rename
  - `test_ace_atoms`: H1/H2/H3 → HH31/HH32/HH33; CH3/C/O pass through
  - `test_nme_atoms` (新規): C → CH3、H1/H2/H3 → HH31/HH32/HH33、universal
    H → HN
  - `test_4char_residue_not_truncated` (新規): POPC / TIP3 が POP / TIP に
    切り詰められない regression test
  - `test_skip_list_contains_caps_and_ions`: NME が SKIP_BACKBONE_RENAME に
    含まれていないことを assert
  - `test_full_translation_pipeline`: ACE 6 atoms (H1/H2/H3/CH3/C/O) を含む
    fixture に拡張、HH31/HH32/HH33 すべての翻訳を verify
  - 全 26 件 PASS / 0.31 s

### Documentation

- 業界実態を `memory/` に文書化 (`project_membrane_charmm36_gotchas.md`):
  - ~95% の academic は CHARMM-GUI を使い pdb2gmx を完全回避
  - 少数派は psfgen + charmm2gmx (Wacha & Lemkul 2023 *JCIM*)
  - 完全 commercial-OK の membrane builder は packmol-memgen のみ
  - **本パッケージの修正は世界で数十人しか踏まない pain points**
- AMBER TIP3P と CHARMM-modified TIP3P の H atom LJ 違い (σ=0/ε=0 vs
  σ≈0.4/ε≈0.046 kcal/mol) を memory に記載。abmptools.membrane では
  backend ごとに自動で適切な TIP3P が選ばれる (AMBER: tleap、
  CHARMM36: pdb2gmx の ff dir tip3p.itp)。

### Known limitations

- 実機 MD 完走 (production-leaning run, ~3 時間 GPU) はまだ未検証。
  smoke (grompp pass) までは PASS。
- 長期的には "正攻法 (a)": protein のみ pdb2gmx + water/lipid/ion は
  ff dir の itp を直接 `#include` する方式へ refactor 余地あり (現状は
  業界正攻法の (b) post-process 方式)。

## [1.17.2] - 2026-05-03

### Added

- **`abmptools.membrane` の `DEFAULT_LIPID_APL` を 14 → 60 entries に拡張**。
  Lipid21 の標準 lipid を網羅 (PC: 11, PE: 10, PG: 10, PS: 10, PA: 9,
  SM: 7, sterol: 3)。命名は packmol-memgen / Lipid21 規則に従う:
  - **PC**: DLPC/DMPC/DPPC/DSPC/POPC/PMPC/SOPC/DOPC/DAPC/DHPC/AHPC
  - **PE**: DLPE/DMPE/DPPE/DSPE/POPE/PMPE/SOPE/DOPE/DAPE/DHPE
  - **PG**: DLPG/DMPG/DPPG/DSPG/POPG/SOPG/DOPG/DAPG/DHPG/AHPG
  - **PS**: DLPS/DMPS/DPPS/DSPS/POPS/SOPS/DOPS/DAPS/DHPS/AHPS
  - **PA**: DLPA/DMPA/DPPA/DSPA/POPA/DOPA/DAPA/DHPA/AHPA
  - **SM** (sphingomyelin、raft component、Lipid21 only):
    LSM/MSM/PSM/SSM/OSM/ASM/HSM
  - **sterol**: CHL1 / CHOL / CHL (alias)
  - APL 値は文献ベース (Kučerka 2011 *BBA*, Marsh 2013 *BBA*, Lipid21
    MD literature) の Lα 相 ~310 K 平均
- **検索 / ロード helper** (`bilayer.py`):
  - `list_known_lipids(head_group=None)` — curated table を head 別 sort
    + filter ("PC"/"PE"/"PG"/"PS"/"PA"/"SM"/"sterol")
  - `_classify_lipid_head(resname)` — resname の suffix / prefix から
    head group 推定
  - `query_packmol_memgen_lipids(packmol_memgen_path)` — packmol-memgen
    の全 259 種を runtime で取得 (curated table 外も含む)
- **CLI**: `python -m abmptools.membrane.lipid_info`
  - `--known` (default) / `--known --head SM`: curated table 表示
  - `--available`: packmol-memgen 全 259 種
  - `--apl RESNAME`: 単一 lipid の APL 解決 (table miss → 65.0 fallback 表示)
- **23 ユニットテスト追加** (`tests/test_membrane_mixed_lipid.py`):
  table size sanity (≥50)、head group 完全性、`_classify_lipid_head`
  の各 head + sterol + other ケース、`list_known_lipids` の filter /
  sort / value 整合性。`v1.17.1` の 22 テストと合わせて 46 件 / 0.23 s

### Documentation

- `docs/membrane.md` の `LipidSpec` セクションを書き換え:
  - 旧 14-row APL 表 → head 別の概観表 (各レンジで代表 5-10 entries) +
    "完全表は CLI で表示" の案内
  - "既知 lipid のロード / 検索" 新サブセクション (CLI 使用例 + Python
    API 使用例)

## [1.17.1] - 2026-05-03

### Added

- **`abmptools.membrane` で混合脂質をサポート**。`MembraneConfig.lipids` に
  複数の `LipidSpec` を並べるだけで、packmol-memgen に
  `--lipids POPC:CHL1 --ratio 4:1` のような mole-ratio (gcd 約分) と
  `--distxy_fix sqrt(sum(n × APL))` の per-lipid 面積総和に基づく box
  サイズが自動で渡されるようになった。
  - `LipidSpec` に `apl_angstrom2: float = 0.0` フィールド追加。0.0 なら
    `bilayer.DEFAULT_LIPID_APL` テーブルから自動 lookup
    (POPC=67 / DOPC=72 / DPPC=63 / POPE=56 / CHL1=38 等、共通 14 種)、
    explicit に値を指定すれば override (低温 gel 相での DPPC=49 等)。
  - `bilayer.estimate_distxy_angstrom` を per-lipid 計算に変更。引数
    `apl_angstrom2: float = 65.0` は **未知残基への fallback** として
    残るが、`DEFAULT_LIPID_APL` に載っている脂質では参照されない。
  - 22 ユニットテスト追加 (`tests/test_membrane_mixed_lipid.py`):
    table coverage、`_resolve_apl` precedence (explicit > table >
    fallback)、binary / ternary mixture の `estimate_distxy_angstrom`、
    multi-lipid `assemble_packmol_memgen_cmd` の `--lipids` /
    `--ratio` / `--distxy_fix` 出力。
  - 実機 build smoke (POPC 24 + CHL1 6 / leaflet、4:1) で
    52 POPC + 12 CHL + ~3000 water + 16 ions の bilayer が
    packmol-memgen で作成できることを確認。

### Documentation

- `docs/membrane.md` の `LipidSpec` セクションに混合脂質の例 (binary /
  ternary) と `DEFAULT_LIPID_APL` の脂質→APL 対応表を追加。
- `docs/tutorial_membrane_us.md` §2.1 の config 例に混合脂質テンプレートを
  コメントブロックで追加。

## [1.17.0] - 2026-05-03

### Added

- **`abmptools.membrane` — peptide-bilayer umbrella-sampling builder** (新サブパッケージ、Phase A〜D)。
  ペプチドの脂質膜透過 PMF 計算用の GROMACS 入力一式を生成する end-to-end ビルダー。
  - `MembraneUSBuilder.build()` で 6 stage を 1 呼び出し:
    bilayer 構築 (packmol-memgen) → AMBER パラメータ化 (tleap + parmed) → 平衡化 MDP
    (em / nvt / npt-semiisotropic、2-group thermostat) → 反応座標生成 pulling
    MDP → US window MDP 一括 → top-level `run.sh`
  - 生成された `run.sh` で `gmx grompp + mdrun` をシーケンシャル実行し、最終的に
    `gmx wham` で `analysis/pmf.xvg` (PMF[z]、kJ/mol) を出力
  - データクラス: `MembraneConfig` / `LipidSpec` / `PeptideSpec` / `IonSpec` /
    `USProtocol` / `EquilibrationProtocol` / `PullingProtocol` (JSON 往復可)
  - **商用利用可な力場のみ**: AMBER ff19SB + Lipid21 + TIP3P + Joung-Cheatham ions
    (AmberTools 配布、free incl. commercial)。CGenFF Web server / CHARMM-GUI に
    依存しない設計
  - **CHARMM36 backend** (Klauda lab GROMACS port、CGenFF 不使用) を
    `parameterize_charmm.py` で実装済 (Phase C):
    - `--charmm` フラグで packmol-memgen に lipid/water/ion を CHARMM 命名で
      出力させ、peptide は `translate_pdb_amber_to_charmm` で AMBER→CHARMM
      残基/原子名翻訳 (HIE→HSE / NME→CT3 / 末端 cap atom names / 標準 AA の
      amide H→HN 等) を当ててから `gmx pdb2gmx -ff charmm36-jul2022 -water
      tip3p` で top/gro を生成
    - **ion auto-translate**: `IonSpec(cation='Na+', anion='Cl-')` を AMBER
      形式で書けば CHARMM backend が内部で `SOD`/`CLA` 等にマップ
    - 25 ユニットテスト (`tests/test_membrane_charmm_translate.py`) で
      翻訳ロジックを検証
    - Klauda 研の `charmm36-jul2022.ff` 配布を `MembraneConfig.charmm_ff_dir`
      で参照する形 (本パッケージは未同梱、license 配布元の差し替えを許容)。
      取得手順は `docs/membrane.md` の "CHARMM36 GROMACS port の取得" 参照
  - smoke test: `tests/integration/run_membrane_us_smoke.sh`
    (poly-Ala 5-mer + POPC 32/leaflet + 7 windows × 1 ns)。
    16 ファイル生成 + 11 MDP すべて `gmx grompp` 通過を ~30 秒で検証
  - GPU 加速対応: `MDRUN_OPTS` env hook を `run.sh` に追加。NVIDIA + CUDA
    で 18k atom 系 ~640 ns/day (CPU 4-core 比 ~4-5×)。WSL2 環境では
    side-env (`gmxcudaenv`) パターンを `docs/tutorial_membrane_us.md` で案内
  - 詳細: [`docs/membrane.md`](docs/membrane.md) (reference) /
    [`docs/tutorial_membrane_us.md`](docs/tutorial_membrane_us.md) (step-by-step ops)

### Notes

- **packmol-memgen 2023.2.24 + NumPy ≥ 1.24 互換性パッチ**: bundle 内
  `pdbremix/v3numpy.py` が削除済 `np.float` を参照しているため、env 内で
  `np.float` → `float` の 2 行 sed パッチが必要。詳細は
  [`docs/membrane.md`](docs/membrane.md) のインストールセクション

### Documentation

- 横断ドキュメント監査の結果を反映 (commit `60faa43`):
  - `README.md` に membrane の Features / Quick Start / Documentation 追加
  - `docs/dependencies.md` / `docs/licenses_third_party.md` に
    `amorphous` / `membrane` 専用セクション新設
  - `docs/overview.md` / `docs/architecture.md` / `docs/directory_structure.md`
    の subpackage / 機能リストを membrane に対応
  - `pyproject.toml` に `membrane` extras を追加 (`parmed` 等の pip 依存のみ;
    gromacs / ambertools は conda)
  - その他陳腐化表記 (Phase 計画形 → 実装済) を解消

## [1.16.0] - 2026-05-01

(release commit lives on `develop` only; not yet merged to `main` /
tagged / uploaded to PyPI. Per the user, the next release window
will pick this up after additional integration testing.)

### Added
- `abmptools.core.system_model` に COGNAC 固有情報を保持するデータクラスを追加:
  - `ClusterData` (cluster 配置 xyz / n_per_cluster / cluster_file)
  - `FixedLabel` (固定原子 atom_indices / label)
  - `SystemModel.cluster_data` / `fixed_labels` / `ensemble_family` フィールド
  - 判定ヘルパー `classify_ensemble(algorithm)` と `COGNAC_ONLY_ALGOS` 定数
    (既定: `NPT_Andersen_Kremer_Grest` / `NPT_Andersen_Nose_Hoover`)
- `abmptools.amorphous.system_model_adapter.from_interchange(interchange, ...)` を新設:
  OpenFF Interchange → (一時 .gro/.top 経由) TopModel → SystemModel の最小充填。
  `ensemble_family="gromacs_ok"` を付与して返す。`mol_topologies` は意図的に空のまま
  (GROMACS .top 出力は `interchange.to_top()` を直接使う想定)
- テスト 23 件追加 (`tests/test_system_model_extensions.py` / `tests/test_interchange_adapter.py`、
  integration は `@pytest.mark.slow`)
- **`abmptools.amorphous` で stacked force fields (water FF override)** (`9e4ee46`)
  - `BuildConfig.forcefield: Any` (旧 `str`)。`str` 単一 FF または `list[str] | tuple[str]`
    で SMIRKS-overlay (後ろが前を上書き) として `OpenFF ForceField(*names)` に展開
  - `create_interchange.forcefield_name` 引数も同様に `str | Sequence[str]` 受付
  - 典型用途: `['openff_unconstrained-2.1.0.offxml', 'tip3p.offxml']` で organic は
    OpenFF organic、water 分子のみ TIP3P。GAFF/openff-water の repulsive σ/ε で
    純 water 系が 1.0 → 0.26 g/cm³ に膨張する問題への根本対策
  - 実機検証: pure water (200 TIP3P) anneal stage で density=0.991 g/cm³
    (literature ~0.99 g/cm³ at 298 K) を 2026-04-29 確認
  - tests/test_parameterizer.py (新規、5 件): str / list / tuple / 空 list / default 単一 FF
  - tests/test_amorphous_models.py: BuildConfig roundtrip with list[str] forcefield (3 件追加)
- **`ComponentSpec.pdb_path`** (Phase 9-a、2026-04-30): SMILES / SDF に加えて
  pre-built oligomer PDB を直接 component 入力にできる経路を追加。
  `smiles` / `sdf_path` / `pdb_path` は exactly one を `__post_init__` で
  validation。`prepare_molecule` は pdb_path 経路で
  `Molecule.from_polymer_pdb` → `Molecule.from_file` fallback で OpenFF
  Molecule を構築 (simple sp3/sp2 oligomer がターゲット、複雑 polymer
  は FF assignment が fail する可能性あり)。fcews-manybody 側の
  `_setup_amorphous_gromacs` から polymer/pdb/<label>.pdb を流す経路で利用
  - tests/test_amorphous_models.py: pdb_path / 二重指定エラー / 全 empty
    エラーの 3 件追加

### Changed
- `abmptools.udf2gro.gromacs.writers` の GroWriter / TopWriter / MdpWriter / ItpWriter で
  `ensemble_family == "cognac_only"` を検出した場合に `ValueError` を送出
  (共有ヘルパー `_validator.raise_if_cognac_only`)。COGNAC 固有アンサンブル
  (`NPT_Andersen_Kremer_Grest` 等) を GROMACS 形式で誤って書き出すことを防ぐ

### Fixed
- **pure-component pair (同名 component 2 個) で発生していた 3 件の runtime bug** (`d385fb3`)
  - `builder._tc_grps_string()`: component name が dedup されず `tc-grps = B_water B_water`
    (2 group) を出していた。`mdp_protocol` は `ref-t` / `tau-t` を 1 値しか書かないため
    grompp が `Invalid T coupling input: 2 groups, 1 ref-t values` で fatal。
    first-seen-order を保つ dict dedup で修正
  - `ndx_writer.write_ndx()`: `groups[comp_name] = comp_indices` という上書き書きで、
    pure-pair の 2 周目が 1 周目を silently overwrite していた。例: 200 H2O が
    `[ B_water ]` group には atom 301..600 しか入らず、grompp の `tc-grps`
    と atom 数で整合しない。aggregate-by-name (`groups[name].extend(...)`) に修正
  - `mdp_protocol.write_run_script()`: 既定で `gmx mdrun` (= 3 tMPI ranks) を発行。
    凝縮系 (water 等) は NPT-high の 600 K で気化 → anneal 冷却で液相に再凝集 →
    box が ~1.8 nm まで縮み、`box size in direction X is too small for a cut-off
    of 1.214 nm with 3 domain decomposition cells` で停止。`MDRUN_OPTS="${MDRUN_OPTS:--ntmpi 1}"`
    を default に変更し DD 自体を無効化 (OpenMP は引き続き効くので速度低下なし)。
    ユーザー側で `MDRUN_OPTS=...` で override 可能
- **mixed-component pair の thermostat / annealing schedule cardinality 不整合** (`c443b12` + `a66cf6e`)
  - `mdp_protocol._thermostat_block`: 多成分系 (例: `tc-grps = A_methanol B_water`、
    2 group) でも `ref-t` / `tau-t` が単一値で書かれており grompp が
    `Invalid T coupling input: 2 groups, 1 ref-t values and 1 tau-t values` で
    fatal していた。`tc-grps.split()` の token 数に合わせて scalar を複製。
    Pure-pair (1 group) は従来通り単一値出力で挙動不変
  - `mdp_protocol.generate_anneal_mdp`: `annealing` / `annealing-npoints` /
    `annealing-time` / `annealing-temp` が GROMACS の per-group 仕様に従わず
    1 group 想定で書かれていた。多成分系では stage 4 (anneal) の grompp で
    `Inconsistent number of components in annealing-time and annealing-temp`
    が出て stop。これも `n_groups` 倍に複製。1 group の場合は token 単一で
    legacy 互換

## [1.15.4] - 2026-04-19
### Added
- `tests/test_builder_mocked.py` (8 tests): `AmorphousBuilder.build()` の 6 stage フローと返り値 dict 構造 (`wrap_script` キー含む)、`config.json` 書き出し、MDP/ndx 生成順序を mock ベースで検証。OpenFF/Packmol/Interchange なしで CI 可
- `tests/test_builder_integration.py` (12 tests, `@pytest.mark.slow`): methane ×10 / box 2 nm の小系で `AmorphousBuilder.build()` を実際に走らせ、成果物 (gro/top/ndx/5 MDP/2 scripts/config.json) を spot check。OpenFF + Packmol + AmberTools + RDKit が揃ったときのみ実行、足りない依存は `pytest.importorskip` / `shutil.which` で自動 skip

### Changed
- `pyproject.toml`: `[tool.pytest.ini_options]` に `markers = ["slow: ..."]` を登録 (integration 系のゲート用)

## [1.15.3] - 2026-04-19
### Added
- `abmptools.amorphous.pubchem` モジュール: PubChem PUG REST API ラッパー
  - `fetch_3d_sdf(query, by)` / `fetch_smiles(query, by)` / `download_3d_sdf(query, path, by)`
  - `by` は `cid`, `name`, `smiles`, `inchi`, `inchikey`
  - 3D conformer が無い場合は `PubChemNo3DError` を明示的に送出
  - `urllib` 標準ライブラリのみ使用 (追加依存なし)
  - CLI: `python -m abmptools.amorphous.pubchem --cid 3825 -o out.sdf` / `--name aspirin --smiles-only`
- `build_amorphous.py` に `--pubchem_cid` / `--pubchem_name` / `--pubchem_cache_dir` オプションを追加
  - 指定された CID/名前から 3D SDF を取得し、そのまま `--mol` 入力として扱う
  - ダウンロード済み SDF はデフォルトで `<output_dir>/input/` にキャッシュ
- `tests/test_pubchem.py`: 11 テスト (network をモックした HTTP 挙動検証)

### Changed
- README.md: Amorphous Features / Quick Start に PubChem 入力 (`--pubchem_cid`) の記述を追加
- docs/:
  - `amorphous.md`: PubChem 自動ダウンロード節を新設、CLI クイックスタートに `--pubchem_cid` 例を追加
  - `dependencies.md`: amorphous セクション末尾に PubChem 追加依存の記述 (urllib 標準のみ、ただし `pubchem.ncbi.nlm.nih.gov` への HTTPS アクセス必須)
  - `faq.md`: "Can I fetch 3D SDFs automatically from PubChem?" を追加
  - `overview.md`: Amorphous Builder 行に PubChem CID 入力を併記
  - `architecture.md`: Subpackages の amorphous 説明に `amorphous.pubchem` の位置づけを追記
  - `dataflow.md`: Amorphous Build Pipeline 図の入力部に PubChem 分岐を追加
  - `ABMPTools-user-manual.md`: amorphous 節に PubChem 入力対応を追記

## [1.15.2] - 2026-04-19
### Added
- amorphous: 自動生成される `md/wrap_pbc.sh` (`gmx trjconv -pbc mol -ur compact` を各 xtc / 最終 gro に適用、VMD で開きやすい `*_pbc.xtc` を生成)
- サンプル `sample/amorphous/ketoprofen_pubchem/`: PubChem 3D SDF (CID 3825, MMFF94 最適化済、水素込み) を `--mol` で読み込む SDF 入力のサンプル一式 (README, run_sample.sh, input SDF 同梱)

### Changed
- README.md: amorphous 機能セクションを拡充 (SMILES/SDF 両対応、Packmol + OpenFF + AM1-BCC、5-stage annealing、`wrap_pbc.sh` 言及)、Quick Start に SDF 入力例を追加、Samples セクションに amorphous サンプル (pentane_benzene / ketoprofen_pubchem) への導線を追加
- docs/:
  - `amorphous.md`: 出力ファイル一覧に `wrap_pbc.sh` / `*_pbc.xtc` / `05_npt_final_pbc.gro` を追加、ビルド後のワークフロー (run_all.sh → wrap_pbc.sh → VMD) セクションと同梱サンプル一覧を新設
  - `dependencies.md`: `abmptools.amorphous` 専用の Optional Dependencies セクションを追加 (必須ランタイム + 電荷バックエンド + `setuptools<81` 注記 + 後処理外部ツール)、Dependency Summary ツリーにも amorphous 行を追加
  - `faq.md`: amorphous ビルダー向けトラブルシュート 3 件を追加 (Packmol `Illegal seek`、`pkg_resources` 消失、WSL2 の NVIDIA OpenCL ICD 不在)
  - `overview.md`: Key Capabilities 表に `Structure Optimization` と `Amorphous Builder`、MD Integration に gro2udf/udf2gro を追記
  - `architecture.md`: サブパッケージ (gro2udf/udf2gro/geomopt/amorphous/core) を紹介する Subpackages セクションを追加
  - `dataflow.md`: Amorphous Build Pipeline の ASCII フロー図を追加
  - `ABMPTools-user-manual.md`: Overview 部に gro2udf/udf2gro、geomopt、amorphous の概要を追加

## [1.15.1] - 2026-04-18
### Fixed
- amorphous/packing.py: packmol 21.2.1 (conda-forge) の stdin シークエラー対応 (`stdin=open(inp_path, "rb")` + pdb/output の絶対パス化)
- anlfmo: `Pool` import の復旧と絞り込み過剰だった `except` 句の修正

### Added
- リグレッションテスト (`tests/test_regression.py`): リファクタリング前の参照出力との比較で挙動ドリフトを検出
  - 51 bundled (`tests/regression/reference/prerefactor/` 同梱) + 9 sample-based (`sample/` 配下の参照との比較) + 16 gated (外部 `abmptools-sample` 依存のため通常は skip)
  - 対象ツール: generateajf, log2cpf, convertcpf, udf2gro, gro2udf, getifiepieda
- `tips/cp_for_dist.sh` の配布物更新

### Changed
- README.md: 回帰テストの説明および developer-only gated テストのセクションを追記、インストール手順を editable (`pip install -e .`) 推奨に変更、Quick Start に amorphous 使用例を追加、テスト件数を 671/30 に更新
- CHANGELOG.md: リリースバージョン毎に整理し直し、旧リリース日付のゼロパディングを統一

## [1.15.0] - 2026-03-21
### Added
- udf2gro サブパッケージ: OCTA UDF → GROMACS (gro/top/mdp/itp) 変換機能
- gro2udf サブパッケージ: GROMACS → OCTA UDF 変換機能
  - `--from-top` モード (topファイルからの変換、NH-Q・Ewald・デフォルトテンプレート対応)
- geomopt サブパッケージ: 構造最適化機能
  - MacePdbOptimizer (MACE/ASE ベースのPDB構造最適化)
  - OpenFFOpenMMMinimizer (OpenFF力場によるPDB構造最小化)
  - QMOptimizerPySCF (PySCF量子化学計算による構造最適化)
- amorphous サブパッケージ: 多成分アモルファス系構造構築機能 (packmol/OpenMM)
- core サブパッケージ: SystemModel 共通データモデル
- 開発者向けドキュメント9件 (architecture, dataflow, dependencies, io_spec, faq 等)
- pytest テストスイート: 28ファイル (全モジュール + 全14 CLIスクリプト、約620テスト)
- Japanese Google-style docstrings: 全公開メソッド/クラス/モジュールに追加
- 型ヒント: abinit_io, anlfmo, pdb_io, readcif 等 (89メソッド)
- CLIスクリプト用 `get_args()` 関数の抽出 (8スクリプト)
- `pyproject.toml` 追加: PyPI publishing 対応 + 全 CLI スクリプトに `main()` エントリポイント

### Changed
- `print()` → `logging` モジュールへ置換 (コアモジュール + gro2udf/udf2gro)
- `exec()`/`eval()` → 安全なデータ読み込みに置換
- bare `open()` → `with` 文コンテキストマネージャに変換
- `subprocess.call(mkdir)` → `os.makedirs()` に置換
- `try/except KeyError` → `dict.get()` パターンに置換
- `is True`/`is False` パターンの修正
- 未使用 import の削除 (9モジュール)
- `doc/` → `docs/` ディレクトリ名変更
- ドキュメント更新: TEST_COVERAGE.md, dev_quickstart.md, directory_structure.md
- README.md を英語に書き換え (全機能を網羅した形に)
- TEST_COVERAGE.md に未テスト関数インベントリを追加
- `.gitignore` に egg-info と `__pycache__` を追加

### Fixed
- cooperative inheritance chain の修復 (MRO関連)
- icflag バグの修正
- エスケープシーケンスの修正
- setup.py: gro2udf/default_template.udf をパッケージデータに含める

## [1.14.6] - 2025-12-21
### Fixed
- log2config (logmanager) の不具合修正
  - 核酸への対応
  - 核酸/タンパク複合体の計算
  - V2 Rev.8 で、CYS架橋がある際にテーブルがずれる例外処理に対応
- fcewsのoutファイル数check機能の修正

## [1.14.5] - 2025-09-22
### Added
- log2configモジュール(nprint=0のログから、fragment configファイルに変換する機能)を追加
- generateajfモジュール: configファイルからajfを生成する機能を追加

## [1.14.1] - 2024-05-17
### Added
- cpf2ifielistモジュール(cpfを読み込んで, 整形されたIFIEリストを出力する機能)を追加
### Fixed
- cpfmanager: CPF Ver.10において、bda-baa 原子が5桁を超えた際の読み込みエラーを修正

## [1.14.0] - 2024-05-12
### Fixed
- [Manualの加筆](doc/ABMPTools-user-manual.md)
- getifiepieda Pandas2系に対応するように修正(append)
- ABINIT-MP Ver.2 Rev.8 対応の一部不具合修正(&CIS等)

## [1.13.5] - 2024-03-14
### Fixed
- 13.4のエラー対応時の二重読み込みエラーの修正

## [1.13.4] - 2024-02-16
### Fixed
- FMOPB 特定のエラー終了時の読み込みエラー対応

## [1.13.3] - 2024-02-07
### Changed
- DIFIE (DIFIE) 出力仕様変更に伴う出力変更

### Fixed
- UDF関連一部修正
- Openbabel対応一部修正

## [1.13.2] - 2023-10-15
### Added
- CPFmanager "CPF ver7.0 (MIZUHO)" 版への対応を追加

## [1.13.1] - 2023-09-27
### Added
- DIFIE 出力機能の並列処理機能 (-np)
- DIFIE のサンプル追加

## [1.13.0] - 2023-09-18
### Added
- Logparser機能(LogManager)
- Logからcpfを作成する機能(log2cpf)
- log2cpfのサンプルを追加

## [1.12.4] - 2023-09-09
### Added
- CPFparser機能(CPFManager)
- DIFIE CPF出力機能(abmptools.generate_difie)
- 一部機能のサンプルを追加

## [1.12.3] - 2023-09-08
### Added
- ABINIT-MP Ver.2 Rev.8 対応('v2rev8')

## [1.12.2] - 2023-08-21
### Fixed
- v1最新版でFMO='OFF'でNFの記載があるとエラーが出るとに合わせたajf出力対応

## [1.12.1] - 2023-07-01
### Added
- 会合体フラグメント情報割り当て機能(pdb2fmo, udf2fmo) のsp2分割への対応

## [1.12.0] - 2023-04-02
### Added
- CHANGELOG.md の追加 (change historyの整理)
- バージョンタグ追加

### future plan
- ~~機能ごとテスト作成~~ → 完了 (658テスト)
- ~~リファクタリング~~ → 完了 (refactor/all ブランチ)
- 権利フリーのサンプルの一般公開
- bsse読み込み機能(developのみ)の影響チェック

## [1.11.3] - 2023-03-15
### Added
- 電荷取得機能
- ログからのオプション取得
- abinitmpテスト一括出力tips
- MD後のtrajectoryを間引くtips
- namd MD後の整形・解析スクリプト改良(rmsd, rdf, dist, autoimage)

## [1.11.0] - 2022-08-25
### Added
- readcifの対応対称性の追加
- gromacs系tips追加
- md-fmo関連tips追加

### Fixed
- getifiepiedaにおけるDimerEnergy取得機能の修正
- N:1出力のバグ修正
- fragidsモードのid selectionの修正

## [1.10.0] - 2022-04-18
### Added
- ABINIT-MP Ver.2 Rev.4対応追加
- ライセンス追加
- autoフラグメント分割の読み取りに対応
- monomer energy, dimer energyの読み取り
- generateajfへのlrdの対応を追加
- 特定の原子の距離情報を測る機能
- md tips群追加
- udfのcutmodeについて、pdbと同一の機能を追加

### Changed
- スクリプト群のモジュール化(-m実行)

## [1.9.0] - 2021-08-11
### Added
- ajfのmldatオプションへの対応
- getifiepiedaへのダイマーES内側のみ取得オプションの追加
- cif構造のpdb, xyz化機能(readcif)の追加
- amber, gromacs, namdのmdtipsの追加

## [1.8.0] - 2021-07-26
### Added
- 連番ajf作成機能
- 簡易マニュアル作成
- rdf算出機能の改良
- 富岳実行スクリプトサンプルの追加

### Changed
cpfのデフォルトを10に変更

## [1.7.0] - 2021-02-13
### Added
- ffmatrix + PB での取得
- getifiepiedaのLRD, HFへの対応
- ajf生成機能のへのcpfバージョン指定機能追加
- ffmatrixで、重複するフラグメントを指定できるように変更

## [1.6.0] - 2020-10-19
### Added
- generateajfを dgemm, mp3, resp, ligand charge に対応するように更新
- ABINIT-MP Open.1 Rev.23に対応
- ifie取得機能のPBへの対応
- pdb2fmoでatomname, residueidの更新を選択できるように拡張

### Changed
- generateajf を argparseによる引数実行に変更


## [1.5.0] - 2020-08-17
### Added
- ajfへのRev22対応, OFP対応
- log or pdbからのフラグメント情報の自動取得
- MP2.5, MP4, MP3.5の情報出力
- 出力csvへの残基情報の書き出し
- READMEへの機能追記 

### Change
- ifie取得機能のgetifiepieda.pyへの集約
- argparseでの引数指定（getifiepieda, pdb2fmo,udf2fmo)

## [1.4.0] - 2020-06-02
### Added
- frag vs fragのmatrix生成
- 時間-フラグメントのマトリックス
- MP3の結果を取得する機能
- svdの実施スクリプト
- ajf生成機能簡易版

## [1.3.0] - 2020-05-18
### Added
- AMBERのMD後構造をcpptrajで解析するスクリプト群

## [1.2.0] - 2020-05-11
### Added
- ifie解析の新機能
    - フラグメントからの距離
    - 分子内の特定のフラグメント番号
    - 分子単位でのifie取得
    - 時系列かつ距離フィルター
    - フラグメントid対を指定した時系列
    - 特定のフラグメントと分子名でのフィルター

### Changed
内部ディレクトリ構成変更

## [1.1.0] - 2020-04-28
### Added
- ABINIT-MPのプリポストツールとして機能統合
    - ajf生成機能
    - ifie取得機能(getifiepieda) 気相中のifie,piedaとMD-FMO用での時間軸での結果取得

## [1.0.8] - 2020-03-05
### Added
- 残基間の動径分布を出す機能
- 分子の並進処理の際、分子idを指定するモードと、既存座標を指定するモード2種を追加

### Fixed
- ajf出力の際の基底関数, cpf名称の軽微なバグ修正


## [1.0.7] - 2019-09-26
### Added
- 周期境界を考慮しpdbの座標を並進させる機能
- PDBの残基間から簡易的な距離分布を出す機能

## [1.0.6] - 2019-08-22
### Added
- 中心座標を出力する機能

### Fixed
- pdb原子上限の対応追加

## [1.0.5] - 2019-08-06
### Added
- pdb読み込みにモードを追加 (rfile, resnum)
- pdbの残基名を変更するtips機能を追加

### Fixed
- pdb読み込みの際の固定長認識を修正

## [1.0.4] - 2019-07-30
### Added
- 切り取らずにフラグメントを割り当てる機能(noneモード)
- 溶質からの距離で切り出す機能(aroundモード)

## [1.0.3] - 2019-07-25
### Changed
- 機能をclass化(pdb_io)

## [1.0.2] - 2019-07-23
### Added
- PDBを読み込む機能を追加

## [1.0.1] - 2019-01-19
### Added
- UDFの任意のレコードをpdbに変換して書き出す機能

## [1.0.0] - 2018-11-16
### Added
- OCTA COGNACのUDFファイルから指定した範囲を切り出す機能
- 切り出した構造にFMOフラグメント情報を割り当てる機能

