# Changelog

## [Unreleased]

(no changes yet — this section accumulates work-in-progress between releases)

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

