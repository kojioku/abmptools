# Changelog

## [Unreleased]
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

