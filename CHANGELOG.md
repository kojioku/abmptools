# Changelog

## [1.12.0] - 2023-04-02
### Added
- CHANGELOG.md の追加 (change historyの整理)
- バージョンタグ追加

### future plan
- 機能ごとテスト作成
- リファクタリング
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

