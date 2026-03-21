# abmptools テストカバレッジマップ

## 概要

| カテゴリ | モジュール数 | テストファイル | テスト数 | ステータス |
|---------|------------|--------------|---------|-----------|
| コアクラス (継承チェーン) | 7 | 7 | 300 | 完了 |
| 複合クラス | 2 | 2 | 72 | 完了 |
| マネージャクラス | 2 | 2 | 64 | 完了 |
| スタンドアロン関数モジュール | 2 | 2 | 61 | 完了 |
| サブパッケージ (core, gro2udf, udf2gro, amorphous, geomopt) | 14 | 14 | 127 | 完了 |
| CLIスクリプト (argparseテスト) | 14 | 1 | 50 | 完了 |
| **合計** | | **28** | **658** | **全PASS** |

---

## テストファイル一覧

### 1. コアクラス (継承チェーン: molcalc → mol_io → abinit_io → pdb_io → anlfmo)

| # | テストファイル | 対象モジュール | 行数 | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|------|---------|--------------|-----------|
| 1 | `test_molcalc.py` | `molcalc.py` | 1170 | 55 | `getdist`, `getangle`, `getcrossprod`, `getCenter`, `getoriginpos`, `gettranspos`, `getmolradius`, `calccellsize`, `moveMolTrans`, `rotate_ardz`, `rotate_ardvec`, `get_element_name`, `get_atomic_radius`, `wrap_to_primary_cell`, `shift_molecule_to_primary_cell`, `getindex`, `getatomisite`, `getmolmass`, `getrenumindex` | [x] 完了 |
| 2 | `test_mol_io.py` | `mol_io.py` | 175 | 27 | `read_mol_name`, `read_xyz`, `getatoms`, `getatom`, `convert_xyz_pdb`, `Exportpospdb` | [x] 完了 |
| 3 | `test_abinit_io.py` | `abinit_io.py` | 2758 | 76 | `__init__` defaults, `chkdepth`, `get_fragsection`, `config_read`, `gen_ajf_body`, `read_ifie`, `getmo_or_fmo`, `flatten`, `functor`, `getfraginfo` | [x] 完了 |
| 4 | `test_pdb_io.py` | `pdb_io.py` | 1875 | 50 | `__init__` defaults, `readpdb`, `getpdbcell`, `exportardpdbfull`, `exportardpdb`, `getpdbinfowrap`, `movemoltranspdb`, roundtrip | [x] 完了 |
| 5 | `test_udf_io.py` | `udf_io.py` | 696 | 7 | `getposatom`, `getposmol`, `getposmolrec`, `getnameAtom`, `getAtomtypename`, `putPositionsMol` (UDFManager mock) | [x] 完了 |
| 6 | `test_udfrm_io.py` | `udfrm_io.py` | 125 | 11 | `__init__`, `getmolname`, `moveintocell_rec`, `convert_udf_pdb` (UDFManager mock) | [x] 完了 |
| 7 | `test_udfcreate.py` | `udfcreate.py` | 1130 | 10 | `__init__`, `setudfparam`, `getconnectdata`, `getbatdata` | [x] 完了 |

### 2. 複合クラス

| # | テストファイル | 対象モジュール | 行数 | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|------|---------|--------------|-----------|
| 8 | `test_setfmo.py` | `setfmo.py` | 637 | 20 | `__init__` defaults, `setrfmoparam` (全キー/部分/デフォルト) | [x] 完了 |
| 9 | `test_anlfmo.py` | `anlfmo.py` | 4410 | 51 | `__init__` defaults, `depth`, `getisdisp`, `getlogmethod`, `getpbflag`, `getifiedf`, `getpiedadf`, `getmomenedf`, `getdimenedf`, `gettgtdf_ff/ffs` | [x] 完了 |

### 3. マネージャクラス

| # | テストファイル | 対象モジュール | 行数 | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|------|---------|--------------|-----------|
| 10 | `test_logmanager.py` | `logmanager.py` | 1313 | 20 | `__init__`, `getversion`, `getcondition`, `readelemslog`, `parse` (サンプルログ統合テスト) | [x] 完了 |
| 11 | `test_cpfmanager.py` | `cpfmanager.py` | 989 | 44 | `__init__`, `flatten`, `functor`, `selectfrag`, `read_header`, `parse`, `write` (往復テスト), `setupfragstr` | [x] 完了 |

### 4. スタンドアロン関数モジュール

| # | テストファイル | 対象モジュール | 行数 | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|------|---------|--------------|-----------|
| 12 | `test_readcif.py` | `readcif.py` | 930 | 36 | `getcartesiancellvec`, `getcartesianmol`, `intocell`, `getsymcoord` (22パターン) | [x] 完了 |
| 13 | `test_getifiepieda.py` | `getifiepieda.py` | 318 | 25 | `get_args` (8パターン), `setupmode` (16パターン) | [x] 完了 |

---

### 6. CLIスクリプト (argparseテスト)

| # | テストファイル | 対象モジュール | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|---------|--------------|-----------|
| 28 | `test_cli_scripts.py` | 全14 CLIスクリプト | 50 | `get_args()` / `parse_args()` のデフォルト値、必須引数、カスタム値 | [x] 完了 |

テスト対象CLIスクリプト: `addsolvfrag`, `ajf2config`, `ajfserial`, `convertcpf`, `cpf2ifielist`, `generate_difie`, `generateajf`, `getcharge`, `getifiepieda`, `log2config`, `log2cpf`, `pdb2fmo`, `pdbmodify`, `udf2fmo`

---

### 5. サブパッケージ

| # | テストファイル | 対象モジュール | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|---------|--------------|-----------|
| 14 | `test_system_model.py` | `core/system_model.py` | 10 | CellGeometry.is_rectangular, dataclass defaults | [x] 完了 |
| 15 | `test_gro_parser.py` | `gro2udf/gro_parser.py` | 7 | parse_frames, _parse_header, _parse_atom_line | [x] 完了 |
| 16 | `test_top_parser.py` | `gro2udf/top_parser.py` | 14 | get_bond_name, get_angle_name, get_torsion_name, is_improper, parse | [x] 完了 |
| 17 | `test_mdp_parser.py` | `gro2udf/mdp_parser.py` | 8 | parse_mdp, MdpParams, load_mdp | [x] 完了 |
| 18 | `test_top_model.py` | `gro2udf/top_model.py` | 6 | mass_to_element, KB constant, dataclass creation | [x] 完了 |
| 19 | `test_gro_adapter.py` | `gro2udf/gro_adapter.py` | 8 | _norm, _dot, molecule boundary, cell conversion | [x] 完了 |
| 20 | `test_gro_writer.py` | `udf2gro/gromacs/writers/gro_writer.py` | 6 | write, _cell_line, rectangular/triclinic | [x] 完了 |
| 21 | `test_top_writer.py` | `udf2gro/gromacs/writers/top_writer.py` | 8 | _strl, _strr, _f2s, full write with sections | [x] 完了 |
| 22 | `test_mdp_writer.py` | `udf2gro/gromacs/writers/mdp_writer.py` | 6 | write, NVT/NPT/vel_gen/constraints | [x] 完了 |
| 23 | `test_density.py` | `amorphous/density.py` | 10 | weight_fractions_to_counts, estimate_box_size_nm | [x] 完了 |
| 24 | `test_amorphous_models.py` | `amorphous/models.py` | 8 | ComponentSpec, BuildConfig, JSON roundtrip | [x] 完了 |
| 25 | `test_ndx_writer.py` | `amorphous/ndx_writer.py` | 5 | write_ndx, group structure, atom IDs | [x] 完了 |
| 26 | `test_mdp_protocol.py` | `amorphous/mdp_protocol.py` | 8 | generate_*_mdp, write_all_mdp, _format_mdp | [x] 完了 |
| 27 | `test_pyscf_parsers.py` | `geomopt/pyscf_optimizer.py` | 11 | _parse_xyz, _parse_pdb, _write_xyz, _infer_element, QMOptimizer init | [x] 完了 |

---

## 未テスト機能一覧

### コアモジュール (abmptools/*.py)

#### molcalc.py (31関数 未テスト / 53関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `Exportardpos` | 座標データのファイル出力 | ファイルI/O | tmpfile + 内容検証 |
| `exportdata` | 汎用データ出力 | ファイルI/O | tmpfile + 内容検証 |
| `exportplus1pos` | +1座標の出力 | ファイルI/O | tmpfile + 内容検証 |
| `exportxyz` | XYZ形式出力 | ファイルI/O | tmpfile + 内容検証 |
| `babelxyzpdb` | OpenBabel経由 XYZ→PDB変換 | obabel必須 | importorskip |
| `getpos` | 座標の取得・配置 | 座標データ依存 | 合成座標で検証 |
| `getrepeatpos` | 繰り返し構造の座標生成 | 座標データ依存 | 合成座標で検証 |
| `getmolnum` | 分子数の取得 | 座標データ依存 | 合成座標で検証 |
| `gettotalmass` | 全質量計算 | 純粋計算 | 合成データで検証 |
| `get2vec` | 2点間ベクトル計算 | 純粋計算 | 合成座標で検証 |
| `label1match` | ラベルマッチング | 純粋計算 | 合成データで検証 |
| `dihed_rotate` | 二面角回転 | 幾何計算 | 既知角度で検証 |
| `dum1_rotate` | ダミー原子回転 | 幾何計算 | 既知角度で検証 |
| `xymatch` | XY座標マッチング | 幾何計算 | 合成座標で検証 |
| `moveMolEuler` | オイラー角による分子回転 | 幾何計算 | 既知回転で検証 |
| `moveMolRotat` | 指定軸まわりの分子回転 | 幾何計算 | 既知回転で検証 |
| `calcLJPairInteraction` | LJペア相互作用エネルギー | 純粋計算 | 既知パラメータで検証 |
| `calcCoulombInteraction` | クーロン相互作用エネルギー | 純粋計算 | 既知電荷で検証 |
| `getcontactlist` | コンタクトリスト生成 | 座標データ依存 | 合成座標で検証 |
| `getcontactfrag` | フラグメントコンタクト | 座標+フラグメント依存 | fixture |
| `getrenumfrag` | フラグメント再番号付け | データ変換 | 合成データで検証 |
| `getchg` | 電荷取得 | ファイルI/O | fixture |
| `getdummyatom` | ダミー原子生成 | 純粋計算 | 合成データで検証 |
| `parse_lammps_data` | LAMMPS dataファイルパース | ファイルI/O | fixture |
| `parse_lammps_trajectory` | LAMMPSトラジェクトリパース | ファイルI/O | fixture |
| `group_atoms_by_molecule` | 分子ごとの原子グループ化 | 純粋計算 | 合成データで検証 |
| `scale_to_real_coords` | スケーリング座標変換 | 純粋計算 | 合成データで検証 |
| `get_radii_for_molecules` | 分子半径取得 | 純粋計算 | 合成データで検証 |
| `getindex` | インデックス取得 (一部テスト済み) | — | — |
| `getatomisite` | 原子サイト取得 (一部テスト済み) | — | — |
| `getmolmass` | 分子量計算 (一部テスト済み) | — | — |

#### abinit_io.py (21関数 未テスト / 31関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `gen_ajf_bodywrap` | AJFボディ生成ラッパー | 複合依存 | gen_ajf_bodyの拡張テスト |
| `saveajf` | AJFファイル保存 | ファイルI/O | tmpfile + 内容検証 |
| `readajf` | AJFファイル読み込み | ファイルI/O | fixture |
| `captfmomp2e` | FMO-MP2エネルギー取得 | ログパース | fixture |
| `captmom_single` | モノマーエネルギー取得 | ログパース | fixture |
| `captpb_single` | PBエネルギー取得 | ログパース | fixture |
| `getTE` | 全エネルギー取得 | ログパース | fixture |
| `getfmopbenergy` | FMO-PBエネルギー取得 | ログパース | fixture |
| `getmomp2ene` | MO-MP2エネルギー取得 | ログパース | fixture |
| `read_ifiepb` | PB-IFIE読み込み | ログパース | fixture |
| `getifiesumpb` | PB-IFIE合計 | ログパース | fixture |
| `del_out` | 出力ファイル削除 | ファイル操作 | tmpdir |
| `unpack_tar` | tarアーカイブ展開 | ファイル操作 | tmpdir + tar fixture |
| `getfragdict` | フラグメント辞書生成 | データ変換 | 合成データで検証 |
| `getmb_frag_seclists` | 多体フラグメントセクション | データ変換 | fixture |
| `writemb_frag_section` | 多体フラグメント書き出し | ファイルI/O | tmpfile |
| `modifyfragparam` | フラグメントパラメータ変更 | データ変換 | fixture |
| `getfraginfo` | フラグメント情報取得 (一部テスト済み) | — | — |
| `read_ifie` | IFIE読み込み (一部テスト済み) | — | — |
| `getmo_or_fmo` | MO/FMO判定 (一部テスト済み) | — | — |
| `chkdepth` | 深さチェック (一部テスト済み) | — | — |

#### pdb_io.py (6関数 未テスト / 12関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `devidepdb` | PDB分割 | 複合データ依存 | fixture PDB |
| `exportardxyzfull` | XYZ全原子出力 | ファイルI/O | tmpfile |
| `getpermol` | 残基ごとデータ分割 | データ変換 | 合成データで検証 |
| `getpermol2` | 分子ごとデータ分割 | データ変換 | 合成データで検証 |
| `getcontact_rmapfmopdb` | FMO近傍分子抽出 | 座標+フラグメント | fixture |
| `moveintocellpdb` | PBC折り返し | 幾何計算 | 合成座標+セルで検証 |

#### udf_io.py (17関数 未テスト / 23関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `Exportpos` | 座標出力 | UDFManager必須 | mock |
| `Exporttgtmolpos` | 対象分子座標出力 | UDFManager必須 | mock |
| `Exportspecificpos` | 指定座標出力 | UDFManager必須 | mock |
| `exportpdb` | PDB出力 | UDFManager必須 | mock |
| `gettotalAtm` | 全原子数取得 | UDFManager必須 | mock |
| `getmolatomnum` | 分子原子数取得 | UDFManager必須 | mock |
| `getcellsize` | セルサイズ取得 | UDFManager必須 | mock |
| `getnamelist` | 名前リスト取得 | UDFManager必須 | mock |
| `getsigmaepsilon` | LJパラメータ取得 | UDFManager必須 | mock |
| `getinteractionsitetable` | 相互作用サイト表 | UDFManager必須 | mock |
| `getatomtype` | 原子タイプ取得 | UDFManager必須 | mock |
| `moveintocell` | PBC折り返し | UDFManager必須 | mock |
| `moveintocell_mol` | 分子PBC折り返し | UDFManager必須 | mock |
| `putnvtnewfile` | NVT新規UDF作成 | UDFManager必須 | mock |
| `putnvtnewfilemb` | NVT多体UDF作成 | UDFManager必須 | mock |
| `putPositionsMol` | 分子座標書き込み (一部テスト済み) | — | — |
| `getposatom` | 原子座標取得 (一部テスト済み) | — | — |

#### udfcreate.py (36関数 未テスト / 39関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `gen_udf` | UDF全体生成 | UDFManager必須 | mock |
| `putheader` | ヘッダー書き込み | UDFManager必須 | mock |
| `putsimulationcondition` | シミュレーション条件 | UDFManager必須 | mock |
| `putstructure` | 構造情報書き込み | UDFManager必須 | mock |
| `putinitialstructure` | 初期構造書き込み | UDFManager必須 | mock |
| `putmolecularattributes` | 分子属性書き込み | UDFManager必須 | mock |
| `putsetofmolecules` | 分子セット書き込み | UDFManager必須 | mock |
| `putatom` | 原子情報書き込み | UDFManager必須 | mock |
| `putbond` | 結合情報書き込み | UDFManager必須 | mock |
| `putangle` | 角度情報書き込み | UDFManager必須 | mock |
| `puttorsion` | 二面角情報書き込み | UDFManager必須 | mock |
| `putatomparam` | 原子パラメータ書き込み | UDFManager必須 | mock |
| `putbondparam` | 結合パラメータ書き込み | UDFManager必須 | mock |
| `putangleparam` | 角度パラメータ書き込み | UDFManager必須 | mock |
| `puttorsionparam` | 二面角パラメータ書き込み | UDFManager必須 | mock |
| `putljparam` | LJパラメータ書き込み | UDFManager必須 | mock |
| `putinteractionsite` | 相互作用サイト書き込み | UDFManager必須 | mock |
| `putinteractionsitetype` | 相互作用サイトタイプ | UDFManager必須 | mock |
| `putinteractions` | 相互作用書き込み | UDFManager必須 | mock |
| `putpointcharge` | 点電荷書き込み | UDFManager必須 | mock |
| `putpos` | 座標書き込み | UDFManager必須 | mock |
| `putclusterpos` | クラスタ座標書き込み | UDFManager必須 | mock |
| `clusterfix` | クラスタ固定 | UDFManager必須 | mock |
| `getfflist` | 力場リスト取得 | 力場ファイル依存 | fixture |
| `gettorsionfflist` | 二面角力場リスト取得 | 力場ファイル依存 | fixture |
| `getffname` | 力場名取得 | 力場ファイル依存 | fixture |
| `getffid` | 力場ID取得 | データ変換 | 合成データで検証 |
| `getbondffparam` | 結合力場パラメータ | 力場ファイル依存 | fixture |
| `getangleffparam` | 角度力場パラメータ | 力場ファイル依存 | fixture |
| `gettorsionffparam` | 二面角力場パラメータ | 力場ファイル依存 | fixture |
| `getljparam` | LJパラメータ取得 | 力場ファイル依存 | fixture |
| `getbatff` | BAT力場パラメータ | 力場ファイル依存 | fixture |
| `getconnectdata` | 結合データ取得 (一部テスト済み) | — | — |
| `getbatdata` | BAT データ取得 (一部テスト済み) | — | — |
| `setudfparam` | UDFパラメータ設定 (一部テスト済み) | — | — |

#### setfmo.py (6関数 未テスト / 7関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `getfragtable` | FMOフラグメントテーブル構築 | UDFManager + PDB依存 | mock + fixture |
| `getpolyconf_rmapfmo` | ポリマー構造設定 | UDFManager依存 | mock |
| `gettermfrag` | 末端フラグメント識別 | UDFManager依存 | mock |
| `getendatom` | 末端原子識別 | UDFManager依存 | mock |
| `getcontact_rmapfmo` | FMO近傍分子抽出 | UDFManager依存 | mock |
| `make_abinput_rmap` | AJFファイル生成 | UDFManager + ファイルI/O | mock + tmpfile |

#### anlfmo.py (35関数 未テスト / 45関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `read_ifiepieda` | IFIE/PIEDA読み込み | ログパース | fixture |
| `read_pieda` | PIEDA読み込み | ログパース | fixture |
| `read_pbifiepieda` | PB-IFIE/PIEDA読み込み | ログパース | fixture |
| `read_ifpif90` | Fortranバイナリ読み込み | f90ライブラリ依存 | importorskip |
| `read_ifpimulti` | 複数構造IFIE読み込み | ログパース | fixture |
| `readlog` | ログファイル読み込み | ファイルI/O | fixture |
| `read_fraginfo` | フラグメント情報読み込み | ファイルI/O | fixture |
| `filterifiewrap` | IFIEフィルタリングラッパー | DataFrame依存 | 合成DataFrameで検証 |
| `getfiltifpiff` | フラグメント-フラグメントフィルタ | DataFrame依存 | 合成DataFrameで検証 |
| `getfiltifpifd` | 距離フィルタ | DataFrame依存 | 合成DataFrameで検証 |
| `getfiltifpifm` | 分子フィルタ | DataFrame依存 | 合成DataFrameで検証 |
| `getifiesummol` | 分子別IFIE合計 | DataFrame依存 | 合成DataFrameで検証 |
| `getsumdf` | DataFrame合計 | DataFrame依存 | 合成DataFrameで検証 |
| `gettgtdf_fd` | 距離ターゲットDF | DataFrame依存 | 合成DataFrameで検証 |
| `gettgtdf_n2ffmatrix` | N:N FFマトリックス | DataFrame依存 | 合成DataFrameで検証 |
| `gettgtdf_n2tfmatrix` | N:N TFマトリックス | DataFrame依存 | 合成DataFrameで検証 |
| `gettgtpidf_n2ffmatrix` | PIEDA FFマトリックス | DataFrame依存 | 合成DataFrameで検証 |
| `gettgtpidf_n2tfmatrix` | PIEDA TFマトリックス | DataFrame依存 | 合成DataFrameで検証 |
| `getlogchg` | ログ電荷取得 | ログパース | fixture |
| `getlogchgall` | 全電荷取得 | ログパース | fixture |
| `getlogresp` | RESP電荷取得 | ログパース | fixture |
| `getlognpa` | NPA電荷取得 | ログパース | fixture |
| `getlogmul` | Mulliken電荷取得 | ログパース | fixture |
| `getlogorpdbfrag` | ログ/PDBフラグメント取得 | ファイルI/O | fixture |
| `getmolfrags` | 分子フラグメント取得 | データ変換 | 合成データで検証 |
| `getallmolfrags` | 全分子フラグメント | データ変換 | 合成データで検証 |
| `getconnect` | 結合情報取得 | PDB依存 | fixture |
| `getbssedf` | BSSE DataFrame | ログパース | fixture |
| `getpbpiedadf` | PB-PIEDA DataFrame | ログパース | fixture |
| `writecsvwrap` | CSV出力ラッパー | ファイルI/O | tmpfile |
| `getifiedf` | IFIE DataFrame (一部テスト済み) | — | — |
| `getpiedadf` | PIEDA DataFrame (一部テスト済み) | — | — |
| `getmomenedf` | モノマーエネルギーDF (一部テスト済み) | — | — |
| `getdimenedf` | ダイマーエネルギーDF (一部テスト済み) | — | — |
| `gettgtdf_ff` | FFターゲットDF (一部テスト済み) | — | — |

#### cpfmanager.py (16関数 未テスト / 23関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `read_atominfo` | 原子情報読み込み | CPFパース | fixture |
| `read_fraginfo` | フラグメント情報読み込み | CPFパース | fixture |
| `read_dimdist` | ダイマー距離読み込み | CPFパース | fixture |
| `read_dipole` | 双極子読み込み | CPFパース | fixture |
| `read_condition` | 条件読み込み | CPFパース | fixture |
| `read_static` | 静的データ読み込み | CPFパース | fixture |
| `read_monomer` | モノマーデータ読み込み | CPFパース | fixture |
| `read_dimer` | ダイマーデータ読み込み | CPFパース | fixture |
| `set_header` | ヘッダー書き込み | 文字列生成 | 既知入力で検証 |
| `set_atom` | 原子データ書き込み | 文字列生成 | 既知入力で検証 |
| `set_mom` | モノマーデータ書き込み | 文字列生成 | 既知入力で検証 |
| `set_dimer` | ダイマーデータ書き込み | 文字列生成 | 既知入力で検証 |
| `set_static` | 静的データ書き込み | 文字列生成 | 既知入力で検証 |
| `set_condition` | 条件データ書き込み | 文字列生成 | 既知入力で検証 |
| `readfragcpf` | CPFフラグメント読み込み | CPFパース | fixture |
| `info` | デバッグ情報出力 | デバッグ用 | 低優先度 |

#### logmanager.py (12関数 未テスト / 16関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `gettotalenergy` | 全エネルギー取得 | ログパース | fixture |
| `getmominfo` | モノマー情報取得 | ログパース | fixture |
| `getdipoleinfo` | 双極子情報取得 | ログパース | fixture |
| `getlogresp` | RESP電荷取得 | ログパース | fixture |
| `getlognpa` | NPA電荷取得 | ログパース | fixture |
| `getlogmul` | Mulliken電荷取得 | ログパース | fixture |
| `getfragchgs` | フラグメント電荷取得 | ログパース | fixture |
| `readifiepieda` | IFIE/PIEDA読み込み | ログパース | fixture |
| `readpdb` | PDB読み込み | ファイルI/O | fixture |
| `getfraginfo` | フラグメント情報取得 | ログパース | fixture |
| `getnr` | レコード数取得 | ログパース | fixture |
| `extend_list` | リスト拡張 | 純粋計算 | 合成データで検証 |

#### udfrm_io.py (1関数 未テスト / 4関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `convert_udf_pdb` | UDF→PDB変換 (一部テスト済み) | 完全な結合テスト未実施 | mock UDFManager |

#### mol_io.py (1関数 未テスト / 7関数中)

| 関数名 | 概要 | 未テスト理由 | テスト方針 |
|--------|------|-------------|-----------|
| `convert_xyzs_pdb` | 複数XYZ→PDB変換 | ファイルI/O | tmpdir + fixture |

### CLIスクリプト メインロジック (全14スクリプト)

現状: argparse (`get_args()` / `parse_args()`) のみテスト済み。メイン実行ロジックは未テスト。

| スクリプト | メインロジック概要 | 依存 | テスト方針 |
|-----------|-------------------|------|-----------|
| `addsolvfrag.py` | 溶媒フラグメント追加 | setfmo, UDFManager | mock + fixture PDB |
| `ajf2config.py` | AJF→config変換 | abinit_io | fixture AJF |
| `ajfserial.py` | AJF連番生成 | ファイルI/O | tmpdir + fixture |
| `convertcpf.py` | CPFバージョン変換 | cpfmanager | fixture CPF |
| `cpf2ifielist.py` | CPF→IFIEリスト | cpfmanager | fixture CPF |
| `generate_difie.py` | DIFIE平均CPF生成 | cpfmanager, multiprocessing | fixture CPF |
| `generateajf.py` | AJF生成 | setfmo | fixture PDB |
| `getcharge.py` | 電荷取得 | anlfmo | fixture ログ |
| `getifiepieda.py` | IFIE/PIEDA取得 (`setupmode`テスト済み) | anlfmo | fixture ログ |
| `log2config.py` | ログ→config変換 | logmanager | fixture ログ |
| `log2cpf.py` | ログ→CPF変換 | logmanager, cpfmanager | fixture ログ |
| `pdb2fmo.py` | PDB→FMO変換 | setfmo | fixture PDB |
| `pdbmodify.py` | PDB編集 | pdb_io | fixture PDB |
| `udf2fmo.py` | UDF→FMO変換 | setfmo, UDFManager | mock |

### サブパッケージ 未テストモジュール

| モジュール | 概要 | 未テスト理由 | テスト方針 |
|-----------|------|-------------|-----------|
| `gro2udf/udf_writer.py` | UDFフレーム書き込み | UDFManager必須 | mock |
| `gro2udf/exporter.py` | GRO→UDFエクスポート | UDFManager必須 | mock |
| `gro2udf/top_exporter.py` | TOP→UDFエクスポート | UDFManager必須 | mock |
| `gro2udf/top_adapter.py` | TOPデータ変換 | UDFManager必須 | mock |
| `udf2gro/udf_adapter.py` | UDFデータ変換 | UDFManager必須 | mock |
| `udf2gro/exporter.py` | UDF→GROエクスポート | UDFManager必須 | mock |
| `udf2gro/gromacs/writers/itp_writer.py` | ITPファイル書き込み | NotImplemented stub | 実装後にテスト |
| `amorphous/builder.py` | アモルファス系構築 | openmm, openff, rdkit, packmol | importorskip |
| `amorphous/molecule_prep.py` | 分子前処理 | rdkit必須 | importorskip |
| `amorphous/packing.py` | Packmolパッキング | packmol必須 | importorskip |
| `amorphous/parameterizer.py` | 力場パラメータ化 | openmm, openff必須 | importorskip |
| `geomopt/mace_optimizer.py` | MACE構造最適化 | ase, mace-torch, torch必須 | importorskip |
| `geomopt/openff_openmm_minimizer.py` | OpenFF最小化 | openmm, openff必須 | importorskip |
| `geomopt/pyscf_optimizer.py` | PySCF最適化 (`optimize`) | pyscf必須 | importorskip |

### テスト優先度

#### 高優先度 (純粋計算・mock可能)

- `molcalc.py`: `calcLJPairInteraction`, `calcCoulombInteraction`, 回転系関数群, LAMMPS系関数群
- `cpfmanager.py`: `read_*` / `set_*` 系パーサー (CPF fixture で検証可能)
- `logmanager.py`: `gettotalenergy`, `getfragchgs`, `extend_list` (ログ fixture で検証可能)
- `anlfmo.py`: `filterifiewrap`, `getfiltifpi*`, `getsumdf` (合成 DataFrame で検証可能)

#### 中優先度 (fixture 作成が必要)

- `abinit_io.py`: エネルギー取得系 (`captfmomp2e`, `getTE` 等)
- `pdb_io.py`: `devidepdb`, `getpermol`, `moveintocellpdb`
- CLIスクリプト: `convertcpf`, `cpf2ifielist`, `ajf2config` (薄いラッパーで fixture で検証可能)

#### 低優先度 (外部ライブラリ必須)

- `udfcreate.py`: UDFManager 全依存 (36関数)
- `udf_io.py`: UDFManager 全依存 (17関数)
- `setfmo.py`: UDFManager 依存 (6関数)
- サブパッケージ: amorphous/builder, geomopt/mace_optimizer 等

---

## 実行方法

```bash
pytest tests/ -v              # 全テスト実行 (658 tests)
pytest tests/ -v -k molcalc   # molcalcのみ
pytest tests/ -v --tb=short   # 簡潔な出力
```
