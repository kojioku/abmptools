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

## テスト対象外サブパッケージモジュール

| モジュール | 理由 |
|-----------|------|
| `gro2udf/udf_writer.py` | UDFManager必須 |
| `gro2udf/exporter.py` | UDFManager必須 |
| `gro2udf/top_exporter.py` | UDFManager必須 |
| `udf2gro/udf_adapter.py` | UDFManager必須 |
| `udf2gro/exporter.py` | UDFManager必須 |
| `udf2gro/gromacs/writers/itp_writer.py` | NotImplemented stub |
| `amorphous/builder.py` | openmm, openff, rdkit, packmol必須 |
| `amorphous/molecule_prep.py` | rdkit必須 |
| `amorphous/packing.py` | packmol必須 |
| `amorphous/parameterizer.py` | openmm, openff必須 |
| `geomopt/mace_optimizer.py` | ase, mace-torch, torch必須 |
| `geomopt/openff_openmm_minimizer.py` | openmm, openff必須 |

---

## 実行方法

```bash
pytest tests/ -v              # 全テスト実行 (658 tests)
pytest tests/ -v -k molcalc   # molcalcのみ
pytest tests/ -v --tb=short   # 簡潔な出力
```
