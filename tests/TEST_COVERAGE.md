# abmptools テストカバレッジマップ

## 概要

| カテゴリ | モジュール数 | テストファイル | テスト数 | ステータス |
|---------|------------|--------------|---------|-----------|
| コアクラス (継承チェーン) | 7 | 7 | 205 | 完了 |
| 複合クラス | 2 | 2 | 41 | 完了 |
| マネージャクラス | 2 | 2 | 64 | 完了 |
| スタンドアロン関数モジュール | 2 | 2 | 61 | 完了 |
| CLIスクリプト | 11 | - | - | テスト対象外 (Phase 3) |
| サブパッケージ (gro2udf等) | 4+ | - | - | 別途管理 |
| **合計** | | **13** | **385** | **全PASS** |

---

## テストファイル一覧

### 1. コアクラス (継承チェーン: molcalc → mol_io → abinit_io → pdb_io → anlfmo)

| # | テストファイル | 対象モジュール | 行数 | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|------|---------|--------------|-----------|
| 1 | `test_molcalc.py` | `molcalc.py` | 1170 | 55 | `getdist`, `getangle`, `getcrossprod`, `getCenter`, `getoriginpos`, `gettranspos`, `getmolradius`, `calccellsize`, `moveMolTrans`, `rotate_ardz`, `rotate_ardvec`, `get_element_name`, `get_atomic_radius`, `wrap_to_primary_cell`, `shift_molecule_to_primary_cell`, `getindex`, `getatomisite`, `getmolmass`, `getrenumindex` | [x] 完了 |
| 2 | `test_mol_io.py` | `mol_io.py` | 175 | 27 | `read_mol_name`, `read_xyz`, `getatoms`, `getatom`, `convert_xyz_pdb`, `Exportpospdb` | [x] 完了 |
| 3 | `test_abinit_io.py` | `abinit_io.py` | 1826 | 39 | `__init__` defaults, `chkdepth`, `get_fragsection` | [x] 完了 |
| 4 | `test_pdb_io.py` | `pdb_io.py` | 943 | 22 | `__init__` defaults, `readpdb` (HETATM/ATOM, 単一/複数残基) | [x] 完了 |
| 5 | `test_udf_io.py` | `udf_io.py` | 696 | 7 | `getposatom`, `getposmol`, `getposmolrec`, `getnameAtom`, `getAtomtypename`, `putPositionsMol` (UDFManager mock) | [x] 完了 |
| 6 | `test_udfrm_io.py` | `udfrm_io.py` | 125 | 11 | `__init__`, `getmolname`, `moveintocell_rec`, `convert_udf_pdb` (UDFManager mock) | [x] 完了 |
| 7 | `test_udfcreate.py` | `udfcreate.py` | 1130 | 10 | `__init__`, `setudfparam`, `getconnectdata`, `getbatdata` | [x] 完了 |

### 2. 複合クラス

| # | テストファイル | 対象モジュール | 行数 | テスト数 | 主要テスト対象 | ステータス |
|---|--------------|--------------|------|---------|--------------|-----------|
| 8 | `test_setfmo.py` | `setfmo.py` | 637 | 20 | `__init__` defaults, `setrfmoparam` (全キー/部分/デフォルト) | [x] 完了 |
| 9 | `test_anlfmo.py` | `anlfmo.py` | 3478 | 21 | `__init__` defaults (全属性検証) | [x] 完了 |

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

## テスト対象外モジュール (Phase 3: 報告)

以下はCLIスクリプト (`if __name__ == "__main__"`) で、テスト可能な関数を持たないため対象外。

| モジュール | 行数 | 理由 |
|-----------|------|------|
| `pdbmodify.py` | 326 | CLIスクリプト。argparse + amptオブジェクト操作 |
| `udf2fmo.py` | 102 | CLIスクリプト。UDFManager + setfmo |
| `generateajf.py` | 272 | CLIスクリプト。setfmoの薄いラッパー |
| `getcharge.py` | 75 | CLIスクリプト。anlfmoの薄いラッパー |
| `log2config.py` | 74 | CLIスクリプト。LOGManagerの薄いラッパー |
| `ajf2config.py` | 39 | CLIスクリプト。極小 |
| `log2cpf.py` | 117 | CLIスクリプト。LOGManager + CPFManager |
| `cpf2ifielist.py` | 84 | CLIスクリプト。CPFManagerの薄いラッパー |
| `ajfserial.py` | 64 | CLIスクリプト。ファイル操作中心 |
| `addsolvfrag.py` | 321 | CLIスクリプト。setfmoベース |
| `pdb2fmo.py` | 103 | CLIスクリプト。setfmoベース |
| `generate_difie.py` | 250 | CLIスクリプト。CPFManagerベース + multiprocessing |
| `convertcpf.py` | 55 | CLIスクリプト。CPFManagerの薄いラッパー |

---

## サブパッケージ (別途管理)

| パッケージ | テストファイル | ステータス |
|-----------|--------------|-----------|
| `gro2udf/` | `tests/test_gro2udf_*.py` | 別途作成予定 |
| `udf2gro/` | `tests/test_udf2gro_*.py` | 別途作成予定 |
| `geomopt/` | `tests/test_pyscf_optimizer.py` | 既存 (空) |
| `amorphous/` | `tests/test_amorphous_*.py` | 別途作成予定 |
| `core/` | `tests/test_core_*.py` | 別途作成予定 |

---

## 実行方法

```bash
pytest tests/ -v              # 全テスト実行 (385 tests)
pytest tests/ -v -k molcalc   # molcalcのみ
pytest tests/ -v --tb=short   # 簡潔な出力
```
