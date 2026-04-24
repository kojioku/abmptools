# amorphous builder — Windows インストール & 実行

メインの [amorphous_tutorial.md](amorphous_tutorial.md) は Linux/macOS を想定。
本ドキュメントは **Windows ユーザーが実行する場合の追加手順** をまとめたもの。

## 要約

| 項目 | Linux/macOS | Windows (native) | Windows (WSL2) |
|---|---|---|---|
| openff-toolkit / interchange | ✅ conda-forge | ✅ conda-forge | ✅ conda-forge |
| OpenMM | ✅ | ✅ | ✅ |
| RDKit | ✅ | ✅ | ✅ |
| Packmol | ✅ | ✅ (conda-forge win-64) | ✅ |
| **AmberTools (`sqm` for AM1-BCC)** | ✅ | ❌ **Windows native なし** | ✅ |
| **openff-nagl (ML 電荷)** | ✅ | ✅ | ✅ |
| GROMACS | ✅ | △ conda-forge あり、CPU のみ | ✅ |

**おすすめは WSL2 経路** (丸ごと Linux stack が動くので既存
[amorphous_tutorial.md](amorphous_tutorial.md) をそのまま適用可)。
どうしても native Windows で動かしたい場合は AmberTools の代替として
`openff-nagl` を使う **パッチ適用が必要**。

---

## 経路 A (推奨): WSL2 + 既存の Linux 手順

WSL2 上の Ubuntu 22.04 LTS に [amorphous_tutorial.md](amorphous_tutorial.md)
の手順をそのまま適用。変更点なし。

```powershell
# PowerShell (管理者)
wsl --install -d Ubuntu-22.04
```

あとは WSL 内で:

```bash
# Ubuntu
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba shell init -s bash

# 以降は amorphous_tutorial.md と同じ
micromamba install -n abmptoolsenv -c conda-forge -y \
    openff-toolkit openff-interchange openmm rdkit packmol ambertools
```

**メリット**:
- AmberTools がそのまま使える (AM1-BCC 電荷)
- Linux native の全スタックが動く
- ファイル I/O は Windows との共有で十分速い (`/mnt/c/...`)

**デメリット**:
- WSL2 では NVIDIA OpenCL ICD が無いため、conda-forge 版 GROMACS の
  OpenCL ビルドは GPU を使えない。小系 (数百原子) の MD は CPU 8 コアで
  10-30 分なので実用上 OK。大系は Linux native の CUDA GROMACS が無難

---

## 経路 B: Native Windows + openff-nagl

AmberTools が使えないので、**ML ベースの openff-nagl で AM1-BCC 電荷を
代替計算**する。`openff-nagl` は学習済みグラフニューラルネットで
AM1-BCC 電荷を 0.01-0.02 e の誤差で再現するライブラリ。

### B-1. 環境構築

```powershell
# PowerShell
# micromamba / mamba を未導入なら:
# 公式インストーラ: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html

# 環境作成
mamba create -n abmptoolsenv -c conda-forge -y python=3.10

# OpenFF stack + Packmol + GROMACS (AmberTools は入れない)
mamba install -n abmptoolsenv -c conda-forge -y `
    openff-toolkit openff-interchange openff-nagl openff-nagl-models `
    openmm rdkit packmol gromacs

# abmptools を editable install (本プロジェクトのルートで)
mamba run -n abmptoolsenv pip install -e .
```

> `openff-nagl` のモデル (.pt ファイル) は `openff-nagl-models` に同梱。
> 初回ロード時に展開される。

### B-2. 動作確認

```powershell
mamba activate abmptoolsenv

# OpenFF/NAGL が読める
python -c "from openff.toolkit import Molecule; from openff.nagl_models import validate_nagl_model_path; print('OK')"

# Packmol バイナリ
where packmol
# 例: C:\Users\...\mamba\envs\abmptoolsenv\Library\bin\packmol.exe

# abmptools.amorphous.cli が読める
python -c "from abmptools.amorphous.cli import main; print('OK')"
```

### B-3. 電荷計算の選び方

abmptools 本体 (main, 2026-04-24 以降) には `--charge_method` フラグが入って
いるため、**追加パッチは不要**。3 つの電荷バックエンドから選ぶ:

| `--charge_method` | 動作 | Windows native |
|---|---|---|
| `""` (default) / `am1bcc` | Interchange の AM1-BCC (AmberTools `sqm` 必須) | ❌ |
| `nagl` | `openff-nagl` の ML AM1-BCC 近似を **各分子に事前焼き込み** し、Interchange には `charge_from_molecules=molecules` で渡す (sqm を呼ばない) | ✅ |
| `gasteiger` | Gasteiger 電荷 (高速、精度は低め、saneness check 向け) | ✅ |

`--nagl_model` でモデルファイル名を差し替え可 (default:
`openff-gnn-am1bcc-0.1.0-rc.3.pt`)。

### B-4. 実行

```powershell
mamba activate abmptoolsenv

cd sample\amorphous\ketoprofen

# run_sample.sh は Linux 向け (AM1-BCC 想定)。Windows では直接 Python 実行:
python ..\..\..\build_amorphous.py `
    --smiles "OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1" `
    --name ketoprofen --n_mol 50 --density 0.8 `
    --temperature 300 --seed 42 `
    --charge_method nagl `
    --output_dir . -v
```

サンプルの `run_sample.sh` を Windows でもそのまま使いたい場合は、スクリプトを
コピーして末尾の `build_amorphous.py ...` 行に `--charge_method nagl` を足す。

#### 補足: `activate` が使えない環境の場合

企業 PC などで `mamba activate` が PowerShell の ExecutionPolicy 等で
使えない場合、`micromamba run -n <env> <command>` でその都度起動する形でも
同じ効果が得られる (env の PATH が CMD プロセスに注入され、`packmol` 等も
自動的に見つかる):

```powershell
cd sample\amorphous\ketoprofen
micromamba run -n abmptoolsenv python ..\..\..\build_amorphous.py `
    --smiles "OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1" `
    --name ketoprofen --n_mol 50 --density 0.8 `
    --charge_method nagl `
    --output_dir . -v
```

動作確認 (B-2) も同様に `micromamba run -n abmptoolsenv python -c ...` /
`micromamba run -n abmptoolsenv where packmol` で代替できる。

MD:

```powershell
cd md
bash run_all.sh      # Git for Windows の bash
bash wrap_pbc.sh
```

### B-5. nagl と AM1-BCC (sqm) の精度比較

| 指標 | AM1-BCC (sqm) | openff-nagl |
|---|---|---|
| 平均絶対誤差 vs QM | reference | ~0.015 e |
| 速度 | 分子 50 個で 約 60 秒 | 約 1-2 秒 |
| Windows native | ❌ | ✅ |

ML の学習範囲外 (荷電種・遷移金属・不自然な原子価) では誤差が大きくなる
可能性があるので、系によっては Linux/WSL2 + AmberTools で確認推奨。

---

---

## 既知の制約

### Packmol の挙動差異

Windows 版 Packmol (conda-forge) はほぼ Linux 版と同じ動作だが、出力 PDB の
改行コードが **CRLF** になる場合がある。MDAnalysis / OpenFF の PDB
パーサーは通常両対応だが、もし `Unexpected line ending` エラーが出たら:

```powershell
# dos2unix で LF 化
dos2unix sample\amorphous\<name>\build\mixture.pdb
```

### GROMACS on Windows

conda-forge 版 `gromacs` は Windows でも動くが、以下の制約あり:

- **CPU only** (CUDA / OpenCL ビルドなし)
- **MPI 並列なし** (thread-MPI のみ)
- 小系 (数百〜数千原子) なら単ノード CPU で実用時間。大系は WSL2 + CUDA 推奨

### 長パス問題

Windows のデフォルトで 260 文字超のパスでエラー。`sample\amorphous\` 以下で
深い階層を作るとぶつかる可能性あり。回避:

```powershell
# レジストリで長パス有効化 (管理者 PowerShell)
New-ItemProperty -Path "HKLM:\SYSTEM\CurrentControlSet\Control\FileSystem" `
    -Name "LongPathsEnabled" -Value 1 -PropertyType DWORD -Force
```

---

## つまずきポイント (Windows 固有)

### `FileNotFoundError: [WinError 2] 指定されたファイルが見つかりません: 'packmol'`

**原因**: PATH に `<env>\Library\bin` が通っていない。

**対処**:
```powershell
# conda/mamba activate 後に確認
where packmol
# 無ければ:
mamba activate abmptoolsenv
```

### `ImportError: DLL load failed while importing openmm`

**原因**: OpenMM の Windows binary が Visual C++ Redistributable に依存。

**対処**: [Microsoft VC++ 2019-2022 Redistributable](https://aka.ms/vs/17/release/vc_redist.x64.exe) をインストール。

### `assign_partial_charges` で `ChargeMethodUnavailableError: openff-gnn-am1bcc-...`

**原因**: `openff-nagl-models` 未インストール、またはモデル名の version 違い。

**対処**:
```powershell
mamba install -n abmptoolsenv -c conda-forge openff-nagl-models
python -c "from openff.nagl_models import list_available_nagl_models; print(list_available_nagl_models())"
# 利用可能なモデル名を確認し、パッチの "openff-gnn-am1bcc-..." 文字列を合わせる
```

### bash スクリプトが CRLF で動かない

**原因**: Git for Windows の `core.autocrlf` で `run_sample.sh` が CRLF に変換される。

**対処**:
```powershell
# 本リポジトリだけ LF 固定
cd abmptools
git config core.autocrlf input
# 既存のチェックアウトを LF に戻す
git rm --cached -r sample\amorphous
git checkout sample\amorphous
```

または `.gitattributes` を追加 (本プロジェクト側で対応推奨):
```
*.sh text eol=lf
```

---

## 次のステップ

- [amorphous_tutorial.md](amorphous_tutorial.md) — 全般的なチュートリアル
  (サンプル実行、出力構造、MD 実行、つまずき一般)
- [amorphous.md](amorphous.md) — CLI / API リファレンス
- `--charge_method` フラグはメインラインに入済み (abmptools 1.15.5+ 相当)。
  本ドキュメントで要求される追加パッチはありません
