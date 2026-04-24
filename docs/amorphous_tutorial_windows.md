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

### B-3. 電荷計算パッチ (必須)

現行の `abmptools/amorphous/parameterizer.py` は `Interchange.from_smirnoff`
に `charge_from_molecules` を渡していないため、OpenFF のデフォルト
AM1-BCC ハンドラが起動し `sqm` を探しに行って **Windows native では失敗**
する。

以下の最小パッチを適用して `--charge_method nagl` を有効化する:

#### 差分 1/2: `abmptools/amorphous/molecule_prep.py`

`prepare_molecule` の末尾に、charges 事前計算ブロックを追加:

```python
def prepare_molecule(
    smiles: str = "",
    sdf_path: str = "",
    name: str = "",
    charge_method: str = "",     # "" (default AM1-BCC in Interchange) / "nagl" / "gasteiger"
) -> Any:
    ...  # 既存コード
    if charge_method == "nagl":
        # openff-nagl で AM1-BCC を ML で近似
        mol.assign_partial_charges(
            partial_charge_method="openff-gnn-am1bcc-0.1.0-rc.3.pt"
        )
    elif charge_method == "gasteiger":
        mol.assign_partial_charges(partial_charge_method="gasteiger")
    # charge_method == "" はそのまま通す (Interchange 側で AM1-BCC を呼ぶ)
    return mol
```

#### 差分 2/2: `abmptools/amorphous/parameterizer.py`

`create_interchange` の `Interchange.from_smirnoff` 呼び出しを条件付きに:

```python
def create_interchange(
    molecules, counts, box_size_nm, mixture_pdb,
    forcefield_name="openff_unconstrained-2.1.0.offxml",
    use_precomputed_charges=False,     # 追加
):
    ...  # 既存の topology 作成まで

    kwargs = {"force_field": ff, "topology": topology}
    if use_precomputed_charges:
        # molecules が事前に partial_charges を持っていればそれを使う
        kwargs["charge_from_molecules"] = molecules
    interchange = Interchange.from_smirnoff(**kwargs)
    return interchange
```

`builder.py` (または `cli.py`) から `charge_method` を受け取り、
`"nagl"` のときだけ両メソッドの新オプションを渡す。

> ** 重要**: このパッチを main に取り込む前に、`--charge_method` CLI フラグを
> 追加し、デフォルト (`""` = AM1-BCC) で既存挙動を保てば、Linux ユーザーには
> 影響しない。パッチ反映の際は [本プロジェクトの issue で相談](https://github.com/kojioku/abmptools/issues)。

### B-4. 実行

パッチを適用した後:

```powershell
cd sample\amorphous\ketoprofen

# run_sample.sh を bash 経由で叩く (Git for Windows の bash 等)
bash run_sample.sh

# あるいは直接 (Windows 用に書き換え)
mamba run -n abmptoolsenv python ..\..\..\build_amorphous.py `
    --smiles "OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1" `
    --name ketoprofen --n_mol 50 --density 0.8 `
    --temperature 300 --seed 42 `
    --charge_method nagl `
    --output_dir . -v
```

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

## 経路 C: Native Windows (PATCH なし、pre-charged SDF を使う)

コードを触りたくない場合の workaround。**別途 Python スクリプトで
openff-nagl 電荷を計算し、SDF に焼き込んで builder に食わせる**:

```python
# charge_with_nagl.py
from openff.toolkit import Molecule

mol = Molecule.from_smiles("OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1")
mol.generate_conformers(n_conformers=1)
mol.assign_partial_charges(
    partial_charge_method="openff-gnn-am1bcc-0.1.0-rc.3.pt"
)
mol.to_file("ketoprofen_charged.sdf", file_format="sdf")
print("Saved with partial charges (property tag 'PartialCharges')")
```

ただし **`Interchange.from_smirnoff` がこの SDF の charges を自動で引く保証は
なく、Interchange 側がデフォルトで AM1-BCC を再計算しに行く** ため、経路 B の
パッチ (特に `charge_from_molecules=molecules` を渡す部分) は結局必要。
よって C は単体では成り立たない — **経路 B がほぼ唯一の native Windows 解**
です。

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
- プロジェクトで `--charge_method` フラグが main に入れば、本ドキュメントの
  §B-3 パッチは不要になる (2026-04-24 時点ではパッチ適用が必要)
