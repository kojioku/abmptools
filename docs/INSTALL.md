# インストール

`abmptools` 本体は PyPI から入る。**手間がかかるのは本体ではなく、使う機能ごとの
外部依存**のほうで、これが OS によって入るもの・入らないものに分かれる。この
ドキュメントが導入手順の唯一の正で、Linux / WSL2 / Windows / macOS を 1 本にまとめてある。

| | |
|---|---|
| ライセンス | Apache-2.0 (自由に使える) |
| 入手 | PyPI (`pip install abmptools`) |
| Python | 3.8 以上。`amorphous` を使うなら 3.10 を勧める (OpenFF の conda パッケージが揃う) |
| 管理者権限 | 不要。venv か conda 環境に入れる |

---

## 1. 本体を入れる

```bash
pip install abmptools
```

重い依存はすべて実行時まで遅延インポートするので、**外部依存が何も無くても
`import abmptools` は成功する**。足りないものは、その機能を呼んだ時点で分かる。

入ったかの確認は **リポジトリの外** で行う。ソースツリーがカレントにあると
そちらが優先されて、確認にならない:

```bash
cd /tmp
pip show abmptools | head -2
python -c "import abmptools; print(abmptools.__file__)"
```

`Location` と `__file__` が食い違っていたら、別のコピーが勝っている。

開発するなら editable で入れる (`--user` は不要):

```bash
pip install -e .
```

---

## 2. 何を追加で入れるか — 機能で決まる

| 使う機能 | 追加で要るもの | 節 |
|---|---|---|
| CPF / IFIE / PIEDA の読み書き、ajf 生成 | **無し** (numpy / pandas のみ) | — |
| `abmptools.amorphous` (アモルファス構造構築) | OpenFF 一式 + Packmol + 電荷バックエンド | §3 |
| `abmptools.geomopt` の `pdbopt` | MACE か OpenFF のどちらか | §4 |
| `abmptools.geomopt` の `qmopt` | PySCF + geomeTRIC | §4 |
| `abmptools.trajectory` | GROMACS (`gmx`) | §5 |
| UDF 変換 (`gro2udf` / `udf2gro` / `udfcharge`) | OCTA 同梱の UDFManager (経路が 3 つある) | §6 |

`pip install 'abmptools[amorphous]'` のような extras もあるが、**extras は宣言だけ**で、
OpenMM や AmberTools のように pip では現実的に入らないものが含まれる。
実際の導入は conda-forge (micromamba / mamba / conda) 経由で行う。

---

## 3. `abmptools.amorphous` の依存

OpenFF でパラメータ化 → Packmol で初期配置 → GROMACS 用の一式を出す、という流れなので、
その 3 つが要る。

### 3.1 Linux / macOS / WSL2 — これ 1 本で終わる

```bash
micromamba create -n abmptoolsenv -c conda-forge -y python=3.10 \
    openff-toolkit openff-interchange openmm rdkit packmol ambertools
micromamba activate abmptoolsenv
pip install abmptools "setuptools<81"
```

micromamba が無ければ先に入れる:

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba shell init -s bash
```

**WSL2 でも中身は Linux と同じ**なので、これがそのまま使える。Windows から使うなら
まずこの経路を勧める (`wsl --install -d Ubuntu-22.04`)。

### 3.2 Windows native — 2 箇所だけ変える

**AmberTools が Windows に無い**ので、AM1-BCC を計算する `sqm` が使えない。
代わりに `openff-nagl` (学習済みグラフニューラルネットで AM1-BCC を 0.01–0.02 e の
誤差で再現する) を入れる。

```powershell
mamba create -n abmptoolsenv -c conda-forge -y python=3.10
mamba install -n abmptoolsenv -c conda-forge -y `
    openff-toolkit-base openff-interchange openff-nagl openff-nagl-models `
    openmm rdkit packmol gromacs
mamba run -n abmptoolsenv pip install abmptools "setuptools<81"
```

Linux との違いはこの 2 点だけ:

- **`openff-toolkit` ではなく `openff-toolkit-base`。** メタパッケージのほうは
  AmberTools に hard-depend していて、Windows では conda が solve できない。
  `-base` は RDKit バックエンドで AmberTools に依存しない
- **`ambertools` の代わりに `openff-nagl` + `openff-nagl-models`。**
  実行時に `--charge_method nagl` を付ける (§3.4)

`mamba activate` が ExecutionPolicy 等で使えない企業 PC では、`micromamba run -n <env> <command>`
で毎回起動しても同じ結果になる (env の PATH がプロセスに入るので `packmol` も見つかる)。

### 3.3 何がどこで要るのか

| パッケージ | 用途 | Linux/macOS/WSL2 | Windows native |
|---|---|---|---|
| `openff-toolkit` | 分子準備・力場適用 (SMILES / SDF) | ✅ | ✅ `-base` を使う |
| `openff-interchange` | GROMACS 形式 (gro/top) への書き出し | ✅ | ✅ |
| `openmm` | Interchange のバックエンド | ✅ | ✅ |
| `rdkit` | SDF 読み込みとサニタイズ | ✅ | ✅ |
| `packmol` | 初期配置 (外部バイナリ) | ✅ | ✅ (conda-forge win-64) |
| `ambertools` | AM1-BCC 電荷 (`sqm`) | ✅ **推奨** | ❌ 無い |
| `openff-nagl` + `-models` | ML の AM1-BCC 近似 | ✅ | ✅ **Windows はこれ** |
| `gromacs` | 生成した MD 一式の実行 | ✅ | △ CPU のみ |

`abmptools.amorphous.pubchem` (PubChem から SDF / SMILES を取る) は Python の追加依存が
無く、`https://pubchem.ncbi.nlm.nih.gov` へ出られることだけが要る。オフラインなら
先に SDF を落として `--mol` で渡す。

### 3.4 電荷バックエンドの選び方

| `--charge_method` | 中身 | Windows native |
|---|---|---|
| 省略 / `am1bcc` | Interchange の AM1-BCC。AmberTools の `sqm` が要る | ❌ |
| `nagl` | `openff-nagl` の ML 近似を各分子に焼き込み、`charge_from_molecules` で渡す (`sqm` を呼ばない) | ✅ |
| `gasteiger` | Gasteiger 電荷。速いが精度は低く、動作確認向け | ✅ |

`--nagl_model` でモデルを差し替えられる (既定 `openff-gnn-am1bcc-0.1.0-rc.3.pt`)。

### 3.5 入ったかの確認

```bash
python -c "from openff.toolkit import Molecule; print('openff OK')"
python -c "from abmptools.amorphous.cli import main; print('amorphous OK')"
packmol < /dev/null | head -3          # Windows: where packmol
gmx --version | head -1                # MD まで回すなら
```

Windows で NAGL 経路を使うなら:

```powershell
python -c "from openff.nagl_models import validate_nagl_model_path; print('nagl OK')"
```

---

## 4. `abmptools.geomopt` の依存

### 4.1 `pdbopt` — どちらか一方

| バックエンド | 入れるもの |
|---|---|
| MACE (ML ポテンシャル) | `pip install ase mace-torch torch` |
| OpenFF | `conda install -c conda-forge openmm openff-toolkit rdkit` |

OpenFF 側は §3 の環境がそのまま使える。

### 4.2 `qmopt` — PySCF

```bash
pip install pyscf geometric        # 構造最適化ドライバは geometric を推奨
pip install simple-dftd3           # D3(BJ) 分散補正 (任意、推奨)
```

- **`pyberny` は勧めない。** 0.6.3 は setuptools 82+ で `pkg_resources` の
  インポートエラーになる (2026-02 確認)。`geometric` を使う
- 分散補正がどれも入っていない場合は、警告を出して dispersion 無しで走る (エラーにはならない)

動作確認済みの組み合わせ (2026-02): pyscf 2.12.1 / geometric 1.1 / numpy 1.26.4。

---

## 5. 外部ツール

| ツール | 要る機能 | 備考 |
|---|---|---|
| **GROMACS** (2020+、2025.4 で確認) | `abmptools.trajectory`、amorphous の MD | Windows native installer もあるが、CUDA GPU を使うなら WSL2 か Linux |
| **UDFManager** | UDF 変換 (`gro2udf` / `udf2gro` / `udfcharge`) | OCTA / GOURMET 同梱で PyPI には無い。取り方は §6 |
| **VMD** | 生成された `*_pbc.xtc` などの可視化 | 任意 |

**WSL2 の GPU について**: WSL2 には NVIDIA の OpenCL ICD が無いので、conda-forge 版
GROMACS (OpenCL ビルド) は GPU を使えない。数百原子なら CPU 8 コアで 10–30 分なので
実用上は問題ないが、大きい系は Linux native の CUDA ビルドを使う。

---

## 6. UDFManager が要る場合 — どこから調達するか

COGNAC の UDF を扱う機能 (`gro2udf` / `udf2gro` / `udfcharge`) だけが `UDFManager` を
要る。**GROMACS の gro/top/xtc しか触らないなら本章は不要。**

PyPI には無く C 拡張なので GOURMET のビルドが要ると思われがちだが、**再ビルドせずに
済む経路が 3 つある**。

| 経路 | 走らせる場所 | UDFManager をどこから取るか | 向いている場面 |
|---|---|---|---|
| **A-1** | WSL / Linux | OCTA 同梱のプリビルド `.so` | 系の構築から変換まで Linux で通す (**推奨**) |
| **A-2** | WSL から Windows の python を呼ぶ | Windows の J-OCTA 同梱 `.pyd` | Linux 側に OCTA が無い |
| **B** | Windows (cmd / PowerShell) | J-OCTA の Python が最初から通している | J-OCTA で可視化しながら回す |

### 6.1 経路 A-1 — WSL でプリビルドの `.so` を使う (推奨)

**GOURMET のビルドは要らない。** OCTA 同梱の
`$OCTA/GOURMET/lib/linux_64/UdfManagerPython.so` が Ubuntu でもそのまま動く。

落ちる原因は**共有ライブラリ 2 本だけ** (Red Hat 系の soname で Ubuntu の既定と違う):

| 必要 | Ubuntu の既定 | パッケージ |
|---|---|---|
| `libjpeg.so.62` | `libjpeg.so.8` しかない | `libjpeg62` |
| `libGLU.so.1` | 未インストールのことが多い | `libglu1-mesa` |

症状は `ImportError: libjpeg.so.62: cannot open shared object file`。

```bash
# root がある場合
sudo apt install libjpeg62 libglu1-mesa

# root が無い場合 (システムに何も入れない)
apt-get download libjpeg62 libglu1-mesa
for f in *.deb; do dpkg-deb -x "$f" ~/octa-libs/sysroot; done
export LD_LIBRARY_PATH="$HOME/octa-libs/sysroot/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH"
```

```bash
export OCTA85_HOME="$HOME/OCTA85"
export UDF_DEF_PATH="$OCTA85_HOME/ENGINES/udf:$UDF_DEF_PATH"
export PYTHONPATH="$OCTA85_HOME/GOURMET/python3:$OCTA85_HOME/GOURMET/lib/linux_64:$PYTHONPATH"

python -c "from UDFManager import UDFManager; print('OK')"
```

> **`libjpeg.so.8` を `libjpeg.so.62` に symlink してはいけない。** 6.2 と 8 は ABI が
> 違う。上の手順は本物の libjpeg 6.2 を置いている。

Ubuntu 22.04 / Python 3.10.19・3.10.20・3.11.11 の 3 つで import できることを確認済み
(`.so` は `libpython3.10.so.1.0` にリンクしているが 3.11 でも動いた)。

> **マニュアルの読み違いに注意**
> `INSTALL_LINUX_v8.5_jp.pdf` の「Red Hat 系以外は Python および GOURMET を
> コンパイルする必要がある」は **GOURMET GUI 一式**の話。`.so` 1 個を使うだけの
> 用途には当たらない。

### 6.2 経路 A-2 — WSL から Windows の J-OCTA python を呼ぶ

Linux 側に OCTA が無い場合。**Windows の J-OCTA には `.pyd` が同梱**されており、
ビルドも不足ライブラリも要らない。

**WSL は環境変数を Windows プロセスへ素通ししない。** `PYTHONPATH` を普通に export
しても届かず `ModuleNotFoundError: No module named 'UDFManager'` になる。`WSLENV` に
載せて初めて渡る。

```bash
PY='/mnt/c/J-OCTA-11.1/bin/win64/Python310/python.exe'
PP="C:\J-OCTA-11.1\python;C:\J-OCTA-11.1\lib\win64"
PP="$PP;C:\OCTA8.5\ENGINES\python3;C:\OCTA8.5\GOURMET\python3;C:\OCTA8.5\ENGINES\lib\win64"
UD='C:\J-OCTA-11.1\def_udf\jp;C:\OCTA8.5\ENGINES\udf'

WSLENV=PYTHONPATH:UDF_DEF_PATH PYTHONPATH="$PP" UDF_DEF_PATH="$UD" \
  "$PY" -m abmptools.gro2udf "$(wslpath -w in.gro)" ...
```

| 項目 | 内容 |
|---|---|
| `WSLENV=PYTHONPATH:UDF_DEF_PATH` | **これが無いと届かない**。`:` 区切りで複数指定 |
| `PYTHONPATH` | **Windows 形式** (`\` 区切り、`;` 連結)。`wslpath -w` で変換 |
| `UDF_DEF_PATH` | 忘れると `No data definition[EngineType]` が大量に出る |
| 入出力パス | `wslpath -w` で Windows 形式に |

J-OCTA python の同梱: numpy 1.23.5 / pandas 1.5.3 / scipy 1.10.0 / matplotlib 3.6.3 /
rdkit 2022.09.5 / networkx 3.1。**MDAnalysis は入っていない。**

### 6.3 経路 B — Windows (J-OCTA の Python)

J-OCTA 同梱の Python (`C:\J-OCTA-11.1\bin\win64\Python310\python.exe`、3.10.9) は
**UDFManager も `UDF_DEF_PATH` も最初から通っている**。ここに abmptools を足す。

> **環境は J-OCTA 起動時に組み立てられ、そこから開いたコンソールが引き継ぐ。**
> 素の cmd には `PYTHONPATH` も `UDF_DEF_PATH` も入っていない (レジストリにも無い)。
> 素の cmd / PowerShell から使うなら先に `jocta_env.bat` を読ませる (§6.3.2 の手順 3)。
> `where python` が J-OCTA の Python310 を先頭に返せば正しい場所にいる。

入れ方は 2 通り。**既定は `--user --no-deps`** —— 有効化が要らず、J-OCTA から
コンソールを開けばそのまま使える。依存を pip に解決させたい場合だけ venv 経路に移る。

| | **`--user --no-deps` (既定)** | venv + constraints (サブルート) |
|---|---|---|
| 使うまでの手間 | **入れて終わり。有効化は要らない** | 毎回 `jocta_env.bat` の取り込み + activate の 2 手 |
| 依存の解決 | されない。足りなければ手で足す | **される** (pin したものだけ守る) |
| 同梱への影響 | 無し | 無し |
| 壊しかけたとき | そもそも解決しない | `ResolutionImpossible` で**止まる** |
| 後片付け | `pip uninstall` | venv のフォルダを消すだけ |

> **venv 経路は順序が決まっている。** `jocta_env.bat` は `PATH` を組み立て直すので、
> activate を先にすると venv が押し出され、`python` が同梱のインタプリタに戻る。
> 必ず bat → activate の順。

#### 6.3.1 入れる — `--user --no-deps` (既定)

```cmd
python -m pip install --user --no-deps "abmptools>=2.9.0"
```

**これで終わり。有効化は要らない。** J-OCTA からコンソールを開けばそのまま使える
(`--user` で入れたものは user site に入り、同梱の Python はそこを自分の
site-packages より先に読む)。

> **moldeck も使うなら、その配布 zip に同じことをする `.bat` が入っている**
> (`install_minimum.bat` が moldeck と abmptools を、`install_for_hbond.bat` が
> MDAnalysis を足す)。abmptools は PyPI 配布なのでスクリプトを同梱していない。

- **`--no-deps` を外すと壊れる。** pip が numpy 2.x を引き、user site は J-OCTA の
  site-packages より**優先される**ので、同梱の scipy 1.10.0 が
  `requires numpy<1.27.0,>=1.19.5 ... incompatible` になる。同梱の numpy 1.23.5 /
  pandas 1.5.3 のままで abmptools は動く
- **`--user` を使ってよい数少ない場面。** J-OCTA のインストール先に書き込まず、
  管理者権限も要らない (§1 の「`--user` は不要」は通常環境の話)
- **依存は解決されない。** 足りないものが出たら、その都度 `--no-deps` 付きで
  手で足すことになる。そうなったら §6.3.2 の venv 経路に移る
- **同梱の abmptools が古いことがある。** J-OCTA 11.1 には 2.6.0 が入っていた。
  確認は `importlib.metadata.version('abmptools')` と `abmptools.__file__` の両方

---

#### 6.3.2 別環境に作る場合 — venv + constraints (サブルート)

**通常は 6.3.1 で足りる。** こちらは、依存を pip に解決させたい場合や、
同梱の Python を触らずに別の環境を作りたい場合の経路。**使うたびに有効化が
要る**。

##### 先に `--user` で入れたものを消す

user site (`%APPDATA%\Roaming\Python\Python310\site-packages`) は **J-OCTA 同梱より
優先される**ので、古い `--user` の abmptools が残っていると venv を作っても
そちらが混ざる。venv 経路に移るなら必ず先に消す。

```powershell
& 'C:\J-OCTA-11.1\bin\win64\Python310\python.exe' -m pip uninstall -y abmptools
# 残骸の確認 (何も出なければきれい)
Get-ChildItem "$env:APPDATA\Python\Python310\site-packages" -ErrorAction SilentlyContinue |
    Select-Object -ExpandProperty Name
```

##### 手順

**同梱を継承する venv** を作り、同梱の版を constraints で固定してから入れる。
J-OCTA のインストール先には一切書き込まない。

```powershell
# 1. 同梱を継承する venv を作る (最初の 1 回だけ)
& 'C:\J-OCTA-11.1\bin\win64\Python310\python.exe' -m venv --system-site-packages $HOME\joctaenv

# 2. 同梱の版を固定する制約ファイル (最初の 1 回だけ)
@'
numpy==1.23.5
pandas==1.5.3
scipy==1.10.0
matplotlib==3.6.3
networkx==3.1
'@ | Set-Content -Encoding ascii $HOME\joctaenv\constraints.txt

# 3. J-OCTA の環境を PowerShell に取り込む (毎回)
cmd /c "call `"C:\J-OCTA-11.1\bin\jocta_env.bat`" > nul 2>&1 && set" |
  ForEach-Object { if ($_ -match '^([^=]+)=(.*)$') { Set-Item -Path "env:$($matches[1])" -Value $matches[2] } }

# 4. venv を有効化して入れる (毎回 activate、install は要るときだけ)
$env:PYTHONNOUSERSITE = 1
& $HOME\joctaenv\Scripts\Activate.ps1
pip install -c $HOME\joctaenv\constraints.txt abmptools
```

確認:

```powershell
python -c "from UDFManager import UDFManager; import numpy; print('UDFManager OK / numpy', numpy.__version__)"
```

`numpy 1.23.5` のままなら同梱は無傷。`Activate.ps1` が ExecutionPolicy で弾かれる
環境では `powershell -ExecutionPolicy Bypass` で起動する。cmd なら 3 は
`call "C:\J-OCTA-11.1\bin\jocta_env.bat"`、4 は `%USERPROFILE%\joctaenv\Scripts\activate.bat`。

**この構成の要点**:

| | |
|---|---|
| `--system-site-packages` | 同梱の numpy / pandas / scipy を venv から見せる。二重に入れない |
| constraints | 巻き上げが起きたとき、黙って上書きせず `ResolutionImpossible` で止める |
| `PYTHONNOUSERSITE=1` | `--system-site-packages` を付けると **user site も有効になる** (探索順は venv → user site → J-OCTA)。上の「先に `--user` で入れたものを消す」を済ませていれば不要だが、付けておくと確実 |

> **制約に書く版は実機で確認すること。** 上は J-OCTA 11.1 の実測値。版が変われば
> 変わるので、`pip list` を見てから書く。

## 7. 詰まったときの 3 点

**`ModuleNotFoundError: No module named 'pkg_resources'`**
`openff.amber_ff_ports` が `pkg_resources` を import するため、setuptools 82 以降
(2025 以降の既定) で出る。`pip install "setuptools<81"` で通る。上流が移行すれば不要になる。

**Windows の conda が openff-toolkit を solve できない**
メタパッケージが AmberTools に hard-depend している。`openff-toolkit-base` に変える (§3.2)。

**入れたはずのバージョンと違うものが動く**
リポジトリのソースツリーがカレントにあると、そちらが勝つ。`cd /tmp` してから
`python -c "import abmptools; print(abmptools.__file__)"` で実際に読まれている場所を見る。

---

## 8. 更新

```bash
pip install --upgrade abmptools
cd /tmp && python -c "import abmptools; print(abmptools.__version__)"
```

conda 環境側 (OpenFF / GROMACS など) は普段動かす必要はない。動かすときは
環境を作り直すほうが速く、確実。

---

## 9. OS 別の早見表

機能ごとの詳細は [platform_support.md](platform_support.md) にある。要点だけ:

| | Linux | macOS | Windows native | WSL2 |
|---|:---:|:---:|:---:|:---:|
| abmptools 本体 (CPF / IFIE / ajf) | ✅ | ✅ | ✅ | ✅ |
| `amorphous` (AM1-BCC) | ✅ | ✅ | ❌ AmberTools 無し | ✅ |
| `amorphous` (NAGL 電荷) | ✅ | ✅ | ✅ | ✅ |
| `trajectory` (GROMACS 後処理) | ✅ | ✅ | ✅ | ✅ |
| GROMACS の CUDA GPU | ✅ | △ | △ | ✅ |
| UDF 変換 (UDFManager) | ✅ | ✅ | ✅ | ✅ |

**迷ったら WSL2。** Linux のスタックがそのまま動き、この文書の Linux 手順を
1 行も変えずに使える。

---

## 関連

- [amorphous.md](amorphous.md) — `abmptools.amorphous` の CLI / API リファレンス
- [amorphous_tutorial.md](amorphous_tutorial.md) — 手を動かして覚えるチュートリアル
- [dependencies.md](dependencies.md) — 依存パッケージの一覧と、なぜ要るのか
- [platform_support.md](platform_support.md) — OS 別の対応表と設計上の判断
- [faq.md](faq.md) — よくあるエラー
