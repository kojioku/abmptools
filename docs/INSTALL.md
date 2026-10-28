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
| UDF 変換 (`gro2udf` / `udf2gro` / `udfcharge`) | OCTA 同梱の UDFManager (OCTA の Python から使うか、§6.1 の方法で通す) | §6 |

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

UDFManager は **OCTA (GOURMET) に含まれ、PyPI には無い**。使い方は 2 通り:

- **OCTA に同梱の Python から使う。** UDFManager が最初から通っているので、
  abmptools をその Python に入れる (§1)
- **手元の Python から使う。** OCTA 同梱のプリビルドの UDFManager を通す (§6.1)。
  GOURMET を再ビルドする必要は無い

### 6.1 手元の Python から使う — プリビルドの `.so` (Linux / WSL)

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
