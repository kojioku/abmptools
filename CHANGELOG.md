# Changelog

## [Unreleased]

## [2.13.2] - 2026-08-24

### Changed — tips/md-fmo: PBC 復元スクリプトを 1 本にまとめた (配布対象外)

PBC 復元が 2 本のスクリプトに分かれていた。`2_optmask-frame.sh` の 4 段
(`-pbc whole` → `cluster` → `mol/compact` → `fit`) と、`2_optmask-frame_v2.sh` の
`-pbc whole` → cpptraj `autoimage`。**4 段の方は組み直すはずの複合体を壊す**
(200 ns 系の 20 フレーム中 11 で失敗、水和殻が最大 28 Å 開く) ので、選べる形で
残しておく理由が無い。`2_optmask-frame.sh` が autoimage 手順を持つようにした。

- コードは `_v2` が走らせていたものと**バイト単位で同一**。変えたのはヘッダの
  説明だけ (「何かの代替」ではなく、それ自体が何をするかを書く形に)
- toy 系で検証: 液滴 5 つが `_v2` の実行結果と**バイト単位で一致**。
  `sample/run_demo.sh` も文書どおりの比較を再現する (cluster は 13.9 Å の穴、
  autoimage は 3.67 Å)
- `2_optmask-frame_v2.sh` はリポジトリから外した (削除ではなく
  `~/llm-project/bak/2_optmask-frame_v2_20260824.sh` へ移動、履歴にも残る)。
  新しい検査を持つのは統合先だけなので、古い方を置いておくと**検査の無い版を
  走らせる余地を残す**ことになる
- `README.md` / `2b_recover-from-center1.sh` / サンプル docs の参照を更新。
  4 段の手順を説明していた箇所は、もうそう動かないスクリプトを指すのをやめて
  手順そのものを書いた。`2b_recover-from-center1.sh` は**旧版で作った液滴専用**
  である旨をヘッダに明記

### Fixed — tips/md-fmo: autoimage の anchor マスクに溶媒が入るのを検出する (配布対象外)

以前入れた範囲検査は「総残基数を超えるマスク」しか捕まえない。toy 系
(334 残基) では十分だったが、**溶媒つきの系では役に立たない** ——
`NRES` は水も数えるので、溶質だけを指したつもりのマスクは余裕で範囲内に収まり、
検査を通ってしまう。

EGFR (`NRES=22215`) で実測: ヘッダに例として書いてあったマスク `:1-840` /
`:841-1185` は範囲検査を通り、**autoimage の anchor に水 516 残基を入れる**。
anchor は他の全フラグメントを imaging する基準なので**液滴が壊れる**。しかも
**エラーも警告も出ない** —— 範囲検査はまさにそれを防ぐためのものだった。

`%FLAG RESIDUE_LABEL` から残基名を読み、anchor / fixed マスクに
`WAT` / `HOH` / `SOL` / `TIP3` / `T3P` が含まれていたら**残基数を挙げて停止**する
ようにした。toy 系と EGFR の両方で検証済み (正しいマスクは通り、`:1-840` は
「溶媒が 516 残基入っている」で止まる)。

### Added — `gro2udf`: テンプレートが構造を持っていたら警告する

UDF は静的セクションと動的レコードに分かれ、`gro2udf` が座標とセルを書くのは
**動的レコードだけ**である。静的 `Structure` と `Initial_Structure` は
`--template` に渡した UDF の内容がそのまま残り、ここが書き換えられることは無い。

同梱テンプレートはこの部分が空なので何も主張しない。問題になるのは**実在の UDF
をテンプレートに渡したとき**で、J-OCTA や COGNAC の系を往復させるときは自然に
そうなる。出来上がる UDF は静的側に**変換前**の座標と箱、レコードに**変換後**の
それを同時に持つ。**どちらも形式としては正しいので `gmx` も COGNAC も
OCTA viewer も何も言わず**、静的側を初期構造として読む経路に乗せると、変換した
はずの構造ではなく元の構造で走る。

- 静的 `Structure` が空でないテンプレートを渡したら **`logger.warning` を出す**
  (何が起きるか + 空のテンプレートを使う手を添える)。**変換は止めない** ——
  静的側を意図的に残す使い方を壊さないため
- 判定は UDFManager 経由の `size("Structure.Position.mol[]")` で行う。ファイルを
  読んで数えると**バイナリ UDF (`.bdf`) で必ず 0 になる**ので見逃す
- 節が無い / 読めない場合は 0 として扱う。この検査が変換を失敗させる理由には
  しない (同梱の cognac101 テンプレートは静的 `Structure` 自体を持たない)

挙動そのものは変えていない。`docs/gro2udf.md` に「座標が入るのは動的レコード
だけ」の節を追加した。テスト: `tests/test_top_exporter_static_structure.py` (7 件)。


## [2.13.1] - 2026-08-24

### Fixed — `udf2gro` が `Unit_Parameter` の無い UDF を黙って誤変換していた

`udf2gro` は値を `udf.get(..., "[nm]")` のように**単位を指定して読んでおり、変換
ロジック自体は正しい**。実際の換算は UDFManager が UDF の `Unit_Parameter`
(1 sigma が何 nm か、1 epsilon が何 kJ/mol か) を基準に行う。

**問題**: `Unit_Parameter` が無いと換算は**黙って素通りする**。J-OCTA lemon が出す
全原子 UDF はこれを持たないため、Å の座標が nm、kcal/mol の epsilon が kJ/mol と
して書き出されていた。`gmx grompp` は形式が正しいので通してしまい、**箱が 10 倍
(= 密度 1/1000) の系が警告なしに走る**ため極めて気付きにくい。

同じ値 (`1.0, 4.184, 0.1`) は `gro2udf` の template が既に書いている
(上記 OCTA8.4 対応の項)。書き出す側は持っていたが、読む側が確認していなかった。

**修正**: `UdfAdapter._ensure_unit_parameter()` を追加。

- UDF が `Unit_Parameter` を宣言していればそれを使う (従来どおり)
- 宣言が無く指定も無ければ **`RuntimeError` で止める** (対処法を示すメッセージつき)
- `unit_parameter='all_atom'` (長さ Å / エネルギー kcal/mol = `(1.0, 4.184, 0.1)`) か
  `(Mass[amu], Energy[kJ/mol], Length[nm])` のタプルで明示できる

```bash
python -m abmptools.udf2gro in.udf out --unit all_atom
```

```python
Exporter().export("in.udf", "out", unit_parameter="all_atom")
```

**検証**: 修正後の出力を `gaff.dat` の実値と照合。bond C-H `b0=0.10969 nm /
kb=276646.08 kJ/mol/nm²`、angle c3-c3-hc `th0=109.80° / cth=387.4384`、
LJ c3 `σ=0.339967 nm / ε=0.4577296 kJ/mol`、1-4 対 351 組、box 2.29911 nm —
**すべて一致**。PE 20 量体 25 分子の系 5 種を GROMACS 2026.3 で完走させ、
`LJ-14` の項が出る (= 1-4 相互作用が入っている) ことも確認した。

テスト: `tests/test_udf2gro_unit_parameter.py` (7 件)。

**影響範囲**: J-OCTA の分子モデリング API (lemon 経路) が出す UDF は
`Unit_Parameter` を必ず持つ (`Length = 0.1` nm / `Energy = 4.184` kJ/mol = σ 1 Å・ε 1 kcal/mol) ため、
**この経路の挙動は変わらない**。`RuntimeError` で止まるのは手書き UDF や
`Unit_Parameter` を書かない他ツール由来のものに限られる。J-OCTA 生成 UDF での
換算は GAFF の既知値と一致することを別途確認済み (結合 c3-hc `kb` = 330.7 vs
GAFF 330.6 kcal/mol/Å²、角度 hc-c3-hc は UDF が補角 72.42° で保持しているのを
107.58° に正しく復元、電荷は `ES_Element / 18.224159264` で e 単位に一致)。

さらに、修正前の udf2gro で変換済みのポリマーバルク **105 系** (GAFF 103 種 +
DREIDING 2 種、元素 H/C/N/O/S/P/F/Cl/Br/I + Si、原子数 4,040〜8,680、結合次数
1.0/1.5/2.0/3.0 混在) を修正後の udf2gro で再変換し、`.gro` / `.top` / `.mdp` を
照合した。**105 系すべてで完全一致** (改行コードの正規化のみ)。


### Fixed — `generate_difie` が CLI として一度も動いていなかった

`getcpfobj()` と `getavestddf()` が `args` / `outcpf` を**グローバル参照**して
いた。これらは `main()` のローカル変数なので、`python -m abmptools.generate_difie`
は `NameError` で必ず落ちる。`sample/generate_difie/TrpCage/run.sh` も同様。

```
NameError: name 'args' is not defined
```

- `getcpfobj()` は `(intime, input_tmpl, zero_padding, fragments)` の
  タプルを受け取る形に変更 (`Pool.map` に渡すため 1 引数)
- `getavestddf()` は電荷ラベルを引数で受け取る
- 出力は従来どおり `<入力 stem>-DIFIE.cpf` (`gly5-xxx.cpf` → `gly5-xxx-DIFIE.cpf`)、
  対ごとに `M-` (平均) / `S-` (標準偏差) 列を持つ

CPF を 2 つ与えて実際に走らせ、`M-` / `S-` 列と `S-Total ≈ 0` を検査する
回帰テストを追加。ライブラリ経由では動いていたので、**壊れていたのは CLI 経路
だけ**で、静かに落ちるのではなく即座に例外になっていた。


## [2.13.0] - 2026-08-20

### Fixed — 分子に名前を付けると OCTA viewer の元素着色が壊れていた

`gro2udf` は OpenFF interchange の per-atom type 名を `<元素><番号>` に直して
UDF に書く (`MOL0_4` → `C4`)。**OCTA viewer は type 名の先頭 1 文字を元素記号と
解釈する**ため、この変換を怠ると `M` が Mg / Mn と読まれ、CPK 配色が壊れて
atom type ごとにランダムな色が振られる。

判定が `^MOL\d+_(\d+)$` の**リテラル固定**だった。ところが interchange は
**type 名を分子名から作る**。`Molecule.name` が未設定なら `MOL0_4`、設定済みなら
`methane_4` になる (openff-interchange 0.4.2 で確認)。つまり**利用者が分子に名前を
付けた瞬間に変換が黙って止まり**、`methane_0` がそのまま書かれて `m` が元素として
読まれていた。文字列としては妥当なので**エラーにならず、描画だけが劣化する**。

- 判定を「同じ topology 内の moleculetype 名 + `_<数字>`」に変更。この条件は
  OPLS の `opls_267` (prefix が moleculetype でない) や GAFF の `c3` を巻き込まない。
- `mol_names` 未指定時は従来の `MOL<n>_<i>` を維持するので、既存呼び出しは不変。
- テスト +13。`tests/test_interchange_adapter.py` の期待値 `("MOL0", 4)` は
  分子名が捨てられていた頃の記述だったため `("methane", 4)` に訂正
  (fixture は `name="methane"` を明示的に渡しており、保持されるのが正しい)。

### Removed — Fortran リーダを廃止し、読み込みを Python に一本化

`readifiepiedalib.so` (ctypes 経由) を使う経路と `read_ifpif90()` を削除した。
**既定の読み込み経路が Fortran から Python に変わる**が、出力は一致する
(末尾桁 1e-6 の丸めを除く)。

速度上の理由が無くなったため:

| ログ | 行数 | Python | Fortran |
|---|---:|---:|---:|
| 6lu7 164 MB | 1,502,511 | **1.09 s** | 1.95 s |
| 6m0j 43 MB | 372,816 | **0.26 s** | 0.62 s |

Python は行数に対して線形 (0.67 µs/行)。Fortran は行数によらない約 0.18 秒の
固定費 (0.7 GB の固定長配列のゼロクリア) を持つため、小さいログほど不利だった。

機能面でも Python が上回る。Fortran は 1 行につき 6 値しか読まないため MP4 の
`GRIMME-MP4` を落とし、PB-IFIE / BSSE-IFIE / monomer・dimer energy にも
対応していなかった。

- **`-nof90` / `--nof90so` は受け付けるが無視する**。既存のコマンドやスクリプト
  を壊さないために残した。
- `anlfmo` の `f90soflag` / `f90sofile` 属性は削除。代入しても副作用は無い。
- gfortran はビルドにも実行にも不要になった。`Makefile` と `abmptools/f90/` は
  参照されなくなったので削除してよい。

### Added — MP3 / MP4(CCPT) ログを Python リーダで読めるようにした

これまで Python リーダは `## MP2-IFIE` と `## HF-IFIE` の表しか見ておらず、
**MP3 / MP4 のログからは 1 行も読めていなかった** (Fortran リーダのみ対応)。

ABINIT-MP は `Method` によって IFIE 表の**列の顔ぶれ自体**を変える。
単に列が増えるのではない:

| Method | DIMER-ES より後の列 |
|---|---|
| `MP2` | HF-IFIE, MP2-IFIE, PR-TYPE1, GRIMME, JUNG, HILL (6) |
| `MP3` | HF-IFIE, MP2-IFIE, USER-MP2, MP3-IFIE, USER-MP3, PADE[2/1] (6) |
| `CCPT` | HF-IFIE, MP2-IFIE, GRIMME-MP2, MP3-IFIE, GRIMME-MP3, MP4-IFIE, **GRIMME-MP4** (7) |

3 箇所に散っていた列定義を `_IFIE_COLUMNS` に一本化し、見出しの判定・単位換算
(×627.5095)・共有結合対のマスクをすべてこの定義から導くようにした。

**`GRIMME-MP4` は Fortran リーダでは取得できない** (1 行につき 6 値しか
`read` しないため 7 列目が落ちる)。Python リーダはこの列も読む。Fortran 経路
では列を残して欠測にし、警告を出す。

### Fixed — `--ffmatrix` が MP3/MP4 で必ず失敗していた

`ValueError: Length of values (863) does not match length of index (20)`。

相手側フラグメントで絞るべきところが
`I.isin(tgt2frag) | J.isin(tgt2frag)` になっており、**注目フラグメント自身が
対象範囲に含まれると全組が通っていた**。MP2 分岐だけは相手側を `I` に寄せる
入れ替えと自己ペアの補完を行っていて正しかったので、その手順を
`_ffmatrix_partner_rows()` に切り出し、HF / MP3 / CCPT の 3 分岐も同じ経路に
乗せた。

### Fixed — IFIE が 0 行でも黙って通ることがあった

警告は「IFIE も PIEDA も BSSE も読めなかった」ときだけ出ていた。MP3 ログは
PIEDA が読めてしまうため、**IFIE が空でも無警告**で下流へ流れていた。
0 行なら手法名を添えて警告する。

### Fixed — PIEDA 節を持たないログで見出し行が表に混入していた

MP4 のログには PIEDA 節が無く、IFIE 表の終端が `## Mulliken` まで伸びる。
その見出し自体が 4 要素の行として IFIE に追加され `IndexError` になっていた。
`##` で始まる行は表の行として扱わない。

### Fixed — 共有結合フラグメント対で分散補正が捨てられていなかった

IFIE 表の `HF-IFIE < -2 Hartree` は共有結合で繋がった対の目印で、その IFIE は
物理的な相互作用を表さないため落とす規約になっている。Python 側の読み取りは
`HF-IFIE` / `MP2-IFIE` / `PR-TYPE1` の 3 列しか落としておらず、
**`GRIMME` / `JUNG` / `HILL` に生の値が残っていた** (結合対で数十 kcal/mol の
オーダー)。Fortran 側は 6 列すべてを落としており、両者の出力が食い違っていた。

Python 側を Fortran に合わせて 6 列とも落とすように修正。Python 経路が既定
なので、分散補正列を使う解析はこの値の影響を受けていた。フラグメント単位の
出力だけでなく、**合計 (`*sum*.csv`) にも伝播していた**点に注意。

回帰テストの参照 5 ファイルがこの値を焼き込んでいたため再生成した。差分が
分散補正 3 列だけに限られることを確認済み。

### Fixed — Fortran リーダの不具合 (その後この経路ごと廃止)

`--nof90` を外したときに通る経路。**いずれも従来から動作していなかった**。

- **`ctypes` の import 漏れ** — `c_char_p` / `create_string_buffer` が未 import で、
  この経路は呼ぶと必ず `NameError` になっていた。
- **固定長配列が 1 億要素** — Python 側と Fortran 側で合わせて約 25 GB を要求し、
  常用機では OOM で落ちた。実測値 (1734 フラグメント = IFIE 150 万行) に対して
  3 倍の余裕を見て 500 万要素に下げた (約 0.7 GB)。
- **配列あふれの検査が無かった** — 上限に達すると黙って切り捨てられる。
  Python 側で検知して `RuntimeError` にした。

- **`close(17)` が無く 1 プロセスで 1 ファイルしか読めなかった** — 装置番号が
  EOF に居座るため、2 回目の呼び出しは 0 件を返し、3 回目は
  `Sequential READ or WRITE not allowed after EOF marker` でプロセスごと停止
  していた。複数ログを回すモードはこの経路を繰り返し呼ぶので機能していない。
  subroutine の出口で閉じるようにした (同一ログの 4 連続読み出しで確認)。

### Performance — IFIE ログの読み取りを約 3.8 倍高速化

150 万行 (164 MB) のログで **3.85 s → 1.00 s**。行数に対して線形になった
(従来は 80 万行→150 万行でデータ 1.88 倍に対し時間 2.58 倍)。

- 節の見出しは必ず `##` で始まるので、見出し判定をフラグ 1 つで短絡させた。
  従来はデータ行ごとにリストスライスを 5 回生成していた。
- パースループを世代別 GC の停止で囲んだ。1000 万個規模の list/str を作るため
  GC が繰り返し全体を走査していた。循環参照を作らないので参照カウントだけで
  回収できる。

参考: 同じログで Fortran 経路は 1.95 s なので、**Python 単体のほうが速い**。

### Removed — 移設に伴う取り残しの整理

リポジトリを精査して見つかった残骸を落とした。機能への影響は無い。

- **`tips/udftips/` (12 本)** — 移設先に**バイト単位で同一のもの**があり、
  そちらは 19 本に育っていて tutorial から名指しで参照されている。
  こちら側を参照している doc は無かったので、二重管理を解消するために削除。
- **`docs/figures/` (4 ファイル)** — 参照していた doc が移設済みで孤立していた。
  移設先の doc が既にこの 4 つを参照していたので、そちらへ移した
  (リンク切れが 5 箇所解消)。

### Fixed — 公開物から個人環境の情報を除去

- **ハードコードされた個人パス**を除去。`.xvg` の gmx ヘッダ
  (`# Executable: /home/...`) と `.top` の `; input:` 行。
  以前 `*.py` / `*.md` / `*.sh` だけを見て「除去済み」と判断していたが、
  データファイルが対象から漏れていた。**拡張子を限定せず全テキストを走査**して
  0 件を確認。
- **内部の作業メモ名への参照 7 箇所**を落とした (`.gitignore` / CHANGELOG /
  `sample/amorphous/*/run_sample.sh` / `abmptools/udfcharge/core.py`)。
  知見の本文は残し、非公開ノートの名前と保管場所への言及だけを削った。

## [2.12.0] - 2026-08-04

### Removed — ペプチド-脂質膜 PMF を非公開リポジトリへ移設

`formulation` と同じ理由。系を組んで終わりではなく、**2 つの力場ルートを合成**し
(AMBER: ff19SB + Lipid21 + TIP3P / CHARMM36: Klauda port、CGenFF は設計上不採用)、
z-pulling → umbrella window → `gmx wham` で **PMF まで出す**。どの力場をどう
組み合わせ、window をどう配るかがコードに埋まっている。

**この移設で壊れるものは無い**。公開側に当該サブパッケージを import している
モジュールは無く、逆にこちらが公開側の何かを import してもいなかった
(完全に独立した葉)。

- extras `[membrane]` を削除。
- **`tests/regression/reference/main/udf2fmo_membrane*` は残る**。パスに
  membrane を含むが `udf2fmo` (公開ツール) のリファレンスで、無関係。
- 参照の整理: `docs/{overview,architecture,dependencies,licenses_third_party,
  faq,amorphous}.md`、`docs/ABMPTools-user-manual.md`、README。
  `docs/faq.md` の "CG: My pull stage hangs at step 0" は CG ビルダー固有の
  項目なので落とした (移設先の tutorial に同じ内容がある)。

### 注記 — 公開範囲の整理はこれで一区切り

v2.9.0 以降、開発途上の手法を順次非公開へ移した。公開を続けるのは
**確立手法のユーティリティと ABINIT-MP まわりの解析・変換**:
FMO 入力生成 / IFIE・PIEDA 解析 / CPF・log 変換 / 非晶質構築 /
GROMACS ⇄ COGNAC UDF 変換 / 構造最適化 / トラジェクトリ後処理。

## [2.11.0] - 2026-08-03

### Removed — ペプチド製剤ワークフローを非公開リポジトリへ移設

`amorphous` (公開継続) との違いは、**特定論文の再現手順そのもの**である点。
水溶液系を組み、配列からペプチドを建て、タンパク質 + 低分子 + 水 + イオンの
**混成力場**を 1 つの topology に合成し、**解析まで通す** (会合遷移 / 接触マップ /
DSSP / SASA / H 結合 / 放出 PMF)。`amorphous` は密度指定でバルクを詰めて MD 入力を
吐くところまでで、解析を持たない。

移設したもの: サブパッケージ本体、そのテスト 12 本、リファレンスとチュートリアル、
`sample/formulation`。

**この移設で壊れるものは無い**。公開側に当該サブパッケージを import している
モジュールは無かった。逆にこちらが使っていた 5 モジュールは**公開側に残る**:
`amorphous.{packing, molecule_prep, parameterizer}` / `trajectory.postprocess` /
`core.acpype`。

- extras `[formulation]` / `[formulation-analysis]` / `[formulation-openff]` を削除。
- **`.github/workflows/windows-native.yml` を削除**。このワークフローは当該
  サブパッケージの Windows ルートを回すためのもので、他に用途が無かった。
  実行履歴 5 回のうち 4 回は CI 自体の作り込み、実開発フローで走ったのは 1 回だけ
  (`paths` + pull_request トリガに対し main/develop 直 push の運用だったため)。
  これで `.github/` は空になった。
- 参照の整理: `docs/{overview,architecture,trajectory,dependencies,platform_support,
  licenses_third_party}.md`。`platform_support.md` は Windows 対応の記述が当該
  サブパッケージ前提だったため、`trajectory` と module 別対応表を残して整理した。

## [2.10.0] - 2026-08-03

### Removed — 水素結合解析を非公開リポジトリへ移設

公開を続ける線引きは従来どおり「確立した手法のユーティリティと FMO の解析・変換」。
水素結合解析は判定基準・官能基単位の役割の切り方・寿命の取り方が中身で、
**手法のノウハウそのもの**なので公開対象から外した。

移設したもの: サブパッケージ本体、そのテスト、リファレンスとチュートリアル、
および解析例のサンプル 5 件。

**この移設で壊れるものは無い**。当該サブパッケージは公開側の他モジュールから
import されておらず、逆に公開側の何かを import してもいなかった。

- `abmptools/formulation/analysis/hbond.py` は**残る**。名前は同じだが
  MDAnalysis による別実装で、GROMACS の xtc を UDF を介さず直接読む。
  `formulation` の Fig 6 はこれまでどおり動く。
- extras `[hbond]` を削除。`[rdkit]` は残る (formulation / geomopt が使う)。
- 参照が切れる箇所を整理: `docs/overview.md` の一覧と索引、
  `docs/trajectory.md` / `docs/gro2udf.md` のリンク、`docs/formulation.md` の
  対比記述、`docs/ABMPTools-user-manual.md` の該当節、README の該当節。

### Added — 非公開ワークフローの所在を README に明記

`abmptools` に無い機能を探した人が行き止まりにならないよう、別管理の
ワークフロー群がどの領域をカバーしているかを README に一覧で示した
(`## Beyond this package`)。手法の中身・パラメータには触れていない。
依存が一方向であること (あちらがこちらを import する) と、公開配布物では
ないことも明記。

## [2.9.1] - 2026-08-03

パッケージのコードに変更はない (`abmptools/` と `tests/` の差分は 0 ファイル)。
extras 名の整理とドキュメントの追随のみ。

### Added

- **extras `[rdkit]`** — rdkit を要する機能 (hbond の官能基判定 / formulation の
  小分子処理 / geomopt) 用。従来の `[fragmenter]` は**別名として残す**ので、
  既存のレシピや環境はそのまま動く。
  改名の理由は、`[fragmenter]` が名前の由来になったサブパッケージより長生きして
  しまい、インストールの案内が**このリポジトリにもう無いもの**を指していたため。
- **`tips/md/gmx/md-fmo/`** — MD-FMO 前処理の復旧スクリプト
  (`2b_recover-from-center1.sh`) と、つまずきどころ・切り分けの指針を README に
  追記。`cluster` が Fatal error で落ちる / `autoimage` の mobile は対象が空だと
  落ちる / `maskpdb` 出力の連番。検証は接触距離だけを見ず `max_gap` を併読する
  (cluster は接触を戻しつつ水和殻を置き去りにするため)。
  例示のファイル名はプレースホルダに置換済み。

### Docs

- **README にインストール手順を追加**。従来は `pip install -e .` などソースからの
  install しか書いておらず、**PyPI から入れるという最も普通の経路が抜けていた**。
  どのコピーが読まれているかの確認方法 (リポジトリの外で `pip show` と
  `__file__` を突き合わせる) も添えた。
  PyPI の wheel は pure-Python なので、ソース install 時に走る Fortran 共有
  ライブラリのビルドは行われない点も明記 (機能に影響はなく `gfortran` も不要)。
- テスト件数の記載を実測に合わせた (1271 件 / 79 ファイル)。`CITATION.cff` の
  version が 2.2.1 のまま取り残されていたのも更新。

## [2.9.0] - 2026-08-01

### 注記 — 2.8.0 の内容はこのリリースが初出
2.8.0 は取り下げ済み(下記)のため、そこに含まれていた以下の変更は 2.9.0 が
実質的な初出となる。
### Removed — 開発中手法のサブパッケージを非公開リポジトリへ移設

手法として開発途上のものを公開対象から外した。公開を続けるのは、確立した手法の
全原子 MD ユーティリティと FMO の解析・変換ツール。

**移設したもの** (テスト・ドキュメント・サンプルを含む):

| 対象 | 内容 |
|---|---|
| DPD 入力ビルダー | Cognac UDF / OCTA viewer dpm の生成、相互作用パラメータの割り当て、組成と初期構造の構築 |
| FMO 自動フラグメント分割 | 小分子・脂質・ポリマーの自動分割と CG segment 構築 |
| CG (粗視化) ビルダー | Martini 3 のペプチド / 膜モデル構築と PMF |
| FMO 結晶パイプライン | CIF から FMO 入力までの自動化 |
| MD エンジン拡張のワークフロー | 拡張サンプリングと結合自由エネルギー評価の自動化 |

**公開のまま残るもの**: root の FMO ツール群 (ajf 生成・IFIE/PIEDA 解析・
CPF/log 変換)、`amorphous` / `hbond` / `formulation` / `membrane` / `geomopt` /
`gro2udf` / `udf2gro` / `udfcharge` / `trajectory` / `core`。

**利用者への影響**:

- 移設したサブパッケージの import と、対応する CLI エントリは利用できない
- `abmptools/__init__.py` の `FragmenterConfig` / `CutSite` / `MoleculeGroup` /
  `FragmentResult` の re-export を削除
- 移設したサブパッケージ用の extras を削除。`[fragmenter]`
  (rdkit) は**名前ごと残す** — `hbond` の官能基判定、`formulation` の小分子処理、
  `geomopt` が rdkit を使い続けるため
- FMO の入力生成は root の `generateajf` / `pdb2fmo` / `addsolvfrag` に残る
  (自動分割だけが非公開側へ移った)

**公開側に残す必要があった共有部品**:

- `abmptools/core/acpype.py` — acpype (GAFF2/AM1-BCC) ラッパ。`formulation` の
  小分子パラメータ化が使う。従来は移設したサブパッケージ側にあった
- `abmptools/core/_subprocess.py` — 外部コマンド実行の薄いラッパ
  (`CommandError` / `run_command`)。従来は CG ビルダー側の内部ヘルパだった

### Added — hbond: 分子種ごとの donor/acceptor 指定 (`種名:group` + `--split-by-species`)

- **`--donor-groups` / `--acceptor-groups` が `種名:group` 形式を受け付ける**ように
  なった。多成分系で「分子 A の hydroxyl を donor、分子 B の ether を acceptor」の
  ように、官能基の種類だけでなく**どの分子種の官能基か**を限定できる:
  `--donor-groups "MOLA:hydroxyl" --acceptor-groups "MOLB:ether_O"`。
  種名は `Mol_Name`、同名種が複数ある場合 (OpenFF は全分子 `UNL`) は
  `UNL_41atoms` のように `<Mol_Name>_<原子数>atoms` で区別する (2D diagram の
  ラベルと同じ規約)。大文字小文字は区別しない。種名なしの group は従来どおり
  全分子種が対象 (**後方互換**)。
- **`--split-by-species`**: 組み合わせを列挙せずに、検出された全 donor/acceptor
  group を分子種ごとに自動分割する。IMC+PVP の ASD で
  `IMC:carboxyl -> PVP:amide_O` (薬物-ポリマー) と `IMC:carboxyl -> IMC:carboxyl_O`
  (薬物自己会合) が分離集計される。実機確認: 24 IMC + 6 PVP の GROMACS
  trajectory で rec=0 に IMC→PVP 13 本 / IMC→IMC(COOH) 4 本 / IMC→IMC(amide) 3 本。
- 種で修飾されたキーは `_pair_stats.csv` / `_pairs.csv` の `kind` / `_count.png`
  凡例にそのまま反映される。存在しない種名は警告 (`species present: ...` 付き) を
  出して空 group 扱い。`--classify-mode imc` では単一種前提の 4-species 分類のため
  種指定を明示的にエラーにする。
- **`_pairs.csv` に `donor_mol_name` / `acceptor_mol_name` 列を追加**。分割せずに
  解析して後から種別クロス集計する用途でも、分子 index を名前に引き直す手間が不要。
- API: `AnalyzerConfig.split_by_species`、`Analyzer.parse_group_spec()` /
  `Analyzer.known_species()`。tests/hbond/test_species_groups.py 18 件追加、
  **hbond 119 passed**。

### Added — hbond: Jupyter に「分子対を選ぶ」ポップアップ (auto モード)

- **Mode = auto (既定) で Run を押すと、解析する分子対を選択するポップアップが開く**。
  官能基と分子種を自動検出し、`種名:官能基` 単位の donor → acceptor 候補を一覧表示。
  各行に **(種間) / (種内)** を表示し、**全選択 / 全解除 / 種間のみ / 種内のみ**の
  ワンクリック選択を用意した。決定するとその選択のまま解析が走る。
- 選択後はパネルに**同じ解析を再現する CLI**
  (`--classify-mode auto --split-by-species --pairs "..."`) が表示されるので、
  GUI で決めた条件をそのままバッチ実行に持ち込める。「分子対を選び直す…」で再選択可。
- API/CLI 側の対応: `AnalyzerConfig.pair_filter`(評価する (donor, acceptor) の組を限定)
  と **`--pairs "IMC:carboxyl->PVP:amide_O,..."`**、候補列挙用の
  `Analyzer.candidate_pairs()`。`imc` モードでは `pair_filter` は無視される。
- **`open_panel()` が MD 入力に対応**: `open_panel("md.tpr", "md.xtc")` のように
  第 2 引数で trajectory を渡せる (CLI の `--traj` 相当)。従来は UDF/BDF 専用だった。
- Jupyter パネルの Mode 既定を **`imc` → `auto`** に変更 (CLI の既定と一致)。
- tests/hbond/test_notebook_pair_picker.py 7 件 (ipywidgets を使ったヘッドレス
  ウィジェット操作テスト。ipywidgets/UDFManager/サンプル BDF が無い環境では skip)、
  pair_filter/candidate_pairs のテスト 5 件を追加 → **hbond 131 passed**。

### Fixed — hbond: GROMACS 等の MD 入力で BDF 出力が例外になる問題

- `.tpr`/`.xtc` など MDAnalysis 経由の入力で解析すると、最後に UDF 専用の出力
  (`<prefix>.bdf` コピー / Attributes タグ / colorize) が `UDFManager` を import して
  `ModuleNotFoundError` で落ちていた。入力が COGNAC UDF/BDF でない場合はこれらを
  **スキップする** (CSV / プロット解析は従来どおり)。`--no-colorize` 等を毎回
  付ける必要がなくなった。

### Added — tips/md-fmo: PBC 復元を autoimage 化 + 溶媒つき toy サンプル (配布対象外)
- `tips/md/gmx/md-fmo/` に **cpptraj `autoimage` ベースの PBC 復元 + 接触チェッカー**
  (`check_contact.py`) を追加。GROMACS `-pbc whole/cluster/mol/compact` の 4 段
  パイプラインでは、単一 moleculetype に潰れた複合体を `whole` で再結合できず、
  `cluster` は「接触は戻すが水和殻を置き去りにする(溶質が脱水)」失敗をする。
  `autoimage` は接触・水和殻とも 1 手で決定論的に回復する。
- **toy サンプル**(`sample/make_toy.py` / `run_demo.sh` / `README.md`): 実タンパク
  座標を使わない合成系で問題を再現。18 原子版に加え、**溶媒つき ~1000 原子**
  3 フラグメント複合体(tight box で cluster 脱水を誘発)版を同梱。gmx 2021.5 +
  AmberTools cpptraj + parmed で end-to-end 確認済。
- `check_contact.py`: 溶質の接触数が均一なとき水和殻計算が空配列でクラッシュする
  不具合を修正。
- 注: `tips/` は PyPI 配布物に含まれない(`packages.find include = ["abmptools*"]`)。
### Added — formulation route B: octreotide route-B sample + end-to-end build 検証
- **`sample/formulation/octreotide_routeB_smoke/`** 新規: octreotide を GAFF whole-molecule
  でなく **tleap で残基ごと (ff14SB)** に組む config (`residue_rename={DPN:PHE,DTR:TRP}` +
  `strip_input_hydrogens` + THO `custom_residue_libs` + Cys2-Cys7 `disulfide_bonds`)。
- **builder end-to-end 実機検証**: build 完走 (2 peptide + caprate + taurocholate + TIP3P +
  JC ion、 17,799 atoms)。 Stage 5 で **Cα-H 修正が自動適用** (`repositioned 16 Cα-H atoms`
  = 8残基×2copy)、 system.gro で octreotide が **8 残基** (`PHE CYX PHE TRP LYS THR CYX THO`)、
  **res1(D-Phe)/res4(D-Trp) が D・他 L** で D キラリティ保持を確認。
- builder が `custom_residue_libs` を **絶対パス解決** (tleap は nested workdir で実行のため)。
### Added — formulation route B: D-アミノ酸を Cα-H 幾何修正で正しくモデリング
- **route B は D-アミノ酸を正しく扱える** (前版の「D で壊れる」結論を訂正)。 D→L rename 後の
  Cα 反転の原因は **tleap が Cα-H を L-テンプレの内部座標で置くこと 100%** で、 MM 項は
  キラリティ対称なので正しく組んだ D 中心は安定な極小。 **Cα-H を幾何的な四面体頂点
  (N/C/CB の反対側) に置き直せば D が保持**される。
- 新 `peptide_atomistic.fix_ca_hydrogens(gro)`: 各残基の Cα-H を幾何頂点に再配置 (どちらの
  キラリティでも正しい)。 builder が **`residue_rename` 使用時に system.gro へ自動適用**。
- **検証**: HA 修正後、 D-octreotide res1(D-Phe)/res4(D-Trp) が **aggressive 最小化 + 50 ps MD
  で D 保持** (triple product 全フレーム < 0)、 L 残基は L のまま。 修正前は同じ最小化で D→L 反転。
  pin/release テストで「熱力学的 L 選好」説は否定、 原因は H 配置と確定。 test +1。
- threoninol (THO) 残基テンプレートは前項どおり完成・検証済 (`sample/formulation/residue_libs/`)。
### Documentation — formulation route B: THO 残基テンプレート完成 (旧 D-制約記述は上で訂正)
- **threoninol (THO) 残基テンプレートを完成・検証** (`sample/formulation/residue_libs/tho.off`
  + `tho_junction.frcmod`): 主鎖+Thr 側鎖は **ff14SB THR 型+RESP 電荷** (前残基との
  peptide-bond 接合が標準 ff14SB になる)、 末端 CH2OH は ff14SB 型 (`2C/H1/OH/HO`) +
  **AM1-BCC 電荷** (ACE-threoninol を antechamber、 SDF 入力で carbonyl 誤 type 回避、
  net-0 に調整)、 欠落する 1 角度 `3C-CX-2C` を junction frcmod で補完。 **octreotide が
  net charge +2 で params 欠落なく構築、 GB 最小化+10 ps MD 安定 (NaN なし)**。
- **重大な制約を発見・明記**: **`residue_rename` で D-アミノ酸を L-テンプレートに写像しても
  D キラリティは保てない**。 初期座標は D でも L-テンプレートの力場が L を選好し、 **最小化/MD で
  Cα が L に反転** (D-octreotide res1/res4 が sander 最小化で D→L、 heavy-atom 拘束でも反転)。
  正しい D-アミノ酸には D 専用パラメータ (キラリティ補正 improper) が必要で rename では不可。
  → **D-ペプチドは GAFF whole-molecule route (キラリティ保持) + route A (`--peptide-ref`) が正しい**。
  models.py の residue_rename docstring / docs/tutorial の誤った「chirality 保持」記述を修正。
  THO テンプレートは **L-ペプチドの threoninol 末端**用。
### Added — formulation builder: tleap+pdb で非標準ペプチドを per-residue build する基盤 (route B, infra)
- whole-molecule GAFF に頼らず、 **tleap で残基ごとに組む** ための PeptideSpec 拡張:
  - `residue_rename`: 入力 PDB の非標準残基名を tleap テンプレートへ写像 (座標保持なので
    **D-アミノ酸のキラリティを維持**、 例 `{"DPN":"PHE","DTR":"TRP"}`)。
  - `strip_input_hydrogens`: H を全部落として tleap に再付加させる (H 名不一致を解消、
    重原子キラリティは不変)。
  - `custom_residue_libs`: tleap が rename で扱えない真の非標準残基 (threoninol 等) の
    AMBER 残基ライブラリ (.prep/.off/.frcmod) を loadpdb 前に source。
- `normalize_existing_pdb(residue_rename=, strip_hydrogens=)` + `render_unified_tleap_input
  (custom_residue_libs=)` (拡張子で loadAmberPrep/loadOff/loadAmberParams を判別) + builder 配線。
- **実機確認**: 1SOC octreotide を rename(DPN→PHE/DTR→TRP)+strip-H で tleap にかけると
  **残基 1-7 が native に構築**され、 C 末端 threoninol (THO) のみが唯一の要カスタム残基と
  判明。 THO の AM1-BCC 電荷は導出済 (antechamber、 SDF で結合次数を与えて carbonyl 誤 type
  を回避)。 **検証済みの THO 残基テンプレート (prepgen connect 原子 + ff14SB 接合 frcmod +
  MD 安定性確認) は focused follow-up** として残す (rush parametrize は誤った FF を生む恐れ)。
  test +3。
### Added — formulation analyze: whole-molecule peptide を per-residue 解析 (--peptide-ref, route A)
- GAFF whole-molecule route (D-octreotide 等、 1 分子=単一 residue `OCT`) では Fig 4
  (残基別接触) と Fig 7 (DSSP) が「1 残基/peptide」に潰れて無意味だったのを、
  **`--peptide-ref <単一 peptide PDB>` で残基レベルに復元**できるように。 DSSP/接触は
  *幾何*ベースなので **MD 再実行不要** — 解析用に残基ラベルを貼り直すだけ。
- 新 `analysis/residue_ref.py:build_perres_reference`: 参照 PDB (原子順が trajectory の
  1 copy と一致) の残基境界を各 copy に当て、 **全系 per-residue の解析トポロジー
  (`perres_reference.gro`)** を生成。 これを MDAnalysis 解析 (aggregate/contacts/hbond) と
  gmx dssp の `-s` に使う。 SASA は元 tpr + 元 selector のまま (vdW 半径のため)。
- `--peptide-ref-resmap "DPN:PHE,DTR:TRP,THO:THR"`: 非標準残基名→標準名の写像で
  gmx dssp が主鎖を protein 認識。 `parse_resmap` / `load_ref_labels` も公開。
- **実機検証**: D-octreotide が **48 残基** (6×8) 認識、 DSSP が per-residue SS
  (sheet 1.8% / turn 30% / coil 68%、 従来は mean_residues=6 で全 coil)、 Fig 4 も 8 残基別
  (res1 DPhe / res4 DTrp / res5 Lys が enhancer 接触上位)。 L 体 (tleap) と対称に比較可能に。
  test +3 (`parse_resmap` / `load_ref_labels` / gro 行フォーマット)。 docs/tutorial 更新。

---

## [2.8.0] - 2026-07-26 — 取り下げ (withdrawn)

**このリリースは配布されていない。** PyPI から削除済みで、タグと GitHub
Release も削除した。同一バージョン番号は PyPI の仕様上再利用できないため、
次のリリースは 2.9.0 となる。

理由はサブシステムの分離 (2.9.0 の Removed 節)。2.8.0 に含まれていた
それ以外の変更は 2.9.0 に再収録した。

---

## [2.7.0] - 2026-07-13

### Changed — hbond: `--classify-mode` の既定を `imc` → `auto` に

- 汎用 H-bond 解析ツールとして安全な既定に変更。`auto` は系に存在する官能基を自動検出して
  generic donor × acceptor pair 統計を回すため、**COOH を持たない系(PVA / peptide /
  アルコール等)を無指定で流しても意味のある結果**になる(旧既定 `imc` は COOH 中心の
  4-species 分類なので、非 COOH 系だと検出 0 で戸惑う)。
- **IMC の 4-species 分類(dual/chain/single/free)には `--classify-mode imc` を明示**する
  必要。IMC サンプル(`run_cli.sh`)/ `docs/hbond.md` / `tutorial_hbond_imc.md` /
  サンプル README を更新。
- `AnalyzerConfig.classify_mode` 既定も `"auto"`。imc regression テストは
  `classify_mode="imc"` を明示するよう更新。

## [2.6.0] - 2026-07-13

### Added — hbond: per-pair 寿命分布プロット + サンプル拡充 + doc 監査

- multi-record trajectory で `<prefix>_lifetime.png`(occupancy + 連続平均寿命の 2 パネル
  ヒストグラム)を自動出力。全体の C(t)(`_autocorr.png`)を per-pair 分布で補完。
- docs/hbond.md に「寿命 / 自己相関の読み方」節を追加。`abmptools.trajectory` の新規 doc
  (`trajectory.md`)、`formulation.md` の索引化、`overview.md` 読書順の更新、
  `tutorial_hbond_imc.md` の陳腐パス修正・現行機能反映。
- imc サンプルに RDKit 構造式 diagram(`output/imc_hbond_diagram.{png,svg}`)と
  101-record トラジェクトリの時系列結果(`output_trajectory/`、count vs record /
  lifetime / autocorr / distance)を同梱。per-frame 着色データは公開 Release から DL。

## [2.5.0] - 2026-07-13

### Changed — hbond: `--out-prefix` 省略時の既定を入力由来 `<input>_hbond` に

- CLI で `-o` / `--out-prefix` を省略したときの既定を、固定文字列 `hbond_result` から
  **「入力パスの拡張子を `_hbond` に置換したもの」**に変更。例: `result.bdf` →
  `result_hbond`(`result_hbond.bdf` / `result_hbond_summary.csv` / ... が入力の隣に出力)。
  `-o` を明示した場合は従来どおり指定 prefix をそのまま使う(挙動不変)。

### Changed — hbond: `_colored.bdf`(Mol_Name 分子色分け)を opt-in 化(既定 `--colorize-mode` を `action` に)

- `--colorize-mode` の既定を `molname`(v1.25 legacy)→ **`action`** に変更。既定出力は
  `.bdf`(素)/ `_attribute_rec{N}.bdf`(J-OCTA 属性)/ `_action.bdf`(官能基 action)の 3 種で、
  **`_colored.bdf` は既定では出力しない**。J-OCTA は別方法(action / 属性)で分子色分けできるため。
- `_colored.bdf`(imc の Mol_Name → `IMC_DUAL/CHAIN/SINGLE/FREE` リネーム + Draw_Attributes)が
  欲しいときは **`--colorize-mode molname`(colored のみ)または `both`(action + colored)を明示**。
- `AnalyzerConfig.colorize_mode` 既定も `"action"` に。`colorize_udf`(直接呼び出し)は不変。
  smoke: 既定で colored 無し、`molname` で colored 有りを確認。

### Fixed — hbond: per-frame action の色替えを復活(attr 読み取り方式・GOURMET フリーズ回避)

- `--attributes-per-record` の per-frame 実装を **「小さい BDF + `{record→roles}` 埋込テーブル
  `.act`(currentRecord 参照)」**から **「per-record Attribute を焼き込んだ BDF + 小さい attr
  読み取り `.act`」**へ戻した。前者は全 record × 全 tag 原子のテーブルで `.act` が巨大化し(IMC
  125 分子 × 875 原子 × 101 record → 約 1.6 MB)、GOURMET が autorun 毎に parse して **開いた瞬間に
  フリーズ**していた。一旦「既定の静的 action が安全」と回避したが、それでは **時系列で色(結合の
  組み換え)が変わらない**。
- 復活した方式: `colorize_udf_action_per_record` が `write_hbond_attributes_per_record` で
  `_action.bdf` に **record 毎の `hbond` Attribute** を焼き込み(topology 保持のため
  `Set_of_Molecules` を record ごと subtree copy)、`.act` は **毎フレーム `get(Attributes[].Value)`
  で読む数 KB のスクリプト**(`_render_attr_reading_body`)。`.act` が軽いので GOURMET は
  フリーズせず、record スライダで **色が変わる**。代償は `_action.bdf` が大きい(IMC ~131MB /
  APZ ~303MB / ASD ~84MB)。
- API: `colorize_udf_action_per_record` / `write_show_python_script_per_record` は
  `taggable_atoms` 引数を取る(埋込テーブル方式の signature から復帰)。currentRecord 埋込方式
  (`_render_currentrecord_body`)は撤去。test 更新(埋込テーブル assert → attr 読み取り assert)。
- 配布データ(private `abmptools-sample` release `data-hbond-perframe-20260712`)の 3 系統
  (IMC/APZ/ASD)の `_action.bdf` を per-frame attr 読み取り版に再生成(`.act` 数 KB、BDF に
  record 毎属性)。rec0/50/100 で 414 原子が role 変化することを実機確認。

### Changed — hbond: J-OCTA 属性着色を 1 record の `_attribute_rec{N}.bdf` に分離

- J-OCTA は `Attributes[]` を**起動時に 1 回だけ読み、record スライダでは読み直さない**(実機確認済)。
  そのため多 record の `<prefix>.bdf` に静的な per-atom 属性タグを付けると「スライダで色が変わるはず」と
  **誤解を生む**。→ **特定フレーム(既定=最終解析フレーム `N`)を 1 record だけ切り出した
  `<prefix>_attribute_rec{N}.bdf`** に hbond `Attributes[]` を付与するよう変更。ファイル名に record 番号を
  入れることで「これは 1 フレームの静的スナップショット」と一目瞭然。
- 出力 4 種を明確化: `.bdf`(素の多 record=通常の時系列)/ `_attribute_rec{N}.bdf`(1 record 属性)/
  `_action.bdf`(多 record・官能基だけ per-frame 着色)/ `_colored.bdf`(多 record・**OCTA 用の分子ごと
  色分け**=Mol_Name リネーム、属性とは別機能)。`.bdf` は着色タグ無しの素コピー。
- 新規: `colorizer.extract_single_record(src, dst, record)`(`UDFManager.eraseRecord` で 1 record 化、
  topology・座標保持)。CLI `--attribute-record N`(既定=最終フレーム)/ `--no-attribute-single`(従来の
  多 record 付与に戻す)。config `do_attribute_single`(既定 True)/`attribute_record`。test +2。
  実機: imc で `_attribute_rec100.bdf` 1 record・1.6 MB、topology/hbond 属性保持を確認。

### Fixed — gro2udf: 圧力に Pres. DC を加算(barostat 目標・J-OCTA と一致)

- `--energy` 経路の `Statistics_Data.Pressure` が **GROMACS の生 `Pressure`(virial のみ、
  long-range 分散補正なし)**を書いていたのを修正。GROMACS の `Pressure` は分散補正を含まず、
  **barostat 制御の真の圧力 = `Pressure` + `Pres. DC`**。`_XVG_TO_UDF_STATS` に `Pres. DC` を
  追加して `Statistics_Data.Pressure` へ fold(既存の合算機構=Proper+Improper Dih.→Torsion と同じ)。
  これで **J-OCTA の gro→bdf 変換と一致**する。従来コメントの「DC は Pressure に含まれるので足さない」
  という判断は誤りだった。
- 実機検証(DRO10-PVP10 NPT、`ref_p = 1 bar`、100k steps): avg `Pressure` 単独 = **44 MPa**(目標から
  大きく外れる)に対し avg(`Pressure`+`Pres. DC`)= **0.22 MPa ≈ 0.1 MPa** 目標で、J-OCTA-export BDF の
  `Statistics_Data.Pressure`(Total_Average 0.215 MPa)とも一致。test +2。

  - `fragment_protein(mol) -> CutSet`: protein を **CA–C 切断 (SP3) → heavy-atom
    連結成分 = fragment** で自動分割。ペプチド C–N とジスルフィド S–S (< 2.2 Å 距離判定)
    は残すので carbonyl は次 fragment に、C 末端 COOH・側鎖・SS 対は自動でまとまる。
  - 電荷: 側鎖 (ARG/LYS +1、ASP/GLU −1) + **H/O カウント末端** (N 末端 `nH(N)−1`、
    C 末端 `−(carboxyl O で H=0 の数)`)。ABINIT-MP felec 規約を再現。
  - `load_protein_pdb(pdb)`: `MolFromPDBFile(removeHs=False)` で atom 順保存
    (mol idx+1 = PDB serial、AJF ifatom と一致)。
- **回帰検証**: `fragment_protein → cutset_to_segment_data` の出力が ABINIT-MP
  `AutoFrag=ON` 参照 (`&FRAGMENT`) と **完全一致**。
  - gly5 (5-Gly, 38 atom, 5 frag) / TrpCage (20 残基, 304 atom, charge +1, ASP/LYS/ARG +
    NH3+/COO- 末端) の atom/charge/connect/connect_num/seg_info が全一致。
- 公開 API に `fragment_protein` / `load_protein_pdb` を追加。

### Changed — hbond: per-frame 着色を currentRecord() 方式に最終化(Set_of_Molecules は global 維持)

- **`--attributes-per-record` の実装を全面変更**(以下の一連の Fixed を統合・最終化)。COGNAC は
  `Set_of_Molecules` を `\begin{global_def}`(record 非依存)で宣言し 1 回だけ格納する。ここへ
  per-record Attribute を書く従来方式は、ブロックを per-record 化させ **肥大化 + topology 欠落**を
  招き J-OCTA/COGNAC モデラを壊した(全ブロック丸コピーで回避したが今度は BDF が肥大)。
- **最終方式**: `Set_of_Molecules` は **global のまま一切触らない**。per-frame は **GOURMET action**
  が `{record -> roles}` テーブルを埋め込み **`currentRecord()`**(cognac_draw.act と同じ)で現フレーム
  を引いて overlay 描画。`_action.bdf` は入力素コピー(topology 完全・肥大化なし)。`.bdf`/`_colored.bdf`
  の Attribute/Mol_Name タグは **global(`jump(-1)`)へ 1 回だけ**書く=静的(J-OCTA は静的フィルタのみ、
  per-frame は原理的に不可)。BDF サイズが正常化(imc .bdf 47.8→25 MB)。`colorize_udf_action_per_record`
  /`write_show_python_script_per_record` は埋込テーブル+`currentRecord()`、`write_hbond_attributes`
  (+generic)/`colorize_udf` は global-only 書込に変更。test 更新(action が BDF 不変+テーブル埋込を検証)。

### Fixed — hbond: per-record Attributes 書込で Set_of_Molecules(名前/型/bond)が消える

- `--attributes-per-record` で record 可変 `Set_of_Molecules` に per-record Attribute を書くと、
  **その record の Set_of_Molecules から Attribute 以外(`Mol_Name`/`Atom_Name`/`Atom_Type_Name`/
  `bond[]` 等)が丸ごと欠落する**バグを修正。原因は record 非依存の `Set_of_Molecules` ブロックが
  per-record put で per-record コピーへ昇格する際、明示的に再 put しないフィールドが引き継がれ
  ないため。**GOURMET autorun は座標から描くので無症状**だったが、**J-OCTA ポスト描画 /
  COGNAC モデラは atom name/type・分子名・結合を record から直接読むため認識不能**になっていた。
- `write_hbond_attributes_per_record` が **`Set_of_Molecules` ブロック全体を record 不変 topology
  として一度キャッシュ(`get`)し、触れた各 record で丸ごと再 put(`put`)してから Attribute を
  上乗せ**するよう修正。これで各 record の Set_of_Molecules は Attribute 追加分を除き定義と
  一致する(座標は `Structure.Position` 側で不変)。個別フィールド列挙より堅牢で高速
  (~0.17 s/record)。`colorize_udf_action_per_record` も本関数経由で修正。既に per-record な UDF
  (gro2udf 出力)では同値の再書込で無害。gro2udf 参照 UDF で回帰テスト +1(name/type + bond[]
  保持を検証)、hbond 8 passed。

### Added — hbond: per-frame (アニメーション) 着色 `--attributes-per-record`

- 時系列解析の着色を **snapshot ごとに変化**させる経路を追加。従来の `<prefix>.bdf`
  の `Attributes[]` タグは **最終フレームの role を全 record に焼き付ける**(= 再生しても
  色が固定)だったが、`--attributes-per-record`(config `attributes_per_record=True`)を
  付けると **各 record 自身の H-bond role** を `Set_of_Molecules.molecule[].atom[].
  Attributes[]`(この容器は COGNAC/gro2udf スキーマで **record 可変**)に書き込む。
  OCTA viewer / J-OCTA が `hbond` Attribute で色分け・フィルタする経路がそのまま
  per-snapshot 化し、H-bond ネットワークの生成消滅が軌道再生で追える。
- あるフレームで idle な atom は `inactive_value`(既定 `'none'`)に **リセット**して
  書く(前フレームの色が滲まない)。タグ対象は全 record で active になった atom の和集合で
  安定 → Attribute フィルタが常に同じ atom 集合を対象にできる。imc/generic 両モード対応。
- 実装: `colorizer.write_hbond_attributes_per_record` +
  `Analyzer._build_per_record_atom_value`。既定 OFF(従来の静的着色を維持)、multi-record
  時のみ有効。実機: acetaminophen 非晶質 (64 mol × 8 record) で 306 atom 全てが
  record 間で値変化することを確認。test +4。
- **GOURMET action もアニメーション化(viewer 非依存)**: `--attributes-per-record` を
  `--colorize-mode action`/`both` と併用すると、`<prefix>_show.act` が **現在 record の
  `hbond` Attribute を描画時に読んで sphere を塗る** autorun action になる
  (`colorize_udf_action_per_record`)。autorun は record 切替で再実行され `get`/`size` は
  現在 record 参照(座標が動くのと同じ機構)なので **GOURMET 上で色が snapshot ごとに変わる**
  = J-OCTA の Attribute 再読込挙動に依存しない確実な per-frame 経路。従来の静的 action
  (role を最終フレーム固定で埋込)は原子は動くが色不変だった。`_show.py`(J-OCTA Python
  パネル用)も同じ読取ロジック。実機(acetaminophen generic 5 record)で `_action.bdf` に
  per-record Attributes 同梱 + `.act` が Attributes 読取形になることを確認。test +3(計 90 passed)。

### Added — formulation route B: octreotide route-B sample + end-to-end build 検証

- **`sample/formulation/octreotide_routeB_smoke/`** 新規: octreotide を GAFF whole-molecule
  でなく **tleap で残基ごと (ff14SB)** に組む config (`residue_rename={DPN:PHE,DTR:TRP}` +
  `strip_input_hydrogens` + THO `custom_residue_libs` + Cys2-Cys7 `disulfide_bonds`)。
- **builder end-to-end 実機検証**: build 完走 (2 peptide + caprate + taurocholate + TIP3P +
  JC ion、 17,799 atoms)。 Stage 5 で **Cα-H 修正が自動適用** (`repositioned 16 Cα-H atoms`
  = 8残基×2copy)、 system.gro で octreotide が **8 残基** (`PHE CYX PHE TRP LYS THR CYX THO`)、
  **res1(D-Phe)/res4(D-Trp) が D・他 L** で D キラリティ保持を確認。
- builder が `custom_residue_libs` を **絶対パス解決** (tleap は nested workdir で実行のため)。

### Added — formulation route B: D-アミノ酸を Cα-H 幾何修正で正しくモデリング

- **route B は D-アミノ酸を正しく扱える** (前版の「D で壊れる」結論を訂正)。 D→L rename 後の
  Cα 反転の原因は **tleap が Cα-H を L-テンプレの内部座標で置くこと 100%** で、 MM 項は
  キラリティ対称なので正しく組んだ D 中心は安定な極小。 **Cα-H を幾何的な四面体頂点
  (N/C/CB の反対側) に置き直せば D が保持**される。
- 新 `peptide_atomistic.fix_ca_hydrogens(gro)`: 各残基の Cα-H を幾何頂点に再配置 (どちらの
  キラリティでも正しい)。 builder が **`residue_rename` 使用時に system.gro へ自動適用**。
- **検証**: HA 修正後、 D-octreotide res1(D-Phe)/res4(D-Trp) が **aggressive 最小化 + 50 ps MD
  で D 保持** (triple product 全フレーム < 0)、 L 残基は L のまま。 修正前は同じ最小化で D→L 反転。
  pin/release テストで「熱力学的 L 選好」説は否定、 原因は H 配置と確定。 test +1。
- threoninol (THO) 残基テンプレートは前項どおり完成・検証済 (`sample/formulation/residue_libs/`)。

### Documentation — formulation route B: THO 残基テンプレート完成 (旧 D-制約記述は上で訂正)

- **threoninol (THO) 残基テンプレートを完成・検証** (`sample/formulation/residue_libs/tho.off`
  + `tho_junction.frcmod`): 主鎖+Thr 側鎖は **ff14SB THR 型+RESP 電荷** (前残基との
  peptide-bond 接合が標準 ff14SB になる)、 末端 CH2OH は ff14SB 型 (`2C/H1/OH/HO`) +
  **AM1-BCC 電荷** (ACE-threoninol を antechamber、 SDF 入力で carbonyl 誤 type 回避、
  net-0 に調整)、 欠落する 1 角度 `3C-CX-2C` を junction frcmod で補完。 **octreotide が
  net charge +2 で params 欠落なく構築、 GB 最小化+10 ps MD 安定 (NaN なし)**。
- **重大な制約を発見・明記**: **`residue_rename` で D-アミノ酸を L-テンプレートに写像しても
  D キラリティは保てない**。 初期座標は D でも L-テンプレートの力場が L を選好し、 **最小化/MD で
  Cα が L に反転** (D-octreotide res1/res4 が sander 最小化で D→L、 heavy-atom 拘束でも反転)。
  正しい D-アミノ酸には D 専用パラメータ (キラリティ補正 improper) が必要で rename では不可。
  → **D-ペプチドは GAFF whole-molecule route (キラリティ保持) + route A (`--peptide-ref`) が正しい**。
  models.py の residue_rename docstring / docs/tutorial の誤った「chirality 保持」記述を修正。
  THO テンプレートは **L-ペプチドの threoninol 末端**用。

### Added — formulation builder: tleap+pdb で非標準ペプチドを per-residue build する基盤 (route B, infra)

- whole-molecule GAFF に頼らず、 **tleap で残基ごとに組む** ための PeptideSpec 拡張:
  - `residue_rename`: 入力 PDB の非標準残基名を tleap テンプレートへ写像 (座標保持なので
    **D-アミノ酸のキラリティを維持**、 例 `{"DPN":"PHE","DTR":"TRP"}`)。
  - `strip_input_hydrogens`: H を全部落として tleap に再付加させる (H 名不一致を解消、
    重原子キラリティは不変)。
  - `custom_residue_libs`: tleap が rename で扱えない真の非標準残基 (threoninol 等) の
    AMBER 残基ライブラリ (.prep/.off/.frcmod) を loadpdb 前に source。
- `normalize_existing_pdb(residue_rename=, strip_hydrogens=)` + `render_unified_tleap_input
  (custom_residue_libs=)` (拡張子で loadAmberPrep/loadOff/loadAmberParams を判別) + builder 配線。
- **実機確認**: 1SOC octreotide を rename(DPN→PHE/DTR→TRP)+strip-H で tleap にかけると
  **残基 1-7 が native に構築**され、 C 末端 threoninol (THO) のみが唯一の要カスタム残基と
  判明。 THO の AM1-BCC 電荷は導出済 (antechamber、 SDF で結合次数を与えて carbonyl 誤 type
  を回避)。 **検証済みの THO 残基テンプレート (prepgen connect 原子 + ff14SB 接合 frcmod +
  MD 安定性確認) は focused follow-up** として残す (rush parametrize は誤った FF を生む恐れ)。
  test +3。

### Added — formulation analyze: whole-molecule peptide を per-residue 解析 (--peptide-ref, route A)

- GAFF whole-molecule route (D-octreotide 等、 1 分子=単一 residue `OCT`) では Fig 4
  (残基別接触) と Fig 7 (DSSP) が「1 残基/peptide」に潰れて無意味だったのを、
  **`--peptide-ref <単一 peptide PDB>` で残基レベルに復元**できるように。 DSSP/接触は
  *幾何*ベースなので **MD 再実行不要** — 解析用に残基ラベルを貼り直すだけ。
- 新 `analysis/residue_ref.py:build_perres_reference`: 参照 PDB (原子順が trajectory の
  1 copy と一致) の残基境界を各 copy に当て、 **全系 per-residue の解析トポロジー
  (`perres_reference.gro`)** を生成。 これを MDAnalysis 解析 (aggregate/contacts/hbond) と
  gmx dssp の `-s` に使う。 SASA は元 tpr + 元 selector のまま (vdW 半径のため)。
- `--peptide-ref-resmap "DPN:PHE,DTR:TRP,THO:THR"`: 非標準残基名→標準名の写像で
  gmx dssp が主鎖を protein 認識。 `parse_resmap` / `load_ref_labels` も公開。
- **実機検証**: D-octreotide が **48 残基** (6×8) 認識、 DSSP が per-residue SS
  (sheet 1.8% / turn 30% / coil 68%、 従来は mean_residues=6 で全 coil)、 Fig 4 も 8 残基別
  (res1 DPhe / res4 DTrp / res5 Lys が enhancer 接触上位)。 L 体 (tleap) と対称に比較可能に。
  test +3 (`parse_resmap` / `load_ref_labels` / gro 行フォーマット)。 docs/tutorial 更新。

### Added — formulation analyze: SASA/DSSP を「簡単に使える」ように (--tpr / --auto-pbc / 親切エラー / tutorial)

- **`--tpr`**: SASA/DSSP 用の tpr を明示指定可能に。 従来は `--top` の `.gro`→`.tpr`
  sibling を暗黙解決していたため `system.gro` (sibling tpr 無し) を渡すと SASA/DSSP が
  黙って skip していた。 `--tpr prod.tpr` で解消 (未指定時は従来 fallback)。
- **`--auto-pbc`**: 解析前に `gmx trjconv -pbc mol` を自動実行し whole-molecule 軌道
  (`analysis/prod_pbcmol.xtc`、 `--target-n-frames` まで間引き) を生成、 全 Fig をそれで
  解析。 NPT raw `prod.xtc` の PBC 割れで SASA/DSSP が境界で誤る問題の手動前処理を不要に。
  新 `analysis/preprocess.py:make_whole_molecule_traj`。
- **gmx version 親切エラー**: `run_gmx_dssp` が古い gmx (`do_dssp` のみ、 GROMACS < 2023) を
  検出したら *「Fig 7 DSSP requires GROMACS >= 2023」* と明示 RuntimeError (従来は無関係な
  subprocess エラー)。
- **`docs/tutorial_formulation_smoke.md` 新規** (formulation.md の See-also が参照する
  リンク切れを解消): analyze の end-to-end 手順 (入力・1 コマンド・出力・GAFF whole-molecule
  route の `--peptide-selector`・troubleshooting 表・octreotide L worked example)。
- formulation.md quickstart に `--tpr`/`--auto-pbc`/`--gmx` を追記、 SASA/DSSP 注記を刷新。
  test +2 (`make_whole_molecule_traj` / 古い gmx の親切エラー)。 実機で raw prod.xtc +
  system.gro + `--tpr` + `--auto-pbc` 一発完走を確認。

### Documentation — formulation: analyze 各図の計算方法

- `docs/formulation.md` に **「解析手法の詳細 (how each figure is computed)」**節を
  追加。 共通前処理 (Universe / stride / atom-index 分割 / PBC・heavy-atom 0.5 nm) と、
  Fig 1b/1c/2 (heavy-atom 接触グラフ→connected components→co-cluster 行列)、
  Fig 4 (per-residue の接触相手原子数、 cap 除外)、 Fig 5 (gmx sasa DCL 法)、
  Fig 6 (MDAnalysis HBA、 3.5 Å/150°、 Pep-Pep は分子間のみ)、 Fig 7 (gmx dssp
  num.xvg → H/E/turn/coil) の各アルゴリズムを実装準拠で記述。

### Fixed — formulation analyze: DSSP (Fig 7) 配線を GROMACS 2023+ 形式に修正

- `run_analysis(run_dssp=True)` が **無効な kwarg** (`run_gmx_dssp(out_dir=, ndx_path=)`)
  で呼んでいて DSSP が常に例外 → skip されていたのを修正。 正しい `out_dat=` /
  `num_xvg=` / `selection=` を渡すように。
- `secondary_structure.parse_dssp_xpm` は旧 `do_dssp` (GROMACS <= 2022) の `.xpm`
  前提だったが、 **GROMACS 2023+ の `gmx dssp` は `.dat` (per-frame SS 文字列) +
  `-num .xvg` (SS 種別本数の時系列) を出力**する別形式。 num.xvg を読む
  `parse_dssp_num_xvg` を新設し (H=α+3₁₀+π / E=strand+bridge / turn=turn+bend /
  coil=loop+PPII、 鎖区切り `=` は残基数から除外)、 `run_gmx_dssp` を `.dat`+`-num`
  出力に更新 (`out_xpm` は後方互換 alias として残置)。
- `plots.plot_dssp` (SS 割合の stackplot) を追加し `plot_workflow_outputs` に配線。
  SASA step は `-output` group を surface selection に合わせるよう修正
  (whole-molecule route の非 protein selection でも通る)。
- **octreotide L/D 10 mM (peptide×6, 100 ns) で Fig 1-7 を初 end-to-end 検証**
  (SASA/DSSP は本 run が初検証)。 L (all-L, ff14SB) vs D (実薬 1SOC, GAFF): D は
  hexamer +14%・SASA −29% でよりコンパクトに凝集する一方 H-bond は L ≫ D。
  なお D の GAFF whole-molecule route は 1 分子 = 単一残基 `OCT` のため **Fig 4
  (残基別) と Fig 7 (DSSP) は L のみ有意**、 凝集 (Fig 1/2) と SASA (Fig 5) が
  route 非依存で比較可能。 test +3 (`parse_dssp_num_xvg` / empty / `plot_dssp`)。

### Changed — amorphous: GROMACS residue 名を `--name` から付与(UNL 廃止)

- amorphous ビルドが各成分の GROMACS residue 名を **`--name`(分子名)から**付けるように
  なった(従来は OpenFF/Interchange 既定の `UNL` で全成分同名 → 混合系で種を区別できず)。
  residue 名は **英数字・大文字・5字**に整形(`IMC`/`PVP`、`acetaminophen`→`ACETA`、空→`MOL`)。
  これで `abmptools.hbond` の混合系解析が `mol_name` で分子種を直接判別でき、原子数フォールバック
  不要。実機: IMC+PVP ASD で `IMC`/`PVP` resname + `_diagram_IMC`/`_diagram_PVP`。
- **適用は GROMACS export 時**(`parameterizer._restore_residue_names`、`from_smirnoff` 後の
  Interchange topology を各成分にマッチさせて命名)。**prep 時の PDB に非-UNL 名を書くと
  `Topology.from_pdb` の residue マッチングが壊れる**(`UnassignedChemistryInPDBError`)ため。
  `molecule_prep._residue_name_for` / `apply_residue_name`。test +3。

### Fixed — hbond diagram: 同名分子種のファイル名衝突

- 混合系で複数の分子種が **同じ `mol_name`** を持つ場合(OpenFF は全分子を `UNL` と
  命名するため頻出)、diagram のファイル名が `<prefix>_diagram_<molname>` で衝突し、
  片方が上書きされて 1 枚しか残らなかった。**原子数で曖昧性を解消**
  (`<prefix>_diagram_UNL_41atoms` / `_UNL_90atoms` 等、なお衝突するなら連番)するよう修正。
  test +1、hbond suite 96 passed。

### Added — hbond: GROMACS / MD 入力対応 (MDAnalysis reader)

- `abmptools.hbond` が COGNAC UDF/BDF 以外に **GROMACS などの MD 出力**
  (`.tpr` + `.xtc`/`.trr`、`.gro`、`.pdb`、`.psf`、AMBER `.parm7` 等)を直接解析できる
  ようになった: `python -m abmptools.hbond md.tpr --traj md.xtc ...`。新 `--traj` オプション。
- 新 module `abmptools/hbond/mda_reader.py` の `MDATrajectory` が `BDFTrajectory` と同じ
  interface(`molecules`/`n_records`/`get_frame`/`get_cell`)を MDAnalysis で提供。分子は
  結合情報から自動分割(fragments)、元素は topology 取得 or atom 名から推定(`CA`→C 等も正しく)。
  入力は拡張子で自動判別(`.udf`/`.bdf`=UDFManager、それ以外=MDAnalysis)。
- MDAnalysis は **optional 依存**(未導入なら明確な ImportError)。acetaminophen の GROMACS
  `.tpr`+`.xtc` で UDF 版と **byte 一致の結果**を確認(record 0 の 27/23/35、τ_HB 212.0347)。
  test `tests/hbond/test_mda_reader.py`(+4)、hbond suite 95 passed。`MDATrajectory` 公開。

### Fixed — formulation: octreotide_l / hexarelin_l のペプチド正味電荷 (+1 → +2)

- `octreotide_l` sample config (smoke / cluster_smoke / aggregation_10mM /
  aggregation_100ns) の `cap_n` を `ACE` → `""` に修正。ACE cap が N 末端アミンの
  +1 を潰して正味 +1 になっていたが、実薬 octreotide (および `octreotide_d` の
  `net_charge:2`) は **+2** (N 末端遊離アミン + Lys)。遊離 N 末端に戻して +2 に一致。
  disulfide RESID は builder の `cap_offset` 自動補正 (builder.py) で追従するため
  spec 変更は不要。
- `hexarelin_l` sample config の `cap_c` を `""` → `NME` に修正。遊離 C 末端酸 (-1)
  になっていたが、実薬 hexarelin の C 末端は **アミド (中性)**。NME cap で中性化し +2 に。
- 影響: tleap 中和イオン (`addions ... 0`) の本数が正しくなる (誤電荷だと Na+/Cl-
  カウントがずれ、系の ionic 環境が不正になる)。

### Fixed — formulation GAFF/tleap route が背景塩 (`salt_concentration_M`) を添加するように

- 従来 GAFF/tleap route (`topology.py`) は `salt_concentration_M` を無視して
  **neutralize のみ** (該当ブロックが `pass`) だった。**水分子数ベースの 2-pass 法**で
  背景塩を添加するよう実装 (gmx `genion` と同じ規約):
  1. **probe pass**: solvate + neutralize してから `savepdb` (パラメータは書かない)
  2. probe PDB の水分子数を数え、`n = round(salt_M / 55.5 * n_water)` で NaCl ペア数を算定
  3. **final pass**: neutralize 後に `addIonsRand sys Na+ n Cl- n` で背景塩をランダム配置
- 新 API: `render_unified_tleap_input(..., extra_addions_count=, probe_pdb_out=)`、
  `count_solvent_waters(pdb)` (WAT/HOH/SOL の O を計数、欠損時は 0 + warning)、
  `compute_salt_ion_pairs(molarity, n_water)`。OpenFF route (gmx `genion -conc`) は従来通り。
- formulation unit tests **59 passed** (abmptoolsenv、rdkit 込み)。topology に two-pass /
  probe / addIonsRand / 水計数 / 塩ペア算定の 7 test を追加。電荷修正の形式電荷も
  in-process で両系 +2 と確認。

### Added — `abmptools.formulation` 解析ワークフローに H-bond / SASA / co-cluster network (Fig 5/6/2)

- `analysis.run_analysis` を Hossain 2023 の Fig 5-7 まで拡張:
  - **Fig 6 H-bond** (`run_hbond=True`、 既定 on): **MDAnalysis HydrogenBondAnalysis** で
    xtc を直接読み、 peptide ↔ {enhancer / bile salt / peptide / water} の H-bond を
    **種別に**時系列集計 (`analysis.hbond.compute_hbonds`)。 contact_map と同じ species
    selection を再利用。 `.gro` は結合/電荷を持たないので **bond を幾何 guess + name
    ベースの donor/acceptor 明示** で電荷・tpr 不要に (GROMACS 2026 tpr は MDAnalysis
    非対応)。 注: `abmptools.hbond` は UDF 専用のため formulation では使えず MDAnalysis 経路。
    **peptide-peptide は分子間のみ**集計 (donor/acceptor が別 peptide copy、 分子内=二次
    構造 H-bond は除外。 insulin 実機で 22.9 → 1.53 H-bond/peptide)
  - **Fig 5 SASA** (`run_sasa=True`、 gmx 依存で既定 off): `gmx sasa` wrapper を配線
  - **Fig 2 co-cluster network**: aggregate の `transitions.npy` (peptide 間共クラスタ行列)
    を networkx で network 図に (`plots.plot_cocluster_network`、 networkx 無しは heatmap)
- 新規 plot: `plot_hbond_summary` (Fig 6 bar) / `plot_cocluster_network` (Fig 2) /
  `plot_sasa` (Fig 5)。 `run_analysis` の plot 一式に統合
- 実機確認: insulin 10mM (6 peptide) で Pep-Enhancer 1.52 / Pep-BileSalt 0.10 H-bond
  per peptide、 co-cluster network で P0-P1 / P2-P4 / P4-P5 の強い共クラスタを可視化
- test 4 件追加 (plot_hbond_summary / cocluster_network / sasa + compute_hbonds import)

## [2.4.1] - 2026-06-28

### Changed — hbond diagram: donor+acceptor atoms drawn magenta

- 2D H-bond-site diagram で、**donor かつ acceptor の原子**(例: OH 酸素、both 役割)を
  **マゼンタ**で描くようにした(従来は acceptor 色で上書きされ donor 役割が見えなかった)。
  donor-only=赤 / acceptor-only=シアン / both=マゼンタ で、colorize の Donor/Acceptor/Both
  規約と一致。凡例にも「magenta = donor+acceptor」を重複時に追記。
- `draw_hbond_diagram` に `both_note`(default `"donor+acceptor"`)を追加(後方互換)。
  test +1、hbond suite 91 passed。acetaminophen/PVA の phenol/alcohol OH が both で正しく表示。

## [2.4.0] - 2026-06-27

### Added — hbond ether-oxygen acceptor detection (`detect_ethers`)

- `abmptools.hbond.functional_groups.detect_ethers` を追加: **C-O-C エーテル酸素**
  (酸素が炭素 2 個と結合・H 無し)を acceptor として検出する。element + 結合グラフ
  ベースで、ester の単結合 O も拾う。OpenFF SMIRNOFF のように専用 ether atom type が
  無い FF でも動作する(従来 element fallback はエーテル O を `carbonyl_O` と tag して
  いたが、結合次数=1 の真の C=O と区別して degree=2 を ether と判定)。
- `--classify-mode auto` がエーテル O を検出時に自動で `ether_O` を acceptor に追加。
  `_build_acceptor_sites_by_type` の `ether_O` 経路を実 detector に置換(従来は no-op
  stub)。2D diagram もエーテル O を acceptor として色付けする。
- 実機: APZ (aripiprazole) で butoxy エーテル O を検出し、record 0 で
  `amide_donor→ether_O` を 6 本検出(従来は 0)。`EtherGroup` / `detect_ethers` を
  公開 API に追加。test `tests/hbond/test_functional_groups.py` (+3)、hbond suite 90 passed。

## [2.3.0] - 2026-06-26

### Added — hbond per-molecule 2D structure diagram (検出 H-bond サイトの可視化)

- 解析実行時に **分子種ごとの 2D 構造式**を自動出力し、検出した donor (赤) /
  acceptor (シアン) 原子を色分け表示する(`<prefix>_diagram[_<molname>].{png,svg}`)。
  COGNAC UDF の topology(`Atom_Name` の元素 + 結合 + record0 座標)から RDKit で
  分子を構築し、UDF に結合次数が無くても幾何から芳香環 / C=O を知覚する
  (`rdDetermineBonds.DetermineBondOrders`)。分子種は element 署名で重複排除。
- RDKit はオプション依存: 未導入 / 知覚失敗(電荷不一致等)なら警告で skip。
  `--no-diagram` で抑制、`--diagram-charge` で正味電荷を指定(default 0 = 中性)。
- 新 module `abmptools/hbond/diagram.py`(`element_from_name` /
  `build_mol_from_topology` / `draw_hbond_diagram`)+ `tests/hbond/test_diagram.py`
  (4 件)。hbond suite 87 passed。APZ 実機で lactam の amide N-H → C=O を自動色分け。

### Fixed — hbond `_count.png` が generic/auto mode でも COOH 凡例を出していた

- count プロットが classify_mode で分岐せず、常に imc の "dual/chain/single/free
  COOH-COOH mols" 凡例(generic では 0 線)を描いていた。generic/auto では
  **donor_type→acceptor_type ペア種別ごとの H-bond 本数**を描くよう修正。

## [2.2.1] - 2026-06-25

### Added — hbond `--classify-mode auto` (官能基の自動認識)

- `abmptools.hbond` に **`--classify-mode auto`** を追加(imc / generic と並列の
  3 つ目の選択肢)。系に存在する官能基(carboxyl / amide / 2 級アミド N-H /
  アミン N-H / hydroxyl)を `detect_carboxyls` / `detect_amides` /
  `detect_amine_donors` / `detect_hydroxyls` で自動検出し、対応する
  donor / acceptor group を列挙して generic pair-stats 経路で実行する。
  COOH 無し系(APZ 等)でも `--donor-groups` / `--acceptor-groups` の手動指定が
  不要になる(明示指定すればそちらが優先)。
- helper `Analyzer._auto_groups()`(検出官能基 → donor/acceptor group リスト)。
  auto は内部で `classify_mode='generic'` に解決。**`ether_O` は自動検出対象外**
  (信頼できる ether 検出器が無いため。必要なら `--acceptor-groups ether_O` を明示)。
- 実機確認: APZ を `--classify-mode auto`(group 指定なし)→ `amide_donor→amide_O`
  を自動選択、手動指定時と同一結果。test `tests/hbond/test_analyzer.py`
  (carboxyl-only / amide-only / mixed / none / run→generic 解決、5 件)。

### Added — hbond `--record-stride` (trajectory 間引き)

- `abmptools.hbond` に **`--record-stride N`**(`AnalyzerConfig.record_stride`、default 1)
  を追加。`run()` の record ループを N record ごとに間引く。大規模 multi-record UDF
  (例: 1001 record)を `--record-stride 10 --record-end 1000` で ~100 frames に削減でき、
  原子ごと UDFManager 読みがボトルネックになる長尺 trajectory の解析を高速化する。
  lifetime/τ_HB を使う場合は `--dt` も同倍率でスケール(間引くと sampled frame が離れ、
  連続性指標 `continuous` は意味が薄くなる)。hbond 専用オプション。
- test `tests/hbond/test_analyzer.py`(default / stride / start-offset / guard、4 件)。
  `docs/hbond.md` の CLI オプション一覧、`sample/amorphous/acetaminophen_amorphous`
  README に記載。commit `e31c57e`。

### Added — `abmptools.udfcharge` 形式電荷の復元 (`restore_formal_charge`)

- 中和 (Σq≈0) された 1 分子 UDF の電荷を、 指定**形式電荷 (整数)** になるよう
  逆変換して別 UDF に出力する。 MD 用に過剰電荷を分散して中和した系から、 元の
  per-atom 電荷 (Σ=形式電荷) を復元する用途 (FMO 解析等)。 **中和ルール `mode`** を
  選択 (中和方法に一致させる):
  - **`mode="proportional"`** (既定): 過剰分を `|q|` 比例で分散。
    forward `B_i = A_i − S·|A_i|/Σ|A|`、 reverse `A_i = B_i/(1∓λ)`。 λ は
    `S·λ² + (P−N)·λ + (P+N−S) = 0` (P=Σ_{B>0}B, N=Σ_{B<0}B) の |λ|<1 の根
    (`SI/A列再現方法.md` の方法)。 `|S| ≥ Σ|q|` (符号反転) は `ValueError`
  - **`mode="uniform"`**: 過剰分を全原子に均等分散。 forward `B_i = A_i − S/N`、
    reverse `A_i = B_i + S/N`。 二次方程式・符号問題なしで**常に厳密・一意**
- `restore_formal_charge(udf, formal_charge, out=None, mol_index=0, mode=...)` → `RestoreResult`
- CLI を **サブコマンド化**: `transfer` (従来の template→bulk 転写) + `restore`
  (新規)。 旧フラット呼び出し (`--template ...`) は後方互換で `transfer` 扱い
  - `python -m abmptools.udfcharge restore --udf mol.udf --formal-charge 12 [--mode uniform] --out out.udf`
- 更新は `electrostatic_Site` のみ (座標・bond/angle/torsion 無改変)
- テスト `tests/udfcharge/` を 11 → 21 件に (proportional/uniform round-trip 復元 /
  負電荷 / S=0 no-op / 符号反転 ValueError / mode 不一致検出 / 座標・結合保持 / CLI)。
  実データ (`reverse-charge.xlsx` の 864 atom、 形式電荷 12) で proportional 復元
  誤差 5.8e-10 を確認。 サンプル `sample/udfcharge/restore_example.py`
  (methylammonium +1)、 docs 追記

## [2.2.0] - 2026-06-24

### Added — `abmptools.udfcharge` (new sub-package)

- 単分子 UDF (電荷あり) の per-atom partial charge を抽出し、 バルク系 UDF の
  **同名分子すべて**へ転写する OCTA/COGNAC UDF 電荷割り当てモジュール
  (例: 量子化学 / FMO で求めた電荷を MD バルク系へ反映)
- `read_molecule_charges(udf, mol_name=/mol_index=)` → `MoleculeChargeTemplate`
  (mol_name / n_atoms / charges[e] / atom_type_names / net_charge)。
  `assign_charges_to_bulk(bulk, template, out)` → `AssignResult`
- 電荷規約は gro2udf / udf2gro と共通: `electrostatic_Site[].Type_Name="POINT_CHARGE"`、
  `.ES_Element = charge[e] × 18.224159264`、 `.atom[0] = atom index`
  (`UDFManager.put` の numpy silent-zero を避けるため `float()` cast)
- 割り当ては **atom index 対応**。 atom 数 (必須) と `Atom_Type_Name` 列
  (`verify_atom_types`) を検証してから書き込み、 不一致は strict で例外 / 非 strict で skip。
  出力は別ファイル (入力 bulk 無改変)
- 更新するのは `electrostatic_Site` のみ。 `Structure.Position` (座標) / `Unit_Cell` /
  `bond`・`angle`・`torsion` (結合トポロジー) は **無改変** (座標・トポロジー保持を回帰テストで固定)
- CLI `python -m abmptools.udfcharge --template mol.udf --bulk bulk.udf --out out.udf`
  (`--mol-name` / `--mol-index` / `--no-verify-types` / `--non-strict`)
- 単体テスト `tests/udfcharge/test_udfcharge.py` 11 件 (default_template.udf から
  methanol 系をプログラム生成、 座標保持 + トポロジー保持テスト含む)。 サンプル
  `sample/udfcharge/` は **topology + 電荷 + 実座標 + bond/angle/dihedral + Unit_Cell**
  を持つ OCTA viewer で分子の形ごと開ける UDF (methanol 単分子 + 2×2×2 格子 bulk の
  end-to-end)、 docs `docs/udfcharge.md` + `docs/tutorial_udfcharge.md`

### Added — `abmptools.formulation` Phase 2: Windows native OpenFF route

- `FormulationBuildConfig.force_field_route` で力場経路を 2 つから選択:
  - `"amber"` (default) — tleap + acpype + parmed。**Linux / macOS のみ** (AmberTools)
  - `"openff"` — **全 OS (Windows native)**。PDBFixer + OpenFF `Topology.from_pdb`
    + ff14SB SMIRNOFF で tleap / acpype を一切使わず peptide formulation を build
- OpenFF route の実装 (`topology_openff.py` / `peptide_atomistic_openff.py` /
  `builder._build_openff`):
  - **PDBFixer** で water 除去 + 欠損補完 + explicit H 付加 (OpenFF 必須)
  - **`Topology.from_pdb`** (≠ `Molecule.from_polymer_pdb`) で multi-chain protein を
    1 分子認識 + **disulfide 自動検出** (insulin 2G4M の S-S を手動宣言なしで処理)
  - protein = ff14SB library charges、small molecule = precomputed (gasteiger/am1bcc)、
    ion = ff tip3p charges を `protein_flags` で振り分け
  - **O(N²) 回避**: full mixture の `Interchange.from_smirnoff` は nonbonded exception
    生成が O(N²) で爆発する (insulin×6 で CPU 2h+) ため、**単一コピーを parametrize
    (~4 s) → `[molecules]` の count を実値に書き換え**て複製。`.gro` は packmol mixture
    から `gmx editconf`、溶媒和は `gmx solvate` + `gmx genion` (Joung-Cheatham)
  - sequence 入力は PeptideBuilder + Biopython で 3D extended chain 生成 (natural L-AA)
- amber route の multi-chain / disulfide 対応 (`topology._patch_mixture_pdb_with_ter_records`
  + `builder._disulfide_tleap_bonds`): packmol が drop する TER record を chain ID 遷移 +
  resnum reset 検出で再挿入、disulfide Cys を CYS→CYX rename + global resid map で
  `bond sys.<i>.SG sys.<j>.SG` 注入
- `ndx.py`: >99999 atom 系 (insulin 156k) で gro の atom-serial 列が 100000 で wrap し
  nvt grompp が "Invalid atom number 0" になる bug を 1-based 連番に修正 + 回帰テスト
- 新規 extras `[formulation-openff]` (openff-toolkit / openff-interchange /
  openff-amber-ff-ports / openmm / pdbfixer / PeptideBuilder / biopython)
- 新規 sample: `octreotide_{l,d}_aggregation_10mM/`、`hexarelin_l_aggregation_10mM/`、
  `insulin_aggregation_10mM/` (論文 Hossain 2023 準拠 10 mM = peptide × 6)
- 実機検証: octreotide L/D 10mM (81k/90k atoms、~10h GPU)、insulin 10mM
  (156k atoms、15.2h / 143 ns/day on RTX 4070 Ti) を 100 ns 完走。insulin は
  tetramer→hexamer の irreversible aggregation、TYR14(A)/PHE22(B-Phe1) が凝集 hotspot
- docs: `platform_support.md` (OS 別対応表 + OpenFF route 実装点)、`formulation.md` /
  `dependencies.md` に OpenFF route 節

### Added — `abmptools.formulation` Phase 2-D: Windows native CI

- `.github/workflows/windows-native.yml` (abmptools 初の CI workflow)。`windows-latest`
  runner で OpenFF route と `abmptools.trajectory` の **Windows native 動作を実機検証**:
  - `windows-pure-python` — pip のみで trajectory + formulation models/ndx/mdp unit test
    (pathlib / subprocess shell=False の OS 依存箇所を検証)
  - `windows-openff-smoke` — micromamba (conda-forge) で OpenFF stack →
    `tests/integration/test_openff_windows_smoke.py` (sequence → PeptideBuilder →
    PDBFixer → `Topology.from_pdb` → ff14SB Interchange の parametrization smoke、@slow)
  - trigger は `workflow_dispatch` (手動) + formulation / trajectory 変更の PR。
    packmol / gmx は外部ツールなので CI では呼ばず AmberTools 非依存の Python 層だけ検証
  - **windows-latest runner で実機検証 green (2026-06-14)**: pure-python 52 unit +
    openff-smoke 2 slow が PASS。 検証で判明 — **Windows conda では `openff-toolkit-base`
    必須** (メタパッケージ `openff-toolkit` は AmberTools に hard-depend、 Windows conda
    build が無く solve 不能)。 base は RDKit backend で ambertools 非依存、 ff14SB
    library charges のみの本 route には十分。 docs/platform_support.md + pyproject の
    `[formulation-openff]` に注記反映
  - **小分子電荷を NAGL で Windows native 化 (2026-06-16 green)**: protein は library
    charges で電荷計算不要だが、 小分子 (caprate 等) の AM1-BCC は `sqm` (AmberTools、
    Windows 無し) が要る。 → **`openff-nagl` (ML graph neural net、 pure-Python、 sqm 不要、
    total charge 保存) で代替**。 `molecule_prep` は `NAGLToolkitWrapper` を明示渡しに
    robust 化 (global registry 非自動登録対策)、 smoke に caprate anion → NAGL → Sage
    Interchange を追加 (3 slow PASS)、 CI env + `[formulation-openff]` extra に
    `openff-nagl` + `openff-nagl-models` 追加

### Added — `abmptools.trajectory` (new sub-package)

- cross-platform GROMACS trajectory post-processor (`thin_nojump` / `nojump` / `thin` /
  `wrap_pbc` / `energy`)。subprocess shell=False + pathlib + shutil.which で **Windows
  native** 動作。旧 bash script (wrap_pbc.sh 等) を置換。`python -m abmptools.trajectory`
  CLI + `--gmx` で gmx version 明示 (GROMACS 2026 tpr は古い gmx で読めない gotcha 対応)

### Added (docs)

- `docs/hbond.md` §7 と `sample/hbond/imc_amorphous/README.md` に距離 / 角度分布
  (`distance_dist`) の出力例画像 (IMC amorphous, `imc_hbond_distplot_distance_*`)
  を埋め込み。サンプル参照出力 (CSV / PNG) 自体は別コミットで追加済み。

## [2.1.0] - 2026-06-08

### Added — `abmptools.hbond` 距離 / 角度分布 (`distance_dist.py`)

- 検出済み H-bond の `d(D...A)` / `∠(D-H...A)` を全 record 横断で集計し、
  3 視点のプロット + 統計 CSV を 1 ランで生成:
  - A: `<prefix>_distance_hist.png` — 全 H-bond の `d_DA` 1-D ヒストグラム
    (mean / peak 縦線注釈付き)
  - B: `<prefix>_distance_by_class.png` — クラス別重ね描き (step + 半透明
    fill)。imc mode は `COOH-COOH (dual)` / `COOH-COOH (chain/single)` /
    `COOH-amide` の 3 群、generic mode では `(donor_type, acceptor_type)`
    pair 別。dual 識別は `FunctionalGroupClassification.carboxyl_roles[]
    .dual_partners` を参照
  - C: `<prefix>_distance_angle_2d.png` — `(d_DA, ∠(D-H...A))` 2-D heatmap
    (default 0.05 Å × 5°)
  - `<prefix>_distance_stats.csv` — クラス別 `n / mean / median / std /
    peak / p25 / p75` 表
  - `<prefix>_distance_hist.csv` — long-form (`label, bin_center_DA, count`)
    で再描画用
- 新規 module `abmptools/hbond/distance_dist.py`。`__init__.py` に export 追加
- `AnalyzerConfig` に `do_distance_plots=True` / `distance_d_min=2.0` /
  `distance_d_max=3.6` / `distance_bin_width=0.05` / `angle_bin_width=5.0`
  フィールド追加。CLI に対応する `--no-distance-plots` / `--distance-d-min`
  / `--distance-d-max` / `--distance-bin-width` / `--angle-bin-width` option
- 単体テスト `tests/hbond/test_distance_dist.py` (8 件) 追加。dual ペア
  識別、空入力時の `None` 返却、generic mode pair-type 分け、CSV
  round-trip を回帰固定
- IMC sample (`SI/IMC_result450.0_out_rec900.bdf`、1 record) 実機検証:
  全 81 H-bonds, mean d_DA=2.82 Å, peak 2.72-2.78 Å。`COOH-amide` (N=50)
  が最も短い側に立ち上がる (peak 2.72 Å)、`COOH-COOH (dual)` (N=10) は
  peak 2.82 Å — Yuan 2015 IMC NMR の cyclic dimer / chain / COOH-amide
  振り分けと整合
- `docs/hbond.md` 出力一覧 + CLI option 一覧 + 「7. 距離 / 角度分布」セク
  ション追加

### Fixed — `addsolvfrag` が `AutoFrag='ON'` を出力していた不具合

- `addsolvfrag` は「溶質を手動フラグメント分割した雛形 AJF に、スナップショット
  ごとの溶媒フラグメントを追加する」ツールだが、既定で `aobj.autofrag = True` と
  なっていたため、出力 AJF が `AutoFrag='ON'` かつ `&FRAGMENT` ブロック空という、
  本来の目的と矛盾した内容になっていた（`-ma/--manual` を明示しないと正しい出力に
  ならなかった）。
- 既定を `aobj.autofrag = False` に修正。`-ma` 無しでも出力は `AutoFrag='OFF'` ＋
  雛形フラグメント表に溶媒（HOH/WAT/NA 等）フラグメントを連結した完全な
  `&FRAGMENT` テーブルとなる。`6lu7-covhip` サンプルで NF=311（溶質）→ 1793
  （+1479 水 +3 Na）を確認。
- `-ma/--manual` は既定挙動と同一になったため後方互換目的で残置（no-op）。
- 回帰テスト参照 (`tests/regression/reference/main/addsolvfrag_{covneu,6lu7-covhip}`)
  はバグ挙動を golden master として固定していたため、修正後の正しい出力で再生成。

### Added — `abmptools.trajectory` (new sub-package)

- Cross-platform Python wrapper around `gmx trjconv` / `gmx energy` for
  trajectory post-process (Linux / macOS / Windows、 旧 bash script の置換)
- Public API: `thin_and_nojump`, `nojump`, `thin`, `wrap_pbc`, `gmx_energy`,
  `run_trjconv` (low-level)
- CLI: `python -m abmptools.trajectory {thin_nojump,nojump,thin,wrap_pbc,energy}`
- `aggregation` 系の基本セット (100 ns prod.xtc → `-pbc nojump -skip 10` で
  100 frame、 1 ns stride、 ~300 MB) を Python 1 行で生成
- `subprocess` を `shell=False` + stdin で呼ぶので Windows でも動作、
  `pathlib.Path` で path 区切り差を吸収
- 18 unit tests (`tests/test_trajectory_postprocess.py`)、 gmx subprocess を
  monkeypatch して引数組立 / output 命名 / center group の 2 行 stdin / ndx flag
  / energy term 列挙 / FileNotFoundError / GmxError を検証

### Added (`abmptools.amorphous` cluster center + posres、 branch `fix/cluster-cut-minimum-image-wrap`)

- `BuildConfig.cluster_pdb_path` (str) — pre-built rigid cluster PDB (e.g.
  water trimer の H-bond triangle)。 packmol input に `fixed <center> 0 0 0`
  constraint で box 中央に rigid block 配置、 first component の `n_mol` を
  cluster 分減算 (UDF route の `cluster_file` 等価)。
- `BuildConfig.frozen_atom_indices` (List[int]) — 1-based global GROMACS ndx
  atom indices。 set すると:
  - `system.ndx` に `[ FrozenAtoms ]` index group 追加
  - first moleculetype に `[ position_restraints ]` (`#ifdef POSRES_TRIMER`
    ガード) 追加、 homo pair で件数が cluster mol 数より多い場合は
    `<name>_TRIMER` (cluster 分) + `<name>` (残り) に moletype split。 split
    時 GROMACS の "Only one moltype with [settles] allowed" 制約回避のため
    TRIMER 側 `[ settles ]` を等価な `[ constraints ]` (LINCS-based) に変換
  - 02_nvt_highT 以降の全 mdp に `define = -DPOSRES_TRIMER` (EM は skip、
    packmol overlap relax での numerical 不安定回避)
- `BuildConfig.posres_force_constant` (float, default 10000 kJ/mol/nm²) —
  harmonic restraint の force constant。
- `udf_io.getcontactstructure` の `cutmode='around'` mode — solute (center
  cluster) atoms から `contact_criteria` Å 以内の atom を持つ mol を neighbor
  とする (hydration shell 切り出し)。 既存 legacy mode (`cutmode='contact'`、
  COM 距離 < 2 (r_i + r_j)) と切替可。
- `udf_io.getcontactstructure` の `contact_criteria` default 4.0 Å (1st
  hydration shell)。
- 中央 cluster (mixflag/clusterflag mode) の周辺 27-cell image-neighbor を
  minimum-image wrap する fix — 修正前は中央 cluster の周りに「+L/-L 方向に
  飛び離れた phantom cluster」が現れ FMO 非収束の主因だった。

### Fixed (`fix/cluster-cut-minimum-image-wrap` branch)

- `trajectory_ingest._settles_per_moltype` — TIP3P/TIP4P 水で `[ bonds ]`
  が空 (settles で剛体拘束) の場合、 MDAnalysis Universe に O-H bond が登録
  されず `make_whole` が境界跨ぎ water を unwrap できない問題。 `[ settles ]`
  block から OW を head に O-H bond を synthesize。
- `molcalc.Exportardpos` を obabel 経由から direct Python PDB writer に置換
  — obabel の XYZ→PDB 自動 bond perception が境界跨ぎ broken water の H を
  orphan (`ATOM/UNK` 行) として誤判定する問題を解消。 nameAtom (atom_type) +
  molnames (Mol_Name) を信頼して PDB strict column layout で write。
- `udf_io.moveintocell` の `np.float32` silent-zero bug を回避するための
  cell float cast cherry-pick (`51cd7b3`、 旧 `fix/trajectory-ingest-float-cast`
  branch)。

### Changed

- `sample/formulation/_postprocess/trajectory_thin_nojump.sh` を DEPRECATED 化
  (機能は `abmptools.trajectory` に移行)
- `sample/formulation/octreotide_{l,d}_aggregation_100ns/README.md` に
  Post-process 節を追加 (Python CLI + API の呼び出し例)
- **`amorphous.mdp_protocol.write_wrap_script`** が生成する script を
  `wrap_pbc.sh` (bash) → **`wrap_pbc.py` (Python)** に変更。
  生成 file は `from abmptools.trajectory import wrap_pbc` を呼び、
  Linux / macOS / Windows のどこでも `python wrap_pbc.py` で実行可能。
- **`amorphous.mdp_protocol.write_udf_export_script`** が生成する script を
  `gen_for_udf.sh` (bash) → **`gen_for_udf.py` (Python)** に変更。
  生成 file は `abmptools.trajectory.gmx_energy` + `nojump` を呼ぶ。
- 全 `sample/amorphous/*/README.md` + `run_sample.sh` + `docs/*.md` の
  `bash wrap_pbc.sh` / `bash gen_for_udf.sh` 表記を
  `python wrap_pbc.py` / `python gen_for_udf.py` に統一。
- `getcontactstructure` の `freezegrps` → `position_restraints` (harmonic)
  に switch — freezegrps は LINCS/SETTLE + Domain Decomposition で
  `determinant = -inf` で abort するため。 ~0.04 nm の wiggle のみで cluster
  geometry を維持する (branch `fix/cluster-cut-minimum-image-wrap`)。
- `mdp_protocol.write_all_mdp` の `freeze_group` 引数を `define_posres`
  に rename (同 branch)。

### Added — `formulation.analysis` workflow (Hossain 2023 Fig 1b/1c/2/4 1-コマンド再現)

- `run_analysis(traj, top, out_dir, n_peptides, enhancer_resnames, bile_salt_resnames)` で
  集合体形成 (Fig 1b/1c/Fig 2) + per-residue contacts (Fig 4) + plot を一括実行
- 新規低レベル API: `per_frame_clusters_heavy_atom` (PBC-aware heavy-atom min distance)、
  `compute_per_residue_contacts` (cap 除外 + per-peptide 平均化)
- `aggregate_transition.compute_aggregate_transitions` を atom-index 分割 + heavy-atom
  cutoff + PBC 補正に refactor (GROMACS の resid 1-N リセット慣習に対応)
- CLI: `python -m abmptools.formulation analyze --traj ... --top ... --out ... --n-peptides 6 --enhancer-resnames CPN,CPC --bile-salt-resnames TCH`
- plots: `plot_aggregate_timeseries` / `plot_max_size_distribution` / `plot_per_residue_contacts`
- 出力: `aggregate/{cluster_states,aggregate_size_timeseries}.csv` + `aggregate_summary.json`、
  `contacts/per_residue_contacts.{npy,csv,json}`、 `plots/*.png`
- L 体 10mM 100 ns で smoke 完了: max size mode=6 (49%)、 pct_agg mean=92.8%、
  Phe3/Trp4 が enhancer 主要 contact (10.0/9.8) で論文と整合

### Added (docs)

- **`docs/platform_support.md`** — OS 別 (Linux / macOS / Windows native /
  WSL2) で各 sub-package + 依存 module がどこまで動くかの早見表 + 用途別
  推奨 setup シナリオ + Phase 1 (`formulation` Windows route) 計画。
- `docs/dependencies.md` と `docs/formulation.md` の冒頭に
  `platform_support.md` への誘導 1 行を追加。

### Added (sample)

- `sample/formulation/octreotide_l_aggregation_10mM/` — Hossain 2023 主要系
  (10 mM peptide、 Fig 1-2 出典) 準拠で `n_copies = 6` の設定。 既存
  `octreotide_l_aggregation_100ns/` (peptide × 2 = 3.3 mM) は論文 Table 1 の
  どの系とも一致しない hybrid 構成だったことが判明し、 訂正用として追加。

### Added — `formulation` OpenFF route (Phase 1 完了)

- `force_field_route: "openff"` config 切替で **AmberTools 依存を完全回避**:
  - peptide: `Molecule.from_polymer_pdb` (要 `pdb_path` 入力)
  - small molecule: SMILES → OpenFF Sage 2.x SMIRNOFF
  - charge: peptide=Gasteiger / small mol=AM1-BCC default、 NAGL 切替可
  - typing: amorphous の `create_interchange` を再利用 (Apache-2.0 互換 OpenFF
    + Interchange 経路)
  - GROMACS export: `Interchange.to_gro()` + `.to_top()`
- 新規 module 3 本 (本実装):
  - `formulation/peptide_atomistic_openff.py` (~135 行)
  - `formulation/small_molecule_openff.py` (~105 行)
  - `formulation/topology_openff.py` (~140 行)
- `builder.py` に `_build_openff` ルート (~110 行) 追加、 既存 amber route と
  独立並行、 stage 1/2/4 のみ差替、 stage 3/5/6/7 (packmol / ndx / mdp /
  run_script) は共有
- 実機検証: kggggg × 2 + caprate × 16 + TC × 1 (box 6 nm、 water なし) で
  end-to-end build PASS (707 atoms)、 grompp em PASS
- Amber route 28 unit tests 回帰なし

### TODO (Phase 2)

- **water solvate + ions balance を OpenFF route 内に追加**: 現状 Phase 1 の
  build 出力は dry-mixed system (water なし)。 OpenMM `Modeller.addSolvent`
  経由で TIP3P 充填 + Joung-Cheatham Na/Cl + 0.15 M NaCl を Interchange に
  combine する経路を追加予定。
- **sequence からの peptide build** (PDBFixer 経由): 現状 OpenFF route は
  `PeptideSpec.pdb_path` 必須。 sequence 単独入力からの 3D build 経路を追加。
- **ff14SB SMIRNOFF for natural L-AA peptide**: 現状 whole-peptide Sage 経路
  のみ。 `openff-amber-ff-ports` で peptide 部分だけ ff14SB SMIRNOFF を適用
  する FF stack 経路を追加。
- **Windows native 環境 (実機 Win10/11) での実走確認 + CI 化**。

## [2.0.0] - 2026-05-28

メジャーリリース。v1.15.4 (2026-04-19) 以降に develop で積み上げた 176 commits 分の機能・修正をまとめて公開。

### Breaking changes

- **gro2udf の UDF atom field semantics 変更** (commit `6b32e92`):
  - `Set_of_Molecules.molecule[].atom[].Atom_Name` = **element symbol** (`C` / `H` / `O` 等、旧版は GROMACS atom name `ca1` / `ha0` 等)
  - `Set_of_Molecules.molecule[].atom[].Atom_Type_Name` = system.top atomtype 列、OpenFF SMIRNOFF の `MOL0_<N>` は `<element><N>` (例: `C4` / `H19`) に rewrite
  - `Set_of_Molecules.molecule[].atom[].Atom_ID` = per-molecule **local 0-indexed** (旧 global counter)
  - 1.x との UDF 出力差は test_regression で確認可能
- **abmptools の license 表記更新**: v1.23 以降は Apache-2.0 + NOTICE + CITATION.cff (≤ v1.22.0 は MIT 固定)

### Added (`abmptools.gro2udf` Time/Pressure/Density/Temperature 同期 — v1.x.x 候補)

- `MdpParams` に新規 accessor 追加: `dt`, `nsteps`, `nstxout_compressed`
  (nstxout-compressed / nstxout_compressed / nstxout / nstvout fallback chain),
  `nstenergy`
- `TopModel` に新 field 追加: `dt_ps`, `nsteps`, `nstxout_compressed`,
  `nstenergy` (defaults: 0.001 / 0 / 0 / 0)
- `TopAdapter` で mdp 値を上記 field に格納
- `top_exporter._set_default_condition()` で
  `Simulation_Conditions.Dynamics_Conditions.Time` の以下を mdp から書込み:
    - `delta_T` ([ps] -> [tau] 変換)
    - `Total_Steps`
    - `Output_Interval_Steps` (nstxout-compressed 優先、fallback nstenergy)
- xvg → UDF Statistics_Data mapping を拡張 (旧 `_XVG_TO_UDF_ENERGY` を
  `_XVG_TO_UDF_STATS` に rename + path/unit tuple 形式に):
  - Total Energy / Hamiltonian は既存 mapping にあったが、複数 LJ / Coulomb
    の合算で潰されていた → "Total Energy" -> Energy.Instantaneous.Total
  - 新規追加: Temperature -> Temperature.Instantaneous ([K])
  - 新規追加: Pressure -> Pressure.Instantaneous ([bar])
  - 新規追加: Density -> Density.Instantaneous ([kg/m^3])
  - 新規追加: Volume -> Volume.Instantaneous ([nm^3])
- `_aggregate_energy_per_frame` を `_aggregate_statistics_per_frame` に rename
- `_append_structure` の energy_values 形式を `{(path, unit): value}` に変更
  (旧 `{field: value}`)、Statistics_Data.{Energy, Temperature, Pressure,
  Density, Volume} 全部に対応
- tests/test_mdp_parser.py に 8 件追加 (dt/nsteps/nstxout_compressed/
  nstenergy の default + parsed + fallback chain)、16/16 PASS
- 実機検証 (ketoprofen amorphous):
    delta_T=0.0205 [tau] (= 0.001 ps)
    Total_Steps=500000
    Output_Interval_Steps=5000
    各 frame で Bond/Total/Temperature/Pressure/Density/Volume が record に embed

### Added (`abmptools.gro2udf` multi-frame trajectory + xvg energy — v1.x.x 候補)

- 新規 `abmptools/gro2udf/trajectory_io.py`:
  - `frames_from_multi_gro(path)` — `gmx trjconv -pbc nojump -o output.gro`
    形式の multi-frame `.gro` を pure-Python で parse (stdlib only)
  - `read_xvg(path)` — gmx energy 出力の xvg を `(times, {legend: values})`
    に parse
- `TopExporter.export()` / `.export_model()` に
  `trajectory_path` / `energy_path` / `energy_times` / `energy_series`
  パラメータを追加。CLI で `--trajectory <.gro/.xtc>` / `--energy <.xvg>`
  を指定すると全 frame の Structure record + Statistics_Data.Energy が
  1 UDF に embed される
- xvg → UDF energy field map (`_XVG_TO_UDF_ENERGY` dict):
  Bond / Angle / Proper Dih.+Improper Dih.→Torsion / LJ-14+LJ(SR)+
  Disper.corr.→Nonbonding / Coulomb*→Electrostatic / Potential / Kinetic /
  Total
- frame time と xvg time grid は nearest-neighbour で照合 (501 row vs
  101 frame など denser xvg からも正しく取れる)
- schema 上に energy field が無い cognac 古版では silently skip (try/except)
- docs/gro2udf.md にトラブルシューティングの前 section
  「Multi-frame trajectory + energy を 1 UDF に embed」を新規追加
- tests: `_capture` helper signature を `energy_values=None` 受け入れに更新
  + 既存 6/6 PASS
- 実機検証 (ketoprofen amorphous 50 mol × 33 atom × 101 frame):
    totalRecord=101、各 frame で Cell.a + Bond + Angle 値が正常に embedding

### Added (`abmptools.gro2udf` OCTA8.4 対応 + 単位 fallback — v1.x.x 候補)

- **`abmptools/gro2udf/default_template_cognac101.udf` (新規 bundle)** —
  OCTA8.4 / OCTA8.4 (cognac10.1 までしか同梱されていない環境)
  向けの schema 互換 minimal template。`\include{"cognac101.udf"}` +
  `Unit_Parameter:{"","",1.0,4.184,0.1}` を含み、cognac10.1 でも
  `[nm]` / `[ps]` / `[kJ/mol]` の unit alias が解決される。cognac11.2 default
  template と **出力データは完全同一** (Cell.a, Unit_Parameter.Length 等)。
- `top_exporter._put_with_unit_fallback()` helper を追加 — unit 引数つきの
  put が失敗した場合に unit なしで retry する fallback。`Unit_Parameter` が
  欠けた template でも動くようにする保険。
- `top_exporter._rewrite_cognac_include()` helper + `--cognac-version` CLI
  option を追加。
- CLI: `--cognac-version 101` / `102` 指定時に bundled の cognac10.x template
  を自動選択 (`--template` を渡さなくても OK)。
- `pyproject.toml` の `package-data` を `gro2udf/default_template*.udf` に
  glob 化、両 template を install パッケージに同梱。
- `UDFExportError` の hint メッセージを cognac10.1 環境向けの具体ガイドに
  刷新 (推奨 `--cognac-version 101`、代替 GOURMET Save As + `--template`)。
- docs/gro2udf.md にトラブルシューティング section
  「OCTA8.4 / OCTA8.4 で `file not found:cognac112.udf` エラー」
  を新規追加、`--cognac-version 101` を推奨対処として記載。
- bundled cognac10.1 template の元 (`A20B40A20_in.udf`) はユーザー提供
  (粗視化 sample、`Unit_Parameter` を追加して unit alias を有効化済み)。

### Added (`abmptools.gro2udf` 診断付きエラー — v1.x.x 候補)

- `top_exporter.UDFExportError` 例外クラス + `_section()` context manager を
  追加。`TopExporter.export_model()` の各 stage (`template-copy`,
  `UDFManager-open`, `erase-existing-records`, `Set_of_Molecules`,
  `Structure[record=N]`, `default_condition`, `Molecular_Attributes`,
  `Interactions`) を context manager で wrap して、UDFManager の cryptic な
  例外 (RuntimeError 等) を **どの section / どの template / どの output / 元の
  underlying error / 対処 hint** を含む `UDFExportError` に再 raise する。
- 主な動機: OCTA84 で gro2udf がエラーになった時、UDFManager がどの field
  で落ちたか分からず原因特定が困難。section 名で範囲を絞り、hint で OCTA
  version 差異を案内 (bundled template の利用 / 当該 OCTA から template を
  再生成等)。
- `import UDFManager` 失敗時は `UDFExportError` で OCTA PATH 設定 hint を表示。
- 存在しない template path 渡し時は明示的に `UDFExportError` を pre-flight
  check で投げる。
- docs/gro2udf.md にトラブルシューティング section 追加。
- tests/test_top_exporter_frames_override.py に 2 件追加 (missing template /
  section failure context)、5/5 PASS。

### Added (`abmptools.amorphous` OCTA viewer export — v1.x.x 候補)

- `mdp_protocol.write_udf_export_script(output_dir, ndx, stage, n_energy_terms)`
  新規 API。OpenFF amorphous protocol 後に OCTA viewer (GOURMET) で読み込み可能な
  energy.xvg + nojump gro を生成する `gen_for_udf.sh` を書き出す。
- 内容:
  - `seq <N> | gmx energy -f <stage>.edr -o <stage>_energy.xvg`
    (default N=50、gmx は存在しない term 番号を silently skip するので
    50 で標準 energy term を網羅)
  - `echo 0 | gmx trjconv -f <stage>.{trr,xtc} -s <stage>.tpr -pbc nojump
    -o <stage>_nojump.gro` (`wrap_pbc.sh` の `-pbc mol -ur compact` と
    違い、分子を box 内に wrap せず PBC を跨いで連続的に追跡。
    OCTA viewer (GOURMET) での軌跡再生に適する)
- `FormulationBuilder.build()` から `wrap_pbc.sh` の隣に自動生成。
  返り値 dict に `udf_script` key を追加。
- 既存 sample (`pva_amorphous/`、`ketoprofen*/`) の `md/` 配下にも script
  配置済み。
- ドキュメント: `docs/amorphous.md` (workflow + tree 図 + 仕様説明)
- テスト: `tests/test_mdp_protocol.py::test_write_udf_export_script_*`
  4 件追加 (default / custom stage + n_energy_terms / omits ndx / executable bit)、
  17/17 PASS。
- 実機検証 (ketoprofen amorphous 50 mol × 33 atom):
  `05_npt_final_energy.xvg` (320 KB) + `05_npt_final_nojump.gro` (7.5 MB,
  101 frames × 1650 atoms) 生成 OK。

### Added — `abmptools.formulation` (v1.30.0 候補)

- AA mixed-solution peptide-formulation builder modeled after Hossain et al.
  2023 (*Nanoscale* 15, 19180-19195) — peptide drug + permeation enhancer +
  intestinal bile salt in a cubic water box. **Commercial-permissive force
  fields only**: AMBER ff14SB + GAFF2 + TIP3P + Joung-Cheatham ions
  (CHARMM36 + CGenFF は学術ライセンスのため不採用).
- `FormulationBuilder.build()` 7-stage orchestrator:
  ff_staging → peptide_atomistic (tleap) → small_mol_parameterize
  (acpype GAFF2/AM1-BCC, lazy RDKit SMILES → 3D) → packmol multi-component
  → solvate_ions (tleap solvatebox + parmed) → index (Python resname table)
  → mdp_render → run_script.
- `python -m abmptools.formulation {example, validate, build, analyze, release_us}`
  argparse CLI; `example` emits a smoke JSON, `release_us` writes pull +
  N-window US MDPs for peptide-from-aggregate PMF.
- Analysis stack (opt-in `[formulation-analysis]` extra): aggregate
  transition matrix + cluster timeseries (MDAnalysis + networkx, lazy
  import), per-residue contact map, `gmx dssp` / `gmx sasa` / `gmx hbond`
  wrappers, matplotlib heatmap + timeseries plots.
- Sample configs in `sample/formulation/`: `kggggg_smoke/` (6 nm box,
  ~8k atoms, 200 ps smoke; runs in ~10 min on 4-core CPU) and
  `insulin_smoke/` (10 nm box, 99k atoms, 1 ns; PDB 2G4M downloaded
  on demand from RCSB, not redistributed).
- License posture: acpype (GPL-3.0) and MDAnalysis (GPL-2.0) are mere
  aggregation — acpype via subprocess only, MDAnalysis via lazy import
  inside `analysis/` modules with `pip install abmptools[formulation-analysis]`
  opt-in install (precedent: `amorphous/trajectory_ingest.py`). The
  abmptools wheel itself remains Apache-2.0 with zero GPL code shipped.
- Tests: 51 (models + mdp + ndx) + 19 (packer + small_molecule + peptide_atomistic)
  + 9 (topology) + 5 (builder) + 11 (analysis) + 5 (umbrella_release)
  = 100 unit tests, all mock-based (no gmx/tleap/acpype required). Slow
  integration smoke is gated by `@pytest.mark.slow + skipif(tool missing)`.

### Added (`abmptools.hbond` element + bond-graph fallback — v1.28.0 候補)

- **`fallback_tag_by_element()` 新規追加** (``func_tags.py``): atom_type が
  ``None`` または per-atom unique な値 (OpenFF SMIRNOFF の ``MOL0_X`` 等) の
  atom について、element + bond graph で機能タグを自動付与:
  - O atom: H と bond → ``hydroxyl_O``、H なし → ``carbonyl_O``
  - H atom: O と bond → ``hydroxyl_H``、N と bond → ``amide_H``
  - N atom: C と bond → ``amide_N`` (amide/amine 区別は後段)
  - C atom: ``carbonyl_O`` (= O without H) と bond → ``carbonyl_C``
- ``AnalyzerConfig.use_element_fallback`` フィールド (default True) +
  CLI ``--no-element-fallback`` (strict mode への切替) 追加
- ``functional_groups._tag_atoms_of_mol`` が None tag のみ fallback で補完
  (FF mapping が hit した atom は touch せず、mapping 優先)
- 効果: **OpenFF Sage 経由の UDF を antechamber 経由の GAFF type patch なし
  で直接解析できる**。PVA amorphous sample の手順から antechamber step が
  消える (build_bdf.py 実行 + hbond CLI 直叩きで完結)
- 実機検証: PVA 10-mer × 30 で element fallback ON / antechamber patch なし
  と OFF / antechamber patch ありで完全同一の結果 (rec=0 で 188 H-bonds、
  ratio_donor_busy=62%、ratio_acceptor_busy=60%)
- IMC は GAFF2 type が既に正しく付いているので fallback は no-op (既存
  baseline dual=10 / chain=41 / single=38 / free=36 を維持、57 → 64 tests
  PASS)
- ``tests/hbond/test_element_fallback.py`` 7 件追加:
  alcohol_OH / carboxyl pattern / amide NH / preserves_existing /
  detect_carboxyls_unknown / detect_hydroxyls_unknown / disabled_no_groups

### Added (`abmptools.hbond` generic mode — v1.28.0 候補)

- **`classify_mode={imc, generic}` を新規追加** (`AnalyzerConfig.classify_mode`、
  CLI `--classify-mode`、Jupyter UI dropdown)。default = imc で既存挙動維持
- **`imc` mode** (既存): COOH 中心 4-species (dual/chain/single/free)
- **`generic` mode** (新規): donor-type × acceptor-type の pair 統計 + atom
  単位の role tag (Donor/Acceptor/Both/Candidate)。PVA / peptide / アルコール /
  混合系等、COOH を持たない任意系で動作
- 新規 module `abmptools/hbond/pair_type_stats.py`:
  `PairTypeStat` + `GenericPairClassification` dataclass、`classify_generic()` /
  `summarize_pair_stats()`
- colorizer に generic 版 3 関数: `write_hbond_attributes_generic` /
  `colorize_udf_action_generic` / `write_show_python_script_generic`
- `DEFAULT_GENERIC_COLORS`: Donor=red / Acceptor=cyan / Both=magenta /
  Candidate=faint gray
- 新規出力: `<prefix>_pair_stats.csv` (generic mode のみ)、
  pairs.csv の `kind` 列は generic では `<donor_type>-><acceptor_type>` 形式、
  Attributes 値は `Donor` / `Acceptor` / `Both` (Candidate は skip)
- tests: `test_pair_type_stats.py` 5 件追加 (no_hbonds / single_pair / both_role /
  multiple_pair / unique_dedup)、計 57/57 PASS
- IMC 系での generic mode 動作確認 (`--donor-groups carboxyl
  --acceptor-groups carboxyl_O,amide_O`): carboxyl→amide_O=50, carboxyl→carboxyl_O=31
  が pair_stats に出る (= imc mode の hb_cc=31 + hb_ca=50 と一致)

### Changed (`abmptools.hbond` 4-species classifier + NMR 比較 plot — v1.27.0 候補)

- **per-COOH 分類を 4 species に拡張**: dual / **chain** / single / free。
  Yuan et al. (2015) *Mol. Pharm.* 12, 4518 (DOI 10.1021/acs.molpharmaceut.5b00705)
  の amorphous IMC 13C SSNMR deconvolution (cyclic dimer ~179 / disordered chain
  ~176 / COOH-amide ~172 / free ~170 ppm) に対応する分類軸。
  - **chain** = cyclic dimer ではない COOH-COOH 一方向 H-bond の参加者 (donor or
    acceptor 側)。論文の "disordered chains having various lengths" + chain end
    + ring-larger-than-dimer を一括して捕捉
  - 優先度: dual > chain > single > free (per-COOH も mol 代表 role も同じ)
  - `CarboxylRole.chain_partners` / `FunctionalGroupClassification.n_carboxyls_chain`
    + `ratio_carboxyl_chain` + `MolRole.n_carboxyls_chain` を新規追加
  - summary.csv に `n_carb_chain` / `ratio_carb_chain` / `n_chain_mols` カラム追加
  - colorizer: Mol_Name リネーム経路に `_CHAIN` (Magenta)、action/script 経路に
    `chain: [0.85, 0, 0.85, 1]` (magenta) を追加、write_hbond_attributes の
    value_map にも `chain: Chain` を追加
- **IMC ベースライン値の再変更**: 旧 (3-species) single=49 / free=66 → 新 (4-species)
  chain=41 / single=38 / free=36 (合計 125 = 旧 single+free 115 のうち chain に
  41 が分流)。dual=10、amide accept=49 / free=76 は不変。
- **`sample/hbond/imc_amorphous/plot_nmr_comparison.py` 新規追加**: MD で計算した
  per-COOH 4-species の割合を棒グラフ化。4-species 分類は Yuan 2015 の固体 NMR
  帰属に倣うが、掲載は MD 値のみ (論文の図・表・数値は非再現)。
  <br>※後日更新: 当初は論文 Figure を同梱していたが著作権上除去。MD 値のみの
  `output/imc_hbond_species_fractions.png` を生成する形に改めた。

### Changed (`abmptools.hbond` per-functional-group 統計 — v1.27.0 候補)

- **分類モデルを per-functional-group に変更**: v1.26 まで「分子単位 1 役割」
  だった分類を「官能基単位 (COOH ごと / amide ごと)」に置き換え。1 分子内に
  複数の COOH や amide がある場合に役割が混在するケースに正しく対応する。
  - 新 dataclass: `CarboxylRole`(role=dual/single/free) + `AmideRole`(role=accept/free)
  - 旧 `MolRole` には per-mol の `n_carboxyls_dual/single/free` + `n_amides_accept/free`
    フィールドを追加。`role` 属性は **分子代表 role** (色付け用、優先度
    dual > single > free)
  - `ClassificationResult` は `FunctionalGroupClassification` の alias
    (backward compat、外部 import パス維持)
- **summary.csv 拡張**: per-functional-group カラム (`n_carb_dual`,
  `n_carb_single`, `n_carb_free`, `n_amides`, `n_amide_accept`,
  `n_amide_free`, `ratio_carb_*`, `ratio_amide_*`) を追加。
  従来の `n_dual_mols/n_single_mols/n_free_mols` (mol-level representative)
  も末尾に残す。
- **`<prefix>_classification.csv` 新規追加**: 全 carboxyl / amide ごとの
  (record, group_type, mol_index, group_index, role, partner_count, partners)
  テーブル。NMR の COOH-C / amide-C 信号分離と直接対応する物理量。
- **`<prefix>.bdf` (OCTA viewer プリ描画用コピー) 新規出力**: Mol_Name 維持 (元の
  `molecular` 等) の単純コピー。`<prefix>_colored.bdf` (Mol_Name リネーム済) は
  OCTA viewer のプリ描画で空表示になる問題への対処。CLI option `--no-copy-uncolored`
  でスキップ可。
- **IMC ベースライン値の変更**: 旧 single=73 (COOH→amide H-bond の両当事者を
  カウント) → 新 single=49 (COOH 状態が single の COOH のみカウント、amide
  acceptor 側 mol は COOH free なら free に入る)。`test_imc_regression.py`
  の許容値を更新。物理的には新値が NMR の COOH 信号と直接対応する。
- **verbose log の表現変更**: `COOH dual/single/free=10/49/66 (8%/39%/53%),
  amide accept/free=49/76 (39%/61%)` のように官能基単位の比率を表示。
- **gourmet 可視化手順の docs 改訂**: Python パネルで
  `show.all("line","mol","molname",...)` に書換が必要なこと、OCTA viewer プリ描画
  には `<prefix>.bdf` を使うことを明記。
- **Python action 経路の追加** (`colorize_mode="action"`): GOURMET
  `Draw_Attributes` schema が Mol_Name 維持での per-functional-group
  色付けに対応していないため、`<prefix>_show.act` (`autorun: showHbond()`)
  + `<prefix>_action.bdf` (Mol_Name 維持コピー + Action ヘッダパッチ) を
  併出する経路を追加。各 carboxyl atoms (c/o/oh/ho) と amide atoms (c/o/n)
  を role に応じた色 (dual=red, single=blue, free=gray, accept=cyan) で
  sphere overlay 描画。1 分子内に複数官能基が異なる役割で参加するケースも
  正しく可視化される。CLI `--colorize-mode {molname,action,both}`、default は
  backward compat の `molname`。`colorize_udf_action()` API 新規 export。
- **`<prefix>.bdf` の Attributes[] に hbond タグを append** (OCTA viewer Attribute
  フィルタ対応): 各 functional-group atom (carboxyl c/o/oh/ho、amide c/o/n)
  の `Set_of_Molecules.molecule[].atom[].Attributes[]` 末尾に
  `Name='hbond' Value='Dual'/'Single'/'Free'/'Accept'` を append。既存
  Attributes (`Name='1' Value='molecular:<id>'` 等の OCTA viewer 内部用) は
  touch せず idempotent (再実行で重複 entry 作らない)。OCTA viewer で
  `<prefix>.bdf` を開いた後、Attribute フィルタで `hbond=Dual` 等の atom
  のみ可視化できる (色分けではなくカテゴリフィルタ)。CLI
  `--no-write-attributes` / `--attributes-name NAME` option、default `hbond`。
  `write_hbond_attributes()` API 新規 export、`AnalyzerConfig` に
  `do_write_attributes=True` / `attributes_name='hbond'` フィールド追加。
- **`<prefix>_show.py` 併出** (OCTA viewer (GOURMET) 対応): OCTA viewer (GOURMET) は
  `<prefix>_show.act` の autorun action 形式で落ちることがあるため、同じ描画
  ロジックを autorun ラッパー無しの **プレーン Python script** として
  `<prefix>_show.py` に出力。OCTA viewer で `<prefix>.bdf` (Mol_Name 維持 copy)
  を開いた後、Python パネルから `Load…` → `Run` で同じ役割色 overlay が描画
  される。`write_show_python_script()` API 新規 export。`colorize_mode in
  {action, both}` で .act / .py 両方を併出。

### Added (`abmptools.hbond` 拡張 — v1.26.0 候補)

- **FF 抽象化** (`func_tags.py`): GAFF2/OPLS-AA/CHARMM36/OpenFF の 4 force field
  に対応。atom type → 機能タグ (`carbonyl_C`/`hydroxyl_O`/`amide_N` 等) の
  マッピング辞書 + 自動 FF 検出 (`detect_force_field`) + ユーザ拡張
  (`add_mapping`)。`functional_groups.py` を tag ベースに refactor。
- **任意官能基対選択**: donor `{carboxyl, amide_donor, amine_donor, hydroxyl}` ×
  acceptor `{carboxyl_O, amide_O, hydroxyl_O, ether_O}` から CLI/Python API/
  Jupyter UI で自由に組み合わせ可能。デフォルトは v1.25.0 互換 (COOH→{COOH O,
  amide O})。
- **Secondary amide donor 対応**: `AmideGroup.tert` フラグで tertiary/secondary
  を区別、`detect_amine_donors()` で N-H を donor 集合に取得可能。
  peptide 主鎖 H-bond network の解析が可能になった。
- **Lifetime + autocorrelation** (`lifetime.py`):
  - Continuous lifetime: 連続存在区間の strict 集計
  - Intermittent lifetime: `gap_tolerance` で許容 gap 指定
  - Luzar-Chandler 自己相関 `C(t) = <h(0)h(t)>/<h(0)>` (unbiased estimator)
  - τ_HB (= ∫C(t)dt) 台形則積分
  - Multi-record CLI: `--gap-tolerance N --dt FLOAT --autocorr-max-lag N`
  - 追加出力: `<prefix>_lifetime.csv` + `<prefix>_autocorr.{csv,png}`
- **Jupyter UI 拡張**: donor/acceptor 官能基チェックボックス, lifetime 設定 box,
  FF 自動検出表示, sec amide N-H donor 数表示。
- **テスト 22 件追加** (合計 42、全 PASS):
  - `test_func_tags.py` (9): 各 FF mapping + auto-detect
  - `test_lifetime.py` (8): continuous/intermittent/autocorr/τ_HB
  - `test_amine_donor.py` (5): synthetic N-methylacetamide (GAFF2/CHARMM36) + IMC 否定確認

### Added (`abmptools.hbond` サブパッケージ — v1.25.0 候補)

- **非晶質 MD トラジェクトリ用 H-bond 解析**: COGNAC UDF/BDF を入力に、
  カルボキシル基 (COOH) 同士の dual H-bond (環状二量体) と COOH→アミド C=O の
  single H-bond を幾何条件で検出・分類し、gourmet で 3 色可視化できる UDF
  を出力する新規サブパッケージ。
  - **官能基自動検出** (`functional_groups.py`): GAFF2 atomtype (`c`/`oh`/`ho`/`o`/`n`)
    + bond graph で carboxyl / amide / hydroxyl を検出。Tertiary amide (N-H なし)
    も `tert=True` でマーキング (インドメタシン対応)。
  - **Luzar-Chandler 幾何判定** (`hbond_detector.py`): `d(D-A) ≤ 3.5 Å` かつ
    `∠(D-H-A) ≥ 120°` を default、`strict` モード (`d(H-A) ≤ 2.5 Å` かつ
    `∠ ≥ 150°`) と `custom` モード (任意閾値) も選択可能。直交 cubic box の
    最短像法で PBC 対応。
  - **dual/single/free 分類** (`classifier.py`): 各分子に 3 役割を割り当て
    (優先度 dual > single > free)。Dual は両方向 COOH↔COOH H-bond が成立する
    分子ペアのみ。
  - **gourmet 色付け** (`colorizer.py`): `Set_of_Molecules.molecule[i].Mol_Name`
    を `IMC_DUAL` / `IMC_SINGLE` / `IMC_FREE` にリネームし、
    `Draw_Attributes.Molecule[]` に Red / Blue / Gray の named color を書き込む。
    **GOURMET の Draw_Attributes color は select 型 (9 色名のみ) で RGBA tuple
    は受け付けない** ことを実機検証で確認。`transparency` は 1.0 = 不透明
    (直観に反する) も確認済。
  - **3 経路インターフェース**:
    - CLI: `python -m abmptools.hbond <bdf> --criteria luzar-chandler -o prefix`
    - Python API: `from abmptools.hbond import Analyzer, AnalyzerConfig`
    - Jupyter UI: `open_panel(bdf_path)` (RDKit 2D 構造図 + ipywidgets コンパネ
      + matplotlib count plot のインライン表示)
  - **出力**: per-record summary CSV + H-bond pair CSV + colored BDF + count PNG
- **IMC amorphous サンプル**: `sample/hbond/imc_amorphous/` に CLI/notebook
  ワンライナーと期待値 (T=450 K, 125 IMC: dual=10, single=73, free=42)
- **テスト**: 20 unit + integration テスト
  - 8 角度パターン + PBC wrap geometry の検出器単体テスト
  - IMC count regression (±2-5 mol tolerance)
  - colorize round-trip 検証
- **依存**: `extras_require['hbond'] = ["matplotlib>=3.5"]` 追加。コア機能は
  `numpy` + `UDFManager` (OCTA 同梱) のみ。Jupyter UI は `[jupyter]` +
  `[fragmenter]` (rdkit) を別途併用。
- **ドキュメント**: `docs/hbond.md`、`README.md` (このセクション)

### Added (`abmptools.membrane` sample)

- **`sample/membrane/amber_phaseD/`** — AMBER backend
  (ff19SB + Lipid21 + GAFF2 + TIP3P) 用 umbrella sampling 入力サンプル。
  既存 `sample/membrane/charmm_phaseD/` (CHARMM36 backend) と完全に同じ
  系構成・プロトコル (poly-Ala 5-mer + POPC 32×2 + 0.15 M NaCl、
  13 windows × 1 ns、`window_spacing_nm=0.25`) で力場のみ差し替えた、
  Phase D = L9 verification の AMBER ベースライン側 reference。
  - `input/config_phaseD.json` — `MembraneConfig` JSON
    (`backend="amber"`, `charmm_ff_dir=""`)
  - `README.md` — 概要 / CHARMM 版との差分表 / 結果サマリ
    (PMF +86.7 kJ/mol、CHARMM +97.9 と Δ-11.3 kJ/mol、典型 FF 差) /
    5-stage パイプライン (`packmol-memgen` → `tleap` → `parmed` →
    GROMACS → `gmx wham`) / 実行例 / 重量級結果保管先
  - 商用利用クリーン: `packmol-memgen` + `tleap` + `parmed` 経路で
    CGenFF / CHARMM-GUI 非依存

### Added (Phase B — regression fixture & smoke sample)

  - `R00001` / `R00002` / `R00004` — 3 構造の入力 CIF + 出力
    (`for_abmp/*layer5Zp1-around_ar6.0.{ajf,pdb}`) を fixture 化
  - 共有ドライバ (`input_param`, `segment_data.dat`, `UNK.ajf`) は
    fixture root に同梱
  - 採取は v1.22.0 の `python -m abmptools.readcif -an 32 -l 5` →
    XYZ を pdb 横に copy →
    `python -m abmptools.pdb2fmo -p input_param -xyz` 経由
    (**直接座標 AJF モード**: `&XYZ` block にフル精度浮動小数点で
    座標を埋め込み、PDB 経由の `%8.3f` 切捨を回避)。Phase C 以降の
    リファクタはこの fixture との byte-equivalence を維持すること
  - `R00001layer5Zp1` の `&XYZ` block (832 行 + ヘッダ) は
    `abmptools-sample/.../csp7_ciftest/cifout/layer5/pdb/for_abmp/` の
    既存出力と byte-equivalent を確認済み (namelist 部分のみ
    abmptools 1.22.0 で `&COUPLING/&CIS/&CISGRD/&GF2` 追加・
    `CPFVER=10` skip 化の正常進化差分あり)
  legacy namespace identity)、いずれも `@pytest.mark.slow` 付き。
  `_compare_output_dir` で for_abmp/ 配下を ajf+pdb 両方 byte-compare
  pipeline (`run.sh` で readcif → pdb2fmo を 5 秒で完走、出力は
  fixture と byte-equivalent)。README で Phase C 以降への移行手順も
  併記

### Added (Phase D-3 — numeric reference & verification doc)

Phase D 完了後、abinitmp v2r8 で **csp7 R00001 layer3 HF/6-31G の
本格 FMO 計算**を実機実行 (1 core で 9h 14min 完走) し、Total energy
+ 24 monomer energies + 95 ESP-AOC IFIE 値を numeric reference として
凍結。

  (21 KB) — abinitmp v2r8 reference snapshot:
  - FMO2-HF / 6-31G、24 fragments / 768 atoms / 276 dimer pairs
    (95 ESP-AOC + 181 ESP-PTC)
  - Total energy = -34450.7633976498 hartree
  - 24 monomer HF energies (~-1435.43 hartree each, ±0.005 hartree
    結晶環境による分散)
  - 95 ESP-AOC IFIE pairs (kcal/mol range -9.47 ~ +1.65)
  - tolerances: total/monomer ±1e-5 hartree, IFIE ±1e-3 kcal/mol
  (10 KB) — full log (322 KB) からの必要セクション抜粋。test fixture
  として extract script の roundtrip 検証に使う
  — abinitmp log から JSON reference を生成する CLI ツール。
  再現用に repo に同梱
  62 passed + 1 gated):
  - `test_extract_script_roundtrip` — excerpt log → extract → JSON
    bit-perfect equality (~0.05 s、常時実行)
  - `test_expected_json_shape` — JSON 必須フィールド検証
  - `test_live_layer3_hf_631g_against_reference` —
    `@pytest.mark.slow` + `ENABLE_FMO_LIVE_REGRESSION=1` gate、
    実 abinitmp で 9h 計算 → 凍結 reference と tolerance 比較
    (default skip)
  test matrix (8 ファイル × カバー機能) / 実行コマンド / 検証済み
  環境 (Python 3.11.11 + ase 3.28.0 + abinitmp v2r8) / numeric
  reference summary / failure mode 逆引き表 / 既知の限界
  (HF/6-31G 1 構造のみ、MP2 未取得、abinitmp v2r8 OMP 並列無効) /
  reference 再現 step-by-step

### Notes (Phase D-3)

- 計算実行時の所要時間 (2026-05-07 → 2026-05-08): 9h 14min 33s
  (33253 sec, exit 0) on WSL2 + 1 core。OMP_NUM_THREADS=4 環境
  変数は abinitmp v2r8 build に効かず、MPI 並列前提と判明
  に退避 (重量級の成果物は repo 外へ)、repo には
  excerpt のみ commit
- `getifiepieda` 経由 CSV 化は今回 skip (log の `## HF-IFIE` block を
  Python regex で直接抽出)。CSV 経路は将来追加予定 (※ 2026-05-09 で
  下記 D-3 revised により実装、Python regex 経路は廃止)

### Notes (Phase D-3 revised)

- legacy `expected_layer3_hf_631g.json` (21 KB) と
  `excerpt_layer3_hf_631g.log` (10 KB) は archival のため残置。
  test 経路は読まない (履歴の参照用に残置)
- production (Fugaku csp7 1500 構造、MP2/6-31Gdag) は影響なし —
  anlfmo.py 修正は HF log 用の防御を追加しただけ、MP2 path は不変
- 1 構造 + HF 6-31G という小規模で getifiepieda 経路が成立
  したことで、今後別 system (R00002/R00004 / 他 space group) も
  同じ extract script を流用可能 (zero-padding 桁違いだけ調整)
  2 skipped (numeric live + 別 1 件)、regression なし

## [1.22.0] - 2026-05-06

`[ENERGY] implicit_solvent=GBSA` を使った protein-ligand 単フレーム

### Deferred to v1.22.x / v1.23.x

- per-residue energy decomposition (POC `4_analyse.py` にも非実装)
- MD ensemble averaging (現状は単 frame minimize-only single-point のみ)
- alanine scanning
- 並列 target 処理 (multiprocessing): 現状は逐次のみ
- `forcefield = CHARMM` 経由 (現状 AMBER 想定で kcal/mol)
- explicit solvent (TIP3P box + PBC + PME) 版 MM/PBSA

## [1.21.0] - 2026-05-06

PDB から自動提案する。タンパク質 / DNA は対象外 (既存 `log2config` 経路を維持)、
小分子 / 脂質 / ポリマーに特化した **canonical SMILES グループ化 + C-C 切断
MW walk + Jupyter / CLI 双方の UI** を提供する。リリース予定は v1.21.0。

## [1.20.0] - 2026-05-06

(generalized Replica-Exchange with Solute Tempering — Solute Side Chain
AMBER ff19SB + TIP3P を tleap で勾配し、prmtop+coor を経由して

### Notes

- `param_type` default は `["C", "L"]` (CHARGE + LJ): POC ログ
  (`Setup_Remd_Solute_Tempering>` 出力の `CHARGE=T LJ=T`) と tutorial 12.3
  step3 の SSCR-typical 設定に合わせた。
- 温度ラダー auto mode は v1.20.0 では `geometric` のみ。Patriksson-van der
  Spoel acceptance-ratio formula は v1.20.x で追加予定 (`NotImplementedError`)。
- POC build caveat: icx で `fileio_data_.c` の `ftello64`/`fseeko64` が
  undeclared になる場合は `sed 's/ftello64/ftello/g; s/fseeko64/fseeko/g'`
  で patch。

## [1.19.0] - 2026-05-05

ペプチド-脂質膜系の **PMF (umbrella sampling 経由)** を end-to-end に組める
`abmptools.membrane` (AA umbrella + WHAM) の generic helper を import 経由で
直接再利用する -- コード重複ゼロ。

### External dependencies (new)

- `insane` (PyPI、**GPL-2.0**): Martini bilayer assembly。
  abmptools は subprocess 経由のみ、ソース改変・同梱なし
  (GPL-2.0 接触なし)。`pip install abmptools[cg]` で自動 install。

## [1.18.0] - 2026-05-04

マルチスケール基盤化に向けた CG (粗視化) 系統 (`abmptools/cg/`) の最初の
モジュール。後続で `cg/polymer/` (polyply 経由) や `cg/smallmol/`
(Auto-Martini 経由) を計画している。

### External dependencies (new, optional via `[cg]` extra)

- `vermouth` (PyPI、Apache-2.0): `martinize2` CLI 提供。
- `gmx` (GROMACS, conda): `solvate` / `genion` / `grompp` / `make_ndx` /
  `insert-molecules`。
- `tleap` (AmberTools, conda; **推奨**): atomistic PDB 生成。**不在時は
  extended-backbone fallback** だが、芳香族残基 (W/F/Y) で CG sidechain
  bead が NaN になる可能性あり。研究用途では tleap 推奨。

### Notes

- Pydantic + Click ベースの非公開 m3-peptide リポジトリ (Martini 2/3 +
  ssDNA 対応) からの再実装。DNA 関連 (`external/martinize-dna.py` 等の
  ライセンス未明 derivative) は除外し、`abmptools.amorphous` /
  `abmptools.membrane` と同じ流儀 (`@dataclass` + `argparse`) に統一。

## [1.17.3] - 2026-05-04

CHARMM36 backend 実機検証で見つかった 7 件の互換性 bug を修正。Klauda lab
GROMACS port (`charmm36-feb2026_cgenff-5.0.ff/` / `charmm36-jul2022.ff/`) で
peptide + bilayer + water + ion の smoke build を end-to-end で完走できるように
なった。

### Fixed

- **PDB 4-char residue 名の切り詰め回避** (`parameterize_charmm.translate_pdb_amber_to_charmm`):
  packmol-memgen `--charmm` 出力は POPC / TIP3 / SOD 等を col 18-21 の 4-char で
  書き込む (chain ID 列にはみ出す)。3-char `line[17:20]` fallback で読むと
  POPC → POP に切り詰められて pdb2gmx が `Residue 'POP' not found` で fatal
  終了していた。`line[17:21].strip()` に統一。
- **Klauda port 規約に合わせた residue rename テーブル縮小** (`AMBER_TO_CHARMM_RESNAME`):
  `NME → CT3`, `CYM → CYS`, `HIP → HSP` を削除。Klauda port は `[ NME ]`
  `[ CYM ]` `[ HIP ]` を直接定義しており、CHARMM オリジナルの慣例に
  rename すると `Residue type 'CT3' not found` で fatal。
- **ACE atom 名 mapping を Klauda 規約に修正** (`PER_RESIDUE_ATOM_MAP['ACE']`):
  Klauda port は AMBER 寄りの命名 (CH3/HH31/HH32/HH33/C/O) を採用しており、
  CHARMM オリジナルの CAY/HY1/HY2/HY3/CY/OY ではない。AMBER tleap 由来の
  H1/H2/H3 を HH31/HH32/HH33 にマップ、CH3/C/O はそのまま通す。
- **NME atom 名 mapping を新規追加** (`PER_RESIDUE_ATOM_MAP['NME']`):
  AMBER NME (`N H C H1 H2 H3`) → Klauda port NME (`N HN CH3 HH31 HH32 HH33`)。
  C → CH3、H1/H2/H3 → HH31/HH32/HH33。`SKIP_BACKBONE_RENAME` から `'CT3'`
  を削除し、NME は universal H → HN を適用。これにより `atom C in residue
  NME 7 was not found in rtp entry NME with 6 atoms` 解決。
- **terminus 'None' index の hardcode** (`KNOWN_CHARMM_TERMINUS_NONE` /
  `_resolve_terminus_none_indices`): pdb2gmx `-ter` のメニュー順は ff の
  `.n.tdb` / `.c.tdb` 読み込み順で決まり、`charmm36-feb2026` /
  `charmm36-jul2022` 系では `None = (8, 7)` (N-term: 8、C-term: 7)。
  デフォルトの `0/0` (`NH3+/COO-`) では ACE/NME cap で
  `atom N not found in buiding block 1ACE` で fatal。subprocess による
  動的 probe は pdb2gmx 2021.3 が無効 input で **99% CPU 無限 spin**
  する未解決 bug を踏むため使わず、ff ごとに手動登録した dict で hardcode。
- **TIP3 spurious O-H-H angles の post-process** (`_strip_water_spurious_angles`):
  CHARMM port `solvent.rtp` の `[ TIP3 ]` には rigid 制約用の H1-H2 "bond"
  が並び、pdb2gmx は通常 bond と解釈して 1 water あたり 3 angles
  (1 real + 2 spurious) を生成。`ffbonded.itp` は `HT-OT-HT` のみ定義
  なので、各 spurious angle が `No default U-B types` error 1 件を出す
  (3500 waters で 7000+ errors)。`run_pdb2gmx` の post-process で全
  chain itp の `[ angles ]` セクションをパースし、middle atom が `OT`
  以外の 3-tuple を削除。除去件数を log 出力。
- **subprocess 自動 probe の禁止 (regression 防止)**: pdb2gmx 2021.3 が
  無効 input (`999\n999\n` 等) で 99% CPU 無限 spin する bug を踏まないよう、
  terminus index 解決から subprocess probe を完全削除し、`KNOWN_CHARMM_TERMINUS_NONE`
  辞書 lookup のみに simplify。新 ff port 追加時はメニューを手動で
  1 度確認して dict に登録する運用。

### Tests

- **`tests/test_membrane_charmm_translate.py` の更新と新規テスト**:
  - `test_histidine_tautomers`: HIP は港 (port) 規約で verbatim、HID/HIE
    のみ HSD/HSE に rename
  - `test_protomers`: ASH/GLH/LYN → ASPP/GLUP/LSN、CYM は verbatim
  - `test_caps`: ACE/NME ともに residue 名は verbatim、atom のみ rename
  - `test_ace_atoms`: H1/H2/H3 → HH31/HH32/HH33; CH3/C/O pass through
  - `test_nme_atoms` (新規): C → CH3、H1/H2/H3 → HH31/HH32/HH33、universal
    H → HN
  - `test_4char_residue_not_truncated` (新規): POPC / TIP3 が POP / TIP に
    切り詰められない regression test
  - `test_skip_list_contains_caps_and_ions`: NME が SKIP_BACKBONE_RENAME に
    含まれていないことを assert
  - `test_full_translation_pipeline`: ACE 6 atoms (H1/H2/H3/CH3/C/O) を含む
    fixture に拡張、HH31/HH32/HH33 すべての翻訳を verify
  - 全 26 件 PASS / 0.31 s

### Documentation

- 業界実態を調査してまとめた:
  - ~95% の academic は CHARMM-GUI を使い pdb2gmx を完全回避
  - 少数派は psfgen + charmm2gmx (Wacha & Lemkul 2023 *JCIM*)
  - 完全 commercial-OK の membrane builder は packmol-memgen のみ
  - **本パッケージの修正は世界で数十人しか踏まない pain points**
- AMBER TIP3P と CHARMM-modified TIP3P の H atom LJ 違い (σ=0/ε=0 vs
  σ≈0.4/ε≈0.046 kcal/mol) を memory に記載。abmptools.membrane では
  backend ごとに自動で適切な TIP3P が選ばれる (AMBER: tleap、
  CHARMM36: pdb2gmx の ff dir tip3p.itp)。

### Known limitations

- 実機 MD 完走 (production-leaning run, ~3 時間 GPU) はまだ未検証。
  smoke (grompp pass) までは PASS。
- 長期的には "正攻法 (a)": protein のみ pdb2gmx + water/lipid/ion は
  ff dir の itp を直接 `#include` する方式へ refactor 余地あり (現状は
  業界正攻法の (b) post-process 方式)。

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
