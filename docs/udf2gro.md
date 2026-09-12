# udf2gro — COGNAC-UDF to GROMACS Converter

## 概要

`abmptools.udf2gro` は COGNAC が生成した UDF ファイルを GROMACS 形式
（`.gro` / `.top` / `.mdp` / `.ndx`）に変換するパッケージです。

---

## モジュール構成

```
abmptools/abmptools/udf2gro/
├── __init__.py         # Exporter / SystemModel を公開
├── __main__.py         # python -m abmptools.udf2gro 用
├── cli.py              # CLI エントリ（main()）
├── exporter.py         # 薄いオーケストレータ
├── system_model.py     # dataclass 中間表現（Gromacs 文字列を持たない）
├── udf_adapter.py      # UDF → SystemModel 変換（Adapter 層）
└── gromacs/
    └── writers/
        ├── gro_writer.py  # .gro 出力
        ├── top_writer.py  # .top 出力
        ├── mdp_writer.py  # .mdp 出力
        └── itp_writer.py  # .itp（スタブ）
```

### 設計方針（Aプラン）

| 層 | クラス/モジュール | 責務 |
|---|---|---|
| Adapter | `UdfAdapter` | UDFManager → `SystemModel` |
| Model | `SystemModel` (dataclasses) | 純粋データ（Gromacs文字列なし） |
| Writer | `GroWriter`, `TopWriter`, `MdpWriter` | 文字列生成・ファイル書き出し |
| Exporter | `Exporter` | Adapter → Writers の薄い協調役 |
| CLI | `cli.py` | 引数解析・`Exporter` 呼び出し |

---

## ★ COGNAC と GROMACS の対応

同じ物理を、2 つのプログラムが**別の量**で表している。ここを取り違えると
「変換は通るが計算がおかしい」という壊れ方をする。

### 熱浴と圧力浴

| 役割 | COGNAC (UDF) | GROMACS (mdp) |
|---|---|---|
| 熱浴 | **`Q`** — 熱浴の質量 `[mass·sigma²]` | **`tau_t`** — 運動エネルギー振動の**周期** [ps] |
| 圧力浴 | **`Cell_Mass`** — セルの質量 `[mass]` | **`tau_p`** — 時定数 [ps] |

**COGNAC は質量で、GROMACS は時間で指定する。** 質量で持つのは COGNAC だけで、
LAMMPS (`Pdamp`) / AMBER (`taup`) / NAMD (`LangevinPistonPeriod`) / HOOMD (`tau`)
はいずれも時間。

### `Q` → `tau_t` の変換

```
tau_t = 2π · √( Q_d / (g · k_B · T) )

  Q_d = Q · Unit_Parameter.Mass · Unit_Parameter.Length²    [amu nm²]
  g   = 3N        (comm-mode = None)
      = 3N − 3    (comm-mode = Linear)
  k_B = 0.0083144626 amu nm² / (ps² K)
```

**この式は両者のマニュアルに書かれているものを解いただけ。** 導出ではない。

GROMACS リファレンスマニュアル (Nose-Hoover の節、**式 51**):

```
Q = τ_T² · N_f · k · T₀ / (4π²)
```

`τ_T` について解くと `τ_T = 2π·√(Q/(N_f·k·T₀))` で、上の式そのもの。
GROMACS が `Q` ではなく周期を入力させる理由も同じ節にある。

> To maintain the coupling strength, one would have to change Q in proportion
> to the change in reference temperature. For this reason, we prefer to let the
> GROMACS user work with the period τ_T of the oscillations of kinetic energy
> between the system and the reservoir instead.

`.mdp` オプションの説明も同じ:

> for nose-hoover ... **tau-t controls the period of the temperature
> fluctuations at equilibrium**, which is slightly different from a relaxation
> time.

COGNAC 側も同じ形で、`Q` は `g·k_B·T` と組で現れる (COGNAC マニュアル
**式 2.6**、`g`: degrees of freedom):

```
dζ/dt = ( Σ pᵢ²/mᵢ − g·k_B·T ) / Q
```

**`2π` は Nose-Hoover のときだけ。** LAMMPS の `Tdamp` や HOOMD の `tau` は
緩和時間なので付かない。

**`g` の決め方。** UDF の `Dynamics_Conditions.Moment`
(`Calc_Moment` / `Stop_Translation`) から `comm-mode` が決まるので、それに
合わせる。重心運動を除くなら 3 を引く。

### `tau_p` の変換

**どのバロスタットでも UDF から変換する。** 経路が 2 つある。

| UDF の algorithm | UDF が持つ量 | 変換 |
|---|---|---|
| `NPT_Berendsen` | `tau_P` `[P*ps]` = 時間 | 単位換算 (`tau_P × unit_P × β`) |
| `NPT_Andersen_Nose_Hoover` | `Cell_Mass` `[mass]` = 質量 | 運動方程式経由（下記） |
| `NPT_Parrinello_Rahman_Nose_Hoover` | 同上 | 同上 |

`--tau-p` で上書きできる。`Cell_Mass` が無い UDF では 2.0 ps に落ちる。

#### `Cell_Mass` → `tau_p`

```
PR       : tau_p = 2π · √( C · β / (3 · L) )
Andersen : tau_p = 2π · √( C · β / L )        = √3 × PR

C = Cell_Mass · Unit_Parameter.Mass   [amu]     ← sigma² は掛けない
L = max(a, b, c)                      [nm]      ← GROMACS の定義に合わせる
β = .mdp に書く compressibility        [nm³ mol/kJ]
```

`Cell_Mass` はスキーマ上 `[mass]`
(`COGNAC1124/def_udf/cognac*.udf:293` ほか、`Cell_Mass:double [mass]`)。

**両側とも `tau_p` は「箱の振動の周期」なので、周期どうしを等置すれば出る。**

*GROMACS Parrinello-Rahman*:

```
W⁻¹ = 4π²β / (3 · τ_p² · L)      L = 最長の箱要素
b̈  = V W⁻¹ b'⁻¹ (P − P_ref)
```

直方体 `b = diag(a,b,c)` で `ΔP = −ΔV/(Vβ)` と線形化すると
`ä = −4π²Δa/τ_p²`。つまり **`τ_p` そのものが周期**。`L` が最長辺なのは
GROMACS の定義どおりなので、立方体でなくてよい。

*COGNAC Parrinello-Rahman* (`COGNAC1124/src/PRsystem.cpp`):

```cpp
cellMass = cellMassFactor;                                    // l.7
h2 = ((-currentStress - press0Tensor)*hinv*volume)/cellMass;  // l.58
```

同じ線形化で係数を突き合わせると `4π²β/(3τ_p²L) = 1/C`。

*COGNAC Andersen* (`COGNAC1124/src/Anphsystem.cpp`):

```cpp
cellMass = cellMassFactor * pow(volume, -4./3.);   // l.17
aVolume  = (currentPress - pressSum)/cellMass;     // l.122
```

体積を座標にした `W V̈ = ΔP`。`V = L³` として `ΔL` で書き直すと
`ω² = L/(Cβ)` で、PR の `3L/(Cβ)` に対して 1/3。**周期は PR の √3 倍**。
言い換えると、同じ周期を PR で出すには `Cell_Mass` を **3 倍**する。

**実際の圧縮率は式から消える。** GROMACS 側の周期は `τ_p·√(β_true/β_mdp)`
になるので、等置すると `β_true` が両辺で相殺し、`.mdp` に書く `β_mdp` だけが
残る。系の本当の圧縮率を知らなくてよい。

実機 (3050 原子、`Cell_Mass` = 14076.4 amu、セル 2.29911 × 2.29911 × 5.20238 nm):

```
Andersen  tau_p = 8.93 ps
PR        tau_p = 5.16 ps
```

#### 実測で確かめたこと

導出だけでなく、COGNAC と GROMACS を実際に回して確認した (2026-09-12)。

| 事項 | 方法 | 結果 |
|---|---|---|
| Andersen は PR の `√3` 倍 | COGNAC で NPT を回し箱の振動周期を測定 | **1.702** vs √3 = 1.732 |
| `tau_p ∝ √Cell_Mass` | 同上、`Cell_Mass` 4 倍 | **2.010** vs 2.000 |
| GROMACS の `tau_p` = 周期 | 水 2165 分子で `tau_p` を 4→16 ps | **0.99〜1.02** |

周期は体積の自己相関の最初のゼロ点と極小から出す。**Welch 平均した
パワースペクトルのピークは使えない** — 長周期側ほど残留トレンドに引きずられ、
比が窓の取り方で ±15% 動く。

**比だけを根拠にしている。** 絶対値は `β` の精度で決まるが、`β` を体積ゆらぎ
(`⟨δV²⟩ = k_B T V β`) から出すと同じ系でも 2.84e-5〜6.79e-5 bar⁻¹ とばらつく。
減衰のないバロスタットでは体積分散がバロスタット自身の振動に支配され、この式が
成り立たないため。周期比は `β` に依存しないので影響を受けない。

> **`cellMass` は初期体積で固定される。** `Anphsystem.cpp:17` は
> `doInitialStep` の中にあり、`W = Cell_Mass · V₀^(-4/3)` が走り始めの体積で
> 一度決まったきり更新されない。変換対象の UDF のセルは常にその run の初期
> セルなので上の式は成立するが、**GROMACS は `W⁻¹` を現在の箱から毎回計算する**
> ので、箱が大きく変わる run では両者が離れていく。

#### 旧実装との違い

旧実装は次の 2 点でずれていた。

1. **`Cell_Mass` に `Q` 用の `[mass·sigma²]` の換算係数
   (`unit_Mass · unit_L²`) を掛けていた。** スキーマは `[mass]` なので
   `unit_L²` が余計で、all-atom (`unit_L` = 0.1) では 0.01 倍になる。
   **現実的な系では下限 2.0 ps に丸められていた**ので、値としては表面化して
   いなかった
2. **Andersen → PR の質量換算の向きが逆だった** (`1/3` 倍)。上の導出では
   **3 倍**。`τ_p` にすると 3 倍の差。これも 1 に隠れて表面化しない

`L` に最長辺を使う点は変えていない。GROMACS の `W⁻¹` の定義どおりなので。

#### ★ 出てくる値は GROMACS の実務より遅い

**`tau_p` は箱の一辺に比例して伸びる。** `Cell_Mass` は
`CognacSystemUtil.setCellMass` が**系の全質量**を入れるだけの慣用値なので、
密度一定なら `C = ρL³`、`tau_p ∝ √(C/L) = L·√(ρβ)`。

密度 0.85 g/cm³ の立方セルで:

| 辺 [nm] | 原子数の目安 | PR | Andersen |
|---|---|---|---|
| 2 | 890 | 4.5 ps | 7.8 ps |
| 3 | 3,000 | 6.7 ps | 11.7 ps |
| 5 | 14,000 | 11.2 ps | 19.4 ps |
| 8 | 57,000 | 17.9 ps | 31.1 ps |
| 20 | 890,000 | 44.9 ps | 77.7 ps |

**GROMACS で普通に使うのは 2〜5 ps** なので、3 nm を超えるとほぼ常に上回る。

これは間違いではなく、**COGNAC の慣用値が意味しているとおり**の値である。
`W` = 全質量は、箱の振動周期を**音波がセルを横断する時間の 1.6〜2.8 倍**に
する選び方になっている:

```
実機の系   ρ = 0.85 g/cm³、β = 4.5e-5 bar⁻¹ → 音速 c = 1617 m/s
           最長辺 5.2 nm の横断時間 L/c = 3.22 ps
           PR 5.16 ps = 1.6 × L/c    Andersen 8.93 ps = 2.8 × L/c
```

圧力が音響的に均される時間より速く箱を動かさない、という意味では下限として
妥当。ただし**密度の緩和には数 × `tau_p` かかる**ので、短い NPT では箱が
まだ動いている。

`2〜5 ps` を超えたら警告を出す。

### `--tau-p-max` — 上限で打ち切る

```bash
python -m abmptools.udf2gro in.udf out --tau-p-max 5
```

**既定は上限なし。** 変換器が黙って値を変えないため。指定したときだけ丸め、
**元の値と一緒に警告に出す**。

```
tau_p 8.95 ps -> 5.00 ps (--tau-p-max). Cell_Mass asks for 8.95 ps;
this is a deliberate cap, not the value the UDF describes.
```

`--tau-p 5` との違いは、**小さい系を引き延ばさない**こと。上限より小さい値は
そのまま通る。`--tau-p` を同時に指定した場合はそちらが勝つ（直接の指示なので）。

> **どちらの経路かで方針が違う。** `udf2gro` は**既にある UDF を変換する**ので
> 既定では忠実に写す。`amorphous` は**プロトコルを組み立てる**ので実用値
> (`tau_p` = 2.0 ps) を選ぶ。GROMACS で走らせるのが目的なら、`udf2gro` にも
> `--tau-p-max 5` を付けるとよい。

#### 注意

- **`pcoupl = MTTK` では近似。** Andersen の既定の行き先だが、MTTK は
  バロスタット質量の定義が違う。厳密に合わせたいなら
  `--barostat Parrinello-Rahman` にする
- **`--barostat C-rescale` は GROMACS 2021 以降が要る。** 2020 系に渡すと
  `Invalid enum 'C-rescale' for variable pcoupl` で `grompp` が止まる。
  2020 系で使える `pcoupl` は No / Berendsen / Parrinello-Rahman /
  Isotropic / MTTK
- **`tau_p ≥ 2·tau_t` を書き出し時に確認する。** `grompp` も言うが、こちらは
  両方の値を持っているので先に言える

### アルゴリズムの対応

COGNAC は `NPT_<バロスタット>_<サーモスタット>` という命名。

| COGNAC | `tcoupl` | `pcoupl` | `pcoupltype` |
|---|---|---|---|
| `NVE` | no | no | — |
| `NVT_Nose_Hoover` | nose-hoover | no | — |
| `NVT_Berendsen` | berendsen | no | — |
| `NPT_Andersen_Nose_Hoover` | nose-hoover | **MTTK** | isotropic |
| `NPT_Parrinello_Rahman_Nose_Hoover` | nose-hoover | Parrinello-Rahman | **anisotropic** |
| `NPT_Berendsen` | berendsen | berendsen | isotropic |

**Andersen と Parrinello-Rahman はどちらもバロスタット**で、動かす自由度が違う
(体積のみ / 体積 + セルの形)。**Nose-Hoover はサーモスタット**。**Berendsen は
温度用と圧力用の両方がある**ので、名前だけでは判断できない。

> **等方圧縮には `NPT_Andersen_Nose_Hoover`。** Parrinello-Rahman を選ぶと
> `pcoupltype = anisotropic` になり、非晶質バルクでは箱が不必要に歪む。
> ただし MTTK は GROMACS で拘束 (LINCS / SETTLE) と併用できないので、
> 水を含む系では `.mdp` を手で `Parrinello-Rahman` + `isotropic`、または
> `C-rescale` に直すこと。

## 使い方

### CLI
```bash
python -m abmptools.udf2gro input.udf output_prefix
```

### Python API

```python
from abmptools.udf2gro import Exporter

Exporter().export("system.udf", "output")
# → output.gro, output.top, output.mdp が生成される
```

---

## テスト手順

`udf2gro/test/` にリファレンス入出力が用意されています。

```bash
cd sample/udf2gro/input

# 1. 変換実行（出力はカレントディレクトリ）
python -m abmptools.udf2gro test.udf test

# 2. 期待値と比較
diff test.top  ../output/test.top
diff test.gro  ../output/test.gro
diff test.mdp  ../output/test.mdp

# 3. 全て差分なしであれば PASS
```


---

## 対応力場・アンサンブル

### 力場（Force Field）

| FF_TYPE | 名称 |
|---|---|
| 1 | AMBER |
| 2 | GAFF |
| 3 | GAFF2 |
| 4 | DREIDING |
| 5 | OPLS |
| 6 | OPLSUA |

### アンサンブル

| UDF アルゴリズム | GROMACS 設定 |
|---|---|
| NVE | integrator=md-vv, tcoupl=no, pcoupl=no |
| NVT_Nose_Hoover | integrator=md-vv, tcoupl=nose-hoover |
| NVT_Berendsen | integrator=md-vv, tcoupl=berendsen |
| NPT_Andersen_Nose_Hoover | integrator=md-vv, tcoupl=nose-hoover, pcoupl=MTTK |
| NPT_Parrinello_Rahman_Nose_Hoover | integrator=md, tcoupl=nose-hoover, pcoupl=Parrinello-Rahman |
| NPT_Berendsen | integrator=md, tcoupl=berendsen, pcoupl=berendsen |
| Kremer_Grest | integrator=sd |

---

## 変更点サマリ

### 変更ファイル一覧

| ファイル | 変更内容 | 理由 |
|---|---|---|
| `udf2gro/udf2gro.py` | `main()` を `abmptools.udf2gro.Exporter` への委譲に変更 | CLI 互換を維持しつつ実装を移管 |
| `abmptools/abmptools/__init__.py` | `Udf2groExporter` の条件付きインポートを追加 | パッケージ統合（UDFManager 未インストール時は無視） |
| `abmptools/docs/udf2gro.md` | 本ドキュメントを新規作成 | テスト手順・設計説明 |

### 新規ファイル一覧

`abmptools/abmptools/udf2gro/` 以下の全ファイル（12ファイル）が新規作成。
元の `exportGromacs()` の責務を以下に分割：

- **UDF読み取り** → `udf_adapter.py`
- **GRO書き出し** → `gromacs/writers/gro_writer.py`
- **TOP書き出し** → `gromacs/writers/top_writer.py`
- **MDP書き出し** → `gromacs/writers/mdp_writer.py`
- **データ定義（各種 dataclass）** → `system_model.py`
- **協調・NDX書き出し** → `exporter.py`
