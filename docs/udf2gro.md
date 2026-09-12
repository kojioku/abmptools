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

**`2π` が要る理由。** COGNAC の `Q = g k_B T τ²` の `τ` は応答時間だが、
GROMACS の `tau_t` は Nose-Hoover では**振動の周期**で、`2π·τ` にあたる。
LAMMPS や HOOMD は緩和時間を取るので `2π` は付けない。

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

#### 旧実装との違い

`Export_GROMACS.py` の移植だったので、次の 2 点を引き継いでいた。

1. **`Cell_Mass` に `Q` 用の `[mass·sigma²]` の換算係数
   (`unit_Mass · unit_L²`) を掛けていた。** スキーマは `[mass]` なので
   `unit_L²` が余計で、all-atom (`unit_L` = 0.1) では 0.01 倍になる。
   **現実的な系では下限 2.0 ps に丸められていた**ので、値としては表面化して
   いなかった
2. **Andersen → PR の質量換算の向きが逆。** `Export_GROMACS.py` は
   `W = W * 1.0/3.0` (コメント `Andersen -> parrinello_Rahman`) としているが、
   上の導出では **3 倍**。`τ_p` にすると 3 倍の差。これも 1 に隠れて
   表面化しない

なお `max_L` を使うのは `Export_GROMACS.py` と同じで、これは GROMACS の
`W⁻¹` の定義どおりなので正しい。

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

`2〜5 ps` を超えたら警告を出す。**生産計算で GROMACS の慣行に合わせたいなら
`--tau-p 2` で決め打ちする**のが素直。

#### 注意

- **`pcoupl = MTTK` では近似。** Andersen の既定の行き先だが、MTTK は
  バロスタット質量の定義が違う。厳密に合わせたいなら
  `--barostat Parrinello-Rahman` にする
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

### ★ J-OCTA の `Export_GROMACS.py` とは値が違う

同じ UDF でも `tau_t` が一致しない。**こちらの式が上記のとおりで、J-OCTA は
分母の `g · k_B` を落としている**ため。

| 原子数 | `abmptools.udf2gro` | J-OCTA `Export_GROMACS.py` |
|---|---|---|
| 110 | **0.6283 ps** | 1.0408 ps |
| 3050 | **0.6282 ps** | 5.4803 ps |

J-OCTA 側は `tau_t ∝ √(3N)` で系のサイズとともに増大する。同じ J-OCTA でも
`Export_LAMMPS.py` / `Export_HOOMD_blue.py` は `3N` で割っており、こちらは
そちらと同じ扱いにしている。詳細は `SI/udfcheck/jocta_tau_report_20260912.md`。

`tau_p` も一致しない。J-OCTA の式の骨格 (`√(4π²βW/(3L))`) は GROMACS の
`W⁻¹` の定義から出ていて正しいが、`[mass]` の `Cell_Mass` に `[mass·sigma²]`
用の係数 (`unit_Mass · unit_L²`) を掛けているため 0.01 倍になり、下限 2.0 ps
に丸められる。Andersen → PR の質量換算も向きが逆に見える (`1/3` 対 `3`)。

| 原子数 | `abmptools.udf2gro` | J-OCTA `Export_GROMACS.py` |
|---|---|---|
| 3050 (Andersen) | **8.93 ps** | 2.0 ps (下限に丸め) |
| 3050 (PR) | **5.16 ps** | 2.0 ps (下限に丸め) |

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
