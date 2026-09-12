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

### `tau_p` は `Cell_Mass` からは変換しない

**COGNAC のどのバロスタットかで扱いが変わる。**

| UDF の algorithm | UDF が持つ量 | `tau_p` |
|---|---|---|
| `NPT_Berendsen` | **`tau_P` `[P*ps]`** = 時間 | **変換する** (既定) |
| `NPT_Andersen_Nose_Hoover` | `Cell_Mass` `[mass]` = 質量 | 変換しない。既定 **2.0 ps** |
| `NPT_Parrinello_Rahman_Nose_Hoover` | `Cell_Mass` `[mass]` = 質量 | 同上 |

**Berendsen は変換する。** COGNAC がもともと時間で持っているので、
`tau_P × unit_P × β` で ps に直すだけの単位換算であり、物理を挟まない。

**Nose-Hoover 系は変換しない。** UDF 側が持っているのは質量なので、時間に
直すには運動方程式を経由する必要がある。ここは **選ぶ量**として扱い、既定
**2.0 ps**（Parrinello-Rahman の実用域 2〜5 ps の下端）を書く。`--tau-p` で
上書きできる。

`Cell_Mass` から逆算しないのは 2 つの理由による。

1. **応答時間は系のサイズで決まらない。** 注目する現象と安定性で決めるもの
2. **COGNAC の `Cell_Mass` は経験則。** `CognacSystemUtil.setCellMass` が
   系の全質量を入れているだけで、物理的に最適な値ではない

参考値は INFO ログに出す。Andersen の運動方程式
(`W_int = Cell_Mass · V^(-4/3)`、`W_int·V̈ = ΔP`、`COGNAC1124/src/Anphsystem.cpp`)
を線形化すると

```
tau_p = 2π · √( Cell_Mass · β / V^(1/3) )
```

Parrinello-Rahman は `V^(-4/3)` の補正が無い (`src/PRsystem.cpp`) ので、この式は
目安にとどまる。

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

`Cell_Mass` の単位換算も、J-OCTA は `[mass]` に `[mass·sigma²]` 用の係数
(`unit_Mass · unit_L²`) を掛けている。こちらは掛けない。

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
