# abmptools.trajectory — GROMACS trajectory 後処理ラッパー

`gmx trjconv` / `gmx energy` を **クロスプラットフォーム(Windows native 対応)の
Python API + CLI** で叩く後処理ユーティリティ。従来 sample / builder に散在していた
bash script(`trajectory_thin_nojump.sh`、amorphous の `wrap_pbc.sh`、`gen_for_udf.sh`
等)を Python に統一したもの。GROMACS 本体は同梱せず `subprocess` で呼ぶだけ。

## CLI

```
python -m abmptools.trajectory <subcommand> --traj <xtc> --tpr <tpr> [options]
```

| subcommand | 内容 | 対応する gmx |
|---|---|---|
| `thin_nojump` | `--skip N` で間引き + `-pbc nojump`(分子を割らない連続座標) | `trjconv -skip N -pbc nojump` |
| `nojump` | `-pbc nojump` のみ | `trjconv -pbc nojump` |
| `thin` | `--skip N` で間引きのみ | `trjconv -skip N` |
| `wrap` | `-pbc mol -ur compact`(+ 任意 `--center`)。VMD 向け compact 表示 | `trjconv -pbc mol -ur compact` |
| `energy` | `.edr` → `.xvg`(エネルギー項抽出) | `energy` |

### 共通オプション(座標系サブコマンド)

| オプション | 既定 | 内容 |
|---|---|---|
| `--traj` | (必須) | 入力 trajectory(`.xtc` 等) |
| `--tpr` | (必須) | reference 構造(`.tpr` 推奨) |
| `--out` | 自動命名 | 出力パス(省略時は `<traj>_<op>_skip<N>.xtc` 等) |
| `--group` | `System` | trjconv の出力グループ |
| `--ndx` | なし | index file |
| `--gmx` | `gmx` | 使用する gmx 実行ファイル(下記 gotcha 参照) |

`energy` サブコマンドは `--edr` / `--out`(`.xvg`)/ `--terms-max`(既定 50)/ `--gmx`。

### 例

```bash
# 間引き + nojump(prod.xtc を 1/10 に、分子を割らずに連続化)
python -m abmptools.trajectory thin_nojump --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10

# VMD 向け compact wrap
python -m abmptools.trajectory wrap --traj 05_npt_final.xtc --tpr 05_npt_final.tpr \
    --out 05_npt_final_pbc.xtc --ur compact

# エネルギー項を xvg に
python -m abmptools.trajectory energy --edr prod/prod.edr --out prod_energy.xvg
```

## Python API

```python
from abmptools.trajectory import (
    thin_and_nojump, nojump, thin, wrap_pbc, gmx_energy, run_trjconv, GmxError,
)

out = thin_and_nojump(trajectory="prod/prod.xtc", tpr="prod/prod.tpr", skip=10)
# -> PosixPath('prod/prod_nojump_skip10.xtc')
```

| 関数 | 内容 |
|---|---|
| `thin_and_nojump(trajectory, tpr, skip=..., ...)` | 間引き + nojump。出力 Path を返す |
| `nojump(...)` / `thin(...)` | それぞれ単独 |
| `wrap_pbc(..., ur="compact", center=None)` | `-pbc mol -ur compact` |
| `gmx_energy(edr, out, terms_max=50, gmx="gmx")` | `.edr` → `.xvg` |
| `run_trjconv(...)` | 低レベル `trjconv` ラッパー(任意フラグ) |
| `GmxError` | gmx 実行失敗時に送出される例外 |

## 注意 (gotcha)

- **tpr の版と gmx の版が合わないと失敗する**。新しい GROMACS で書いた `.tpr`
  (例: 2026 系の tpr v138)は古い `gmx` では読めない。使いたい `gmx` を
  **`--gmx /path/to/gmx` で明示**して、tpr を書いた版と揃える。
- グループ選択は `--group`(既定 `System`)。溶質だけ等にしたい場合は `--ndx` +
  グループ名を指定。

## 関連

- 上流の MD 生成: [`amorphous.md`](./amorphous.md)
  (どちらも `wrap_pbc` 相当の後処理を内部で使う)。
- 生成した `.xtc` を可視化・解析へ: [`gro2udf.md`](./gro2udf.md)(COGNAC UDF 化)。
