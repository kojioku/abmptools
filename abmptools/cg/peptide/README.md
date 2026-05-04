# abmptools.cg.peptide

Martini 3 ペプチドアモルファス系の end-to-end ビルダー。残基配列 (1-letter)・
本数・box サイズを JSON (or YAML) で指定すると、GROMACS MD 入力一式 (CG `.itp`,
`.gro`, `topol.top`, `index.ndx`, `mdp/{em,nvt,npt,md}.mdp`, `run/run.sh`)
を出力する。

## パイプライン

```
PeptideBuildConfig (JSON/YAML)
  ↓ (1) atomistic PDB         tleap (推奨) or extended-backbone fallback
  ↓ (2) martinize2 -ff martini3001
  ↓ (3) gmx insert-molecules  -> packed.gro
  ↓ (4) gmx solvate           -> system_solv.gro (Martini W)
  ↓ (5) gmx grompp + genion   -> system_ions.gro (NaCl 中和 + 0.15 M)
  ↓ (6) MDP + index + run.sh
```

## クイックスタート

```bash
# 1. CLI で example config を吐く
python -m abmptools.cg.peptide example > kgg.json

# 2. Martini 3 force field を別途取得 (パッケージ未同梱)
#    cgmartini.nl は zip 一括配布。必要 3 itp を ff/ に展開する:
mkdir ff
curl -L -o /tmp/m300.zip \
    https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip
unzip -o -j /tmp/m300.zip \
    "martini_v300/martini_v3.0.0.itp" \
    "martini_v300/martini_v3.0.0_solvents_v1.itp" \
    "martini_v300/martini_v3.0.0_ions_v1.itp" \
    -d ff/

# 任意: solvent_enabled=True で gmx solvate する場合は
# Martini 3 water box (martini_v3.0.0_water.gro) も必要。
# cgmartini.nl は water.gro を直接配布していないので、
# gmx insert-molecules で自作するか、Martini 3 tutorial archive
# (例: martinize2 sample, lipid tutorial) に同梱されているものを利用する。

# 3. 依存と .itp の存在確認
python -m abmptools.cg.peptide validate --config kgg.json --ff-dir ./ff

# 4. ビルド
python -m abmptools.cg.peptide build --config kgg.json --ff-dir ./ff -o ./out

# 5. GROMACS で実行
bash out/run/run.sh
```

## 外部ツール (本パッケージ未同梱)

| ツール | 必須/任意 | 用途 | 取得方法 |
|---|---|---|---|
| `martinize2` (vermouth, Apache-2.0) | 必須 | M3 CG mapping | `pip install abmptools[cg]` (PyPI 経由 vermouth) |
| `gmx` (GROMACS) | 必須 | solvate / genion / grompp / make_ndx / mdrun | `mamba install -c conda-forge gromacs` |
| `tleap` (AmberTools) | **推奨** (任意) | atomistic PDB | `mamba install -c conda-forge ambertools` |
| Martini 3 ITP 3 ファイル (`martini_v3.0.0.itp` / `_solvents_v1` / `_ions_v1`) | 必須 | force field | `martini_v300.zip` から取得 |
| Martini 3 water box `.gro` | 任意 (`solvent_enabled=True` のみ) | solvate template | cgmartini 直接配布なし。`gmx insert-molecules` で自作 / tutorial archive 流用 |

## tleap fallback の注意

`tleap` が PATH にない場合、本パッケージは python だけで beta-strand 風の伸びきり
backbone を書き出して martinize2 に渡す。ただし sidechain heavy atom が無いため、
**芳香族残基 (W / F / Y) で martinize2 が CG sidechain bead 座標を NaN にする
可能性**がある (extended-backbone fallback の構造的限界)。

研究品質では tleap の使用を推奨する。`validate` サブコマンドが tleap の有無を
チェックして警告を出す。

## 入力 schema 概要

`PeptideBuildConfig` (`@dataclass`、`to_json/from_json` で JSON 往復可):

```json
{
  "peptides": [
    {"name": "kgg", "sequence": "KGG", "count": 5}
  ],
  "box_size_nm": 8.0,
  "box_lengths_nm": null,
  "temperature": 310.0,
  "solvent_enabled": true,
  "neutralize": true,
  "nacl_molar": 0.15,
  "martini_itp_dir": "./ff",
  "martinize2_path": "martinize2",
  "gmx_path": "gmx",
  "tleap_path": "tleap",
  "output_dir": ".",
  "seed": null,
  "mdp_em": true, "mdp_nvt": true, "mdp_npt": true, "mdp_md": true,
  "em_steps": 50000,
  "nvt_nsteps": 250000,
  "npt_nsteps": 250000,
  "md_nsteps": 5000000,
  "dt_fs": 20.0
}
```

YAML 入力は PyYAML の有無で動的に判定 (extras `[cg]` に含まれる)。

## 出力レイアウト

```
output_dir/
├── config.json              # 入力 config の dump (再現用)
├── topol.top
├── packed.gro
├── system_solv.gro
├── system_ions.gro          # 最終座標 (build() 返値の "gro")
├── ions.{mdp,tpr}
├── index.ndx
├── molecules/
│   └── <name>/
│       ├── <name>_atomistic.pdb
│       ├── <name>_cg.pdb
│       └── <name>.itp
├── mdp/
│   ├── em.mdp / nvt.mdp / npt.mdp / md.mdp
└── run/
    └── run.sh                # bash run/run.sh で MD 4 段階を順次実行
```

## ライセンス・権利

- 本サブパッケージ自体: abmptools 全体と同じ MIT License
- `martinize2` (`vermouth-martinize`): Apache-2.0 (本パッケージは subprocess で
  呼ぶのみ、コードもデータも同梱せず改変なし)
- Martini 3 force field `.itp` (cgmartini.nl 配布): 本パッケージは**未同梱**。
  ユーザーが各自取得すること

## abmptools 内での位置付け

`abmptools/cg/` は MO-AAMD-CGMD マルチスケール基盤の **CG (粗視化) 系統**として
新設された namespace。全原子系の `abmptools/amorphous` (OpenFF/GAFF) や
`abmptools/membrane` (CHARMM36) と同じ流儀 (`@dataclass` + `argparse` +
`__main__.py` entry) で書かれている。
