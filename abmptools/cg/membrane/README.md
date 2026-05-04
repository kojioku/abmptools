# abmptools.cg.membrane

Martini 3 ペプチド-脂質膜系の **PMF (umbrella sampling 経由)** を end-to-end に
組む CG 系統 2 番目のサブパッケージ。`abmptools.cg.peptide` を内部 sub-call
してペプチド CG を生成し、`insane` (GPL-2.0) で POPC bilayer に埋め込み、
13 windows × 1 ns × k=1000 kJ/mol/nm² の umbrella sampling 用 MDP を出力する。

## パイプライン (7 stage)

```
MembraneCGBuildConfig (JSON / YAML)
  ↓ (0) copy 4 Martini 3 ITPs            -> output_dir/
  ↓ (1) cg.peptide sub-call               -> molecules/<name>/{name}_cg.pdb + .itp
  ↓ (2) insane                            -> bilayer.gro + insane_topol.top
  ↓ (3) topology composer                 -> topol.top (4 ITP includes + peptide ITP)
  ↓ (4) normalize NA+/CL- in .gro         -> system_ions.gro
  ↓ (5) write_ndx_from_gro_cg             -> index.ndx (Bilayer/Peptide/W/NA/CL/Non_Bilayer)
  ↓ (6) MDP生成 (em/nvt/npt/pull/windows) -> mdp/, pull/, windows/win_NNN/
  ↓ (7) run.sh                            -> run.sh
```

## クイックスタート

```bash
# 1. Sample config を取得
python -m abmptools.cg.membrane example > kgg_popc.json

# 2. Martini 3 force field を別途取得 (パッケージ未同梱)
mkdir ff
curl -L -o /tmp/m300.zip \
    https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip
unzip -o -j /tmp/m300.zip \
    "martini_v300/martini_v3.0.0.itp" \
    "martini_v300/martini_v3.0.0_solvents_v1.itp" \
    "martini_v300/martini_v3.0.0_ions_v1.itp" \
    "martini_v300/martini_v3.0.0_phospholipids_v1.itp" \
    -d ff/

# 3. insane を install (GPL-2.0、subprocess 呼び出しのみ)
pip install insane

# 4. 依存とファイル存在確認
python -m abmptools.cg.membrane validate --config kgg_popc.json --ff-dir ./ff

# 5. ビルド
python -m abmptools.cg.membrane build --config kgg_popc.json --ff-dir ./ff -o ./out

# 6. GROMACS で実行 (em -> nvt -> npt -> pull -> windows -> wham)
bash out/run.sh
```

## 外部ツール (本パッケージ未同梱)

| ツール | 必須/任意 | ライセンス | 用途 | 取得方法 |
|---|---|---|---|---|
| `insane` | 必須 | **GPL-2.0** | Martini bilayer 構築 | `pip install insane` |
| `martinize2` (vermouth) | 必須 | Apache-2.0 | M3 CG mapping (cg.peptide 経由) | `pip install abmptools[cg]` |
| `gmx` (GROMACS) | 必須 | LGPL | grompp / mdrun / wham | `mamba install -c conda-forge gromacs` |
| `tleap` (AmberTools) | **推奨** (任意) | LGPL | atomistic PDB | `mamba install -c conda-forge ambertools` |
| Martini 3 ITPs (4 files) | 必須 | (cgmartini.nl) | force field | `martini_v300.zip` から抽出 |

## ライセンスと配布

- 本サブパッケージ自体: abmptools 全体と同じ MIT License
- `insane` (GPL-2.0): **subprocess で呼ぶのみ、コードもデータも同梱せず改変なし** -- GPL FAQ における "mere aggregation" 扱い、abmptools の MIT 性は保たれる
- `martinize2` (vermouth-martinize、Apache-2.0): cg.peptide 経由で同様に subprocess only
- Martini 3 force field `.itp` (cgmartini.nl 配布): **本パッケージ未同梱**、ユーザーが各自取得

## Insane の挙動と統合課題 (内部処理が解決済み)

Insane の出力をそのままでは GROMACS で動かないため、`topology_composer` が
以下を自動で post-process する:

| 課題 | 解決 |
|---|---|
| `#include "martini.itp"` (generic 名) | M3 4-ITP の include 群に書き換え |
| `[ molecules ]` で `Protein` (placeholder) | peptide ITP の実 moleculetype 名 (cg.peptide は `molecule_0`) で置換 |
| `NA+` / `CL-` (atom name & resname) | `NA` / `CL` に正規化 (gro + topol.top 両方) |

## 入力 schema 概要

`MembraneCGBuildConfig` (`@dataclass`、`to_json/from_json` で JSON 往復可):

```json
{
  "lipids": [{"resname": "POPC", "n_per_leaflet": 64}],
  "peptide": {"name": "kgg", "sequence": "KGG", "initial_z_offset_nm": 3.0},
  "insane_d_nm": 8.0,
  "box_z_nm": 14.0,
  "nacl_molar": 0.15,
  "equilibration": {"em_steps": 5000, "nvt_nsteps": 25000, "npt_nsteps": 250000, "temperature_K": 310.0, "tau_p_ps": 12.0},
  "pulling":       {"pull_rate_nm_per_ps": 0.001, "pull_force_constant": 1000.0, "nsteps": 250000},
  "umbrella":      {"z_min_nm": -1.5, "z_max_nm": 1.5, "window_spacing_nm": 0.25, "force_constant_kj_mol_nm2": 1000.0, "window_nsteps": 50000}
}
```

YAML 入力は PyYAML の有無で動的に判定 (extras `[cg]` に含まれる)。

## 出力レイアウト

```
output_dir/
├── input/config.json                    # config dump (再現用)
├── martini_v3.0.0*.itp (4 files)        # ff_dir からコピー
├── topol.top                            # 完成 topology
├── bilayer.gro                          # raw insane 出力 (audit用)
├── insane_topol.top                     # raw insane topology (audit用)
├── system_ions.gro                      # 最終座標 (build() 返値の "gro")
├── index.ndx                            # Bilayer / Peptide / W / NA / CL / Non_Bilayer
├── molecules/<name>/
│   ├── <name>_atomistic.pdb
│   ├── <name>_cg.pdb
│   └── <name>.itp
├── mdp/
│   ├── em.mdp / nvt.mdp / npt.mdp       # equilibration
├── pull/pull.mdp                        # constant-rate steered MD (NVT chassis)
├── windows/win_000/window.mdp .. win_012/window.mdp  # 13 static umbrella
└── run.sh                               # 全 stage 実行 (em→nvt→npt→pull→windows→wham)
```

## abmptools 内での位置付け

CG 系統 (`abmptools/cg/`) の 2 番目のモジュール。先行 `abmptools.cg.peptide`
(1.18.0) を内部で sub-call する形で再利用、コード重複ゼロ。AA 系の
`abmptools.membrane` (umbrella sampling + WHAM) と並走する CG 版。
WHAM 解析は `abmptools.membrane.pmf.run_wham` を duck-typing で直接再利用。

## CG umbrella protocol の default 値

| 項目 | 値 | 由来 |
|---|---|---|
| 反応座標 | z 軸 (Bilayer COM → Peptide COM) | Martini AA membrane と同じ |
| z 範囲 | -1.5 〜 +1.5 nm | 13 windows |
| 窓間隔 | 0.25 nm | AA membrane と同じ (CG dt=20 fs に scale) |
| 窓数 | 13 | -1.5 .. +1.5 step 0.25 inclusive |
| force constant | 1000 kJ/mol/nm² | AA membrane と同じ |
| 窓あたり | 1 ns (50,000 steps × 20 fs) | AA は 5 ns at 2 fs ≒ 同じ wall time |
| pulling rate | 0.001 nm/ps (= 1 nm/ns) | AA と同じ |
| pulling 時間 | 5 ns (250,000 steps × 20 fs) | AA は 10 ns at 2 fs |
| pulling geometry | `direction-periodic` (NVT, no Pcoupl) | 大変位対応 |
| 窓 geometry | `direction` (NPT semiisotropic) | 動的 box 互換 |
| 温度 | 310 K (V-rescale, Bilayer / Non_Bilayer 2-group) | 生体系の習慣 |
| 圧力 | 1 bar (c-rescale, semiisotropic, tau_p=12 ps) | Martini 標準 |
| 圧縮率 | 3e-4 bar⁻¹ (xy + z 共通) | Martini 標準 (水の 6×) |

## 後続予定 (v1.20+)

- 複数 lipid 種の混合 (POPC/POPE/CHOL 等) のサポート
- non-amino-acid permeant (small molecule + Auto-Martini) の埋め込み
- post-pull 時の `find_pbc_center_atom` を builder で自動計算 (現状は default で None、安全な小箱で動くが大きい bilayer では要手動指定)
