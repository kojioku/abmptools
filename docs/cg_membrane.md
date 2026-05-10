# abmptools.cg.membrane — Martini 3 Peptide-Membrane PMF Builder

Martini 3 粗視化 (CG) でペプチド-脂質膜系の PMF (Umbrella Sampling 経由) を
end-to-end 計算するビルダー。`insane` で bilayer 配置 → `cg.peptide` で
ペプチド CG 化 → topology compose → 13-31 window umbrella + WHAM までを
1 つの `MembraneCGBuilder.build()` 呼び出しで完結させる。

姉妹サブパッケージ:

- `abmptools.membrane` (AA、CHARMM36 / Lipid21 backend) — 全原子版の同機能
- `abmptools.cg.peptide` — Martini 3 ペプチド単独の CG 系 (本サブパッケージが
  内部で sub-call する)

CG だから AA の 30-100× 速い (KGG-POPC 31 windows × 5 ns で wall ~45 分、
8 thread CPU)。代償として left/right asymmetry など single-direction pulling
の sampling artifact は AA より顕著に出やすい。

## ライセンス上のルール (Martini 3 配布物の取り扱い)

このサブパッケージは **subprocess only** 戦略で、外部ツールのソースを
abmptools に同梱・改変しない設計:

- ✅ **`insane`** (PyPI、**GPL-2.0**) — Martini bilayer assembly。abmptools は
  `subprocess.run(["insane", ...])` で呼ぶだけ、ソース改変・同梱なし
  (GPL FAQ "mere aggregation" の範囲、abmptools 自体は Apache-2.0、v1.23.0+)
- ✅ **`vermouth-martinize`** (PyPI、Apache-2.0) — `cg.peptide` 経由で同様に
  subprocess only
- ✅ **`gmx`** (LGPL) — subprocess only
- ✅ **`tleap`** (LGPL、推奨) — `cg.peptide` 経由
- ⚠️ **Martini 3 force field `.itp`** (cgmartini.nl 配布) — **license 未明記**。
  本サブパッケージは未同梱。ユーザーが各自 `martini_v300.zip` を
  cgmartini-library.s3.ca-central-1.amazonaws.com から取得して `ff/` に
  4 ITP を unzip する設計。`validate` サブコマンドが取得手順を表示

このルールにより、abmptools 本体は **Apache-2.0 (v1.23.0+)、商用利用 OK** を
維持できる (subprocess 起動の GPL-2.0 ツールとの mere aggregation は GPL に
感染しない、FSF の公式見解)。

許可されているもの (GPL-2.0 互換 + Apache-2.0 互換):

- ✅ **insane** (GPL-2.0) — POPC / DOPC / DPPC 等の Martini bilayer 構築
- ✅ **vermouth-martinize** (Apache-2.0) — `martinize2` CLI 提供
- ✅ **AmberTools** `tleap` (LGPL) — atomistic peptide PDB (`cg.peptide` 経由)
- ✅ **GROMACS** (LGPL) — grompp / mdrun / wham / trjconv
- ✅ **Martini 3 force field 値の利用** — cgmartini が「academic + industrial
  community で広く使われる前提」で配布している前提

## 全体パイプライン (7 stage)

```
MembraneCGBuildConfig (JSON / YAML)
  ↓
MembraneCGBuilder.build():
  ↓ (0) _copy_ff_files     -> output_dir/{martini_v3.0.0,_solvents,_ions,_phospholipids}.itp (4 files)
  ↓ (1) _stage1_cg_peptide -> molecules/<name>/{name}_cg.pdb + {name}.itp
                              (cg.peptide sub-call: solvent_enabled=False, mdp_*=False)
  ↓ (2) _stage2_insane     -> bilayer.gro + insane_topol.top
  ↓ (3) _stage3_topology   -> topol.top (M3 4-ITP includes + peptide ITP +
                              Protein -> moleculetype name + NA+/CL- -> NA/CL)
  ↓ (4) _stage4_ions       -> system_ions.gro (insane の NaCl 配置を継承、
                              .gro 中の NA+/CL- atom name のみ正規化)
  ↓ (5) _stage5_index      -> index.ndx (Bilayer / Peptide / W / NA / CL / Non_Bilayer)
  ↓ (6) _stage6_mdps       -> mdp/{em,nvt,npt}.mdp + pull/pull.mdp +
                              windows/win_NNN/window.mdp (×N)
                              + bilayer pbc-atom 計算 (find_pbc_center_atom)
  ↓ (7) _stage7_run_script -> run.sh
```

`run.sh` の MD 実行順:

```
em -> nvt -> npt -> pull -> python -m abmptools.cg.membrane make-windows
   -> per-window grompp + mdrun (×N)
   -> python -m abmptools.cg.membrane wham
```

## API (主要 dataclass)

### `MembraneCGBuildConfig` (top-level)

```python
@dataclass
class MembraneCGBuildConfig:
    lipids: List[LipidMix]                     # v1.19 は 1 種のみ (POPC default)
    peptide: Optional[PeptideMembraneSpec]
    insane_d_nm: float = 8.0                   # insane の -d (xy box criterion)
    box_z_nm: float = 14.0                     # insane の -dz (z 寸法)
    insane_pbc: str = "hexagonal"
    insane_extra_args: List[str] = []          # extra argv to insane
    solvent_type: str = "W"                    # Martini W
    nacl_molar: float = 0.15
    neutralize: bool = True
    equilibration: EquilibrationCGProtocol
    pulling: PullingCGProtocol
    umbrella: UmbrellaCGProtocol
    martini_itp_dir: str                       # 4 ITP の場所 (空なら output_dir/ff)
    insane_path: str = "insane"
    martinize2_path: str = "martinize2"
    gmx_path: str = "gmx"
    tleap_path: str = "tleap"
    seed: Optional[int] = None
    output_dir: str = "."
    grompp_maxwarn: int = 5
```

### サブ dataclass

```python
@dataclass
class LipidMix:
    resname: str = "POPC"           # M3 phospholipid 名
    n_per_leaflet: int = 64         # 片葉あたりの脂質数

@dataclass
class PeptideMembraneSpec:
    name: str = "kgg"
    sequence: str = "KGG"           # 1-letter; or...
    cg_pdb_path: str = ""           # ...pre-built CG PDB
    cg_itp_path: str = ""           # ...と対の ITP
    initial_z_offset_nm: float = 3.0  # peptide-bilayer COM z 距離 (insane -dm)

@dataclass
class EquilibrationCGProtocol:
    em_steps: int = 5000
    nvt_nsteps: int = 25_000        # 0.5 ns at dt=20 fs
    npt_nsteps: int = 250_000       # 5 ns
    dt_ps: float = 0.020            # Martini 3 標準
    temperature_K: float = 310.0
    pressure_bar: float = 1.0
    tau_t_ps: float = 1.0
    tau_p_ps: float = 12.0          # Martini 3 標準

@dataclass
class PullingCGProtocol:
    pull_rate_nm_per_ps: float = -0.001  # peptide DOWN through bilayer
    pull_force_constant: float = 1000.0
    nsteps: int = 250_000           # 5 ns at dt=20 fs

@dataclass
class UmbrellaCGProtocol:
    z_min_nm: float = -1.5
    z_max_nm: float = +1.5
    window_spacing_nm: float = 0.25  # smoke; production は 0.10
    force_constant_kj_mol_nm2: float = 1000.0  # smoke; production は 500
    window_nsteps: int = 50_000     # smoke = 1 ns; production は 250000 = 5 ns
    pull_geometry: str = "direction"     # 動的 box 互換
    pull_vec: str = "0 0 1"
    pull_dim: str = "N N Y"
```

JSON / YAML I/O は `to_json`/`from_json` で対応。

## CLI

```bash
# example config 取得
python -m abmptools.cg.membrane example > kgg.json

# config + 依存性検証 (insane / gmx / martinize2 / tleap + 4 ITP)
python -m abmptools.cg.membrane validate --config kgg.json --ff-dir ./ff

# 7 stage build (MDP / topol.top / index.ndx / run.sh 生成)
python -m abmptools.cg.membrane build --config kgg.json --ff-dir ./ff -o ./out

# pull MD 後、各 window 用 start.gro を抽出 (run.sh が自動で呼ぶ)
python -m abmptools.cg.membrane make-windows \
    --pull-tpr pull/pull.tpr --pull-xtc pull/pull.xtc \
    --pull-xvg pull/pullx.xvg --windows-dir windows \
    --config input/config.json

# 全 window 完了後、WHAM PMF 計算 (run.sh が自動で呼ぶ)
python -m abmptools.cg.membrane wham \
    --windows-dir windows --analysis-dir analysis \
    --config input/config.json
```

## 設計判断とその根拠

### 1. `cg.peptide` を sub-call で再利用 (Stage 1)

CG ペプチド生成は `cg.peptide.PeptideCGBuilder` の責務。membrane builder は
Stage 1 で `solvent_enabled=False, mdp_*=False` で sub-call し、生成された
`<name>_cg.pdb` + `<name>.itp` だけを抽出して `molecules/<name>/` にコピー。

理由:

- DRY: martinize2 の subprocess 呼び出しロジックが 2 箇所に分散しない
- v1.18.0 の cg.peptide は M3 化が成熟しており再利用が無理なく済む
- sub-call は scratch dir (`output_dir/_cg_peptide_subbuild/`) で実行、
  cg.peptide が出す bookkeeping ファイル (`packed.gro`/`topol.top`/`run/run.sh`)
  が membrane builder の同名出力を上書きしないように隔離

### 2. `abmptools.membrane.{pulling,pmf,mdp_us_protocol}` の helpers を import 経由で再利用

CG / AA で本質的に同じ helper:

- `parse_pullx_xvg`, `_read_gro_positions`, `_read_ndx_groups`
- `find_pbc_center_atom`, `extract_window_frames`
- `estimate_initial_pull_coord`, `_gmx_trjconv_dump`, `_override_mdp_field`
- `render_pull_block` (MDP の [pull] / [pull-coord*] block render)
- `run_wham` (`gmx wham` のみで FF 非依存)

`UmbrellaCGProtocol` と AA の `USProtocol` は関連 field 名が一致しており、
duck-typing で AA 関数に CG config を渡せる。コード重複ゼロ。

### 3. MDP chassis: NVT-pull / NPT-window split

- **Pulling stage**: `pcoupl = no` (NVT) + `pull-coord1-geometry = direction-periodic`
- **Window stage**: `pcoupl = c-rescale, semiisotropic` (NPT) + `geometry = direction`

理由:

- `direction-periodic` は peptide が水→bilayer→反対側水を貫通する大変位
  (>0.49 × box) を扱える唯一の geometry
- しかし `direction-periodic` は **dynamic box 不可** (semiisotropic Pcoupl と
  非互換、grompp で fatal)
- そのため pull 段階は box 固定 (NVT)、window 段階は box 緩和を許す (NPT)
- pull は数 ns なので NVT でも bilayer はほぼ変形しない (実測 OK)
- これは AA membrane と同じ split

### 4. pbc-atom は **window MDP のみ** に注入

bilayer xy-extent は box の半分を超えるので、`pull-group1-pbcatom` を明示
しないと grompp が "Pull group 1 is larger than half the box" でエラー。

実装:

- Stage 6 で `pulling.find_pbc_center_atom` を `system_ions.gro` + `index.ndx`
  に対して呼んで bilayer COM 最近接の atom 番号を計算
- **window MDP のみ**にこの atom を `pull-group1-pbcatom` として注入
- **pull MDP には注入しない** — `direction-periodic` は内部で periodic image
  を扱うため pbcatom 不要、かつ CG (dt=20 fs) で両者を併用すると mdrun が
  step 0 で deadlock する経験的観察 (AA dt=2 fs では発生しない)

### 5. Insane の出力 topology を post-process

Insane は `topol.top` を出すが、abmptools 流儀に直接適合しない:

| 課題 | 解決 |
|---|---|
| `#include "martini.itp"` (generic) | M3 4-ITP の include 群に書き換え |
| `[ molecules ]` の `Protein` (placeholder) | peptide ITP の実 moleculetype 名 (`molecule_0`) で置換 |
| `NA+`/`CL-` (atom name + resname) | `NA`/`CL` に正規化 (gro + topol.top 両方) |

`topology_composer` モジュールが純 Python で text rewrite を実行
(re-pattern + column-aware substitution)。

### 6. Index は Python 直書き、`gmx make_ndx` を呼ばない

```
groups = {
    System,
    Bilayer    = (POPC, DOPC, ...)
    Peptide    = (LYS, GLY, ALA, ... 20 standard AA + variants)
    W          = (Martini W bead)
    NA, CL     = ions
    Non_Bilayer = System ∖ Bilayer  (2-group thermostat 用)
}
```

`system_packer.write_ndx_from_gro_cg` が `.gro` の residue 名から分類して
`.ndx` を直接書く。`gmx make_ndx` の stdin 駆動グループ追加 (脆い) を回避。

## 既知の caveat / 改善余地

| 問題 | 根本原因 | 対処 / 緩和策 |
|---|---|---|
| **PMF の左右非対称** (production で ~80 kJ/mol) | Single-direction pulling の sampling history bias | Bidirectional pulling + Bennett Acceptance、Replica-exchange umbrella sampling |
| **Lys 側鎖 dehydration cost が文献より高め** (~70-80 kJ/mol) | Martini 3 の charged residue parametrization の傾向 | 文献 (Vorobyov et al. 2014 ACS) で確認済の M3 系統的傾向、再パラメータ化が必要なら polyply / abmptools.cg.{polymer,smallmol} 経由で対応 |
| **smoke run の twin-peak + dip artifact** (z=-0.05) | k=1000 + Δz=0.25 での σ/spacing = 5σ → overlap 不足 | production protocol (k=500, Δz=0.10, 5 ns/window) に切替。`gmx wham` warning が 19→0 |
| **Multi-lipid 混合非対応** | v1.19 は単一 lipid 種限定 | v1.20+ で `LipidMix` の `ratio` field 追加予定 |
| **Small-molecule permeant 非対応** | `cg.smallmol/` (Auto-Martini) 未実装 | v1.20+ で予定 |

## デフォルト数値の理由

| パラメータ | 値 | 由来 |
|---|---|---|
| `dt_ps` | 0.020 | Martini 3 標準 (20 fs、AA の 10×) |
| `cutoff` | 1.1 nm | Martini 3 標準 (`coulombtype = reaction-field`、`epsilon_r = 15`) |
| `tau_p_ps` | 12.0 | Martini 3 標準 (AA の 5×) |
| `compressibility` | 3e-4 | Martini 3 標準 (水の 6×、CG の密度緩和を反映) |
| `tau_t_ps` | 1.0 | Martini 3 標準 (V-rescale) |
| `temperature_K` | 310.0 | 生体系の慣習 |
| `nacl_molar` | 0.15 | 生理食塩水濃度 |
| `pull_rate_nm_per_ps` | -0.001 | 1 nm/ns、AA membrane と同じ。負号は peptide DOWN sweep |
| `pull_force_constant` (pulling) | 1000.0 | AA membrane と同じ |
| `force_constant_kj_mol_nm2` (umbrella) | 1000 / 500 | smoke / production。500 で σ ≈ 0.07 nm、Δz=0.10 で overlap > 0.1 |
| `n_windows` (= z range / spacing + 1) | 13 / 31 | smoke / production |
| `window_nsteps` | 50,000 / 250,000 | smoke (1 ns) / production (5 ns) |

## ファイル構成

```
abmptools/cg/membrane/
├── __init__.py
├── __main__.py                 # `python -m abmptools.cg.membrane`
├── _subprocess.py              # cg.peptide からの re-export + logger
├── builder.py                  # MembraneCGBuilder.build() (7 stage)
├── cli.py                      # build/validate/example/make-windows/wham
├── forcefield_check.py         # 4 ITP + insane 含む REQUIRED 確認
├── insane_runner.py            # subprocess wrapper (GPL-2.0 respect)
├── mdp_templates.py            # CG MDP renderer (em/nvt/npt + pull + window)
├── models.py                   # 6 dataclass (config 5 段 nested)
├── pmf.py                      # AA membrane.pmf.run_wham に delegate
├── pulling.py                  # AA membrane.pulling re-export + CG pull MDP writer
├── system_packer.py            # add_ions_cg + write_ndx_from_gro_cg
├── topology_composer.py        # insane topology + .gro post-process
├── umbrella.py                 # write_window_mdps + write_run_script
└── README.md                   # ユーザー向け Quick Start
```

## 関連リソース

- `docs/tutorial_cg_membrane_us.md` — step-by-step チュートリアル
  (KGG-POPC で smoke 5-6 分 → production 45 分 → PMF 確認まで)
- `docs/membrane.md` — AA 版の同機能 (CHARMM36 / Lipid21)
- `sample/cg_membrane/kgg_popc_smoke.json` — 5-6 分で完走する smoke
- `sample/cg_membrane/kgg_popc_production.json` — 45 分で完走する production
- `docs/figures/cg_membrane_kgg_popc_pmf_smoke.png` — smoke PMF + histograms
- `docs/figures/cg_membrane_kgg_popc_pmf_production.png` — production PMF + smoke 比較
- `abmptools-sample` repo の `sample/membrane_us/peptide-kgg_POPC64_cg_us_*/` —
  軽量 archive (~35-60 MB、xtc 除外)
- OneDrive `abmptools-dump/membrane-us/peptide-kgg_POPC64_cg_us_*/` — 完全
  archive (~102-405 MB、xtc 含む) for re-analysis

## バージョン履歴

| version | 日付 | 主な変更 |
|---|---|---|
| 1.19.0 | 2026-05-05 | サブパッケージ新設。Initial release with smoke + production validation |
