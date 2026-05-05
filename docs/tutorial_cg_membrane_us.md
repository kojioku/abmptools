# Tutorial: Martini 3 Peptide-Membrane PMF (CG Umbrella Sampling)

ペプチドの脂質膜透過 PMF を `abmptools.cg.membrane` (Martini 3 + insane) で
計算する step-by-step チュートリアル。`abmptools.membrane` (AA、CHARMM36 /
Lipid21) の **30-100× 速い CG 版**。

**対象**: GROMACS と AmberTools が動く Linux / WSL2 環境を持っているユーザー。
Martini 3 の基本概念 (W bead = 4 H₂O、cutoff 1.1 nm reaction-field、
dt = 20 fs) は既知とする。

**所要時間** (KGG ペプチド + POPC 64 lipid/leaflet):

| ステップ | wall-time |
|---|---|
| 1. 環境準備 | 10-30 分 (初回のみ) |
| 2. 系の設定 + build | 1-2 分 |
| 3. **Smoke** MD (CPU 8-thread) | **~5-6 分** |
| 4. Smoke PMF 確認 | 1 分 |
| 5. **Production** MD (CPU 8-thread) | **~45 分** |
| 6. Production PMF 確認 + smoke 比較 | 1 分 |

最初は §1〜§4 を smoke で動作確認、§5〜§6 で production 品質に拡張する流れ。

---

## §0. 全体像

```
[Step 1] 環境準備 (一度きり)
    ├─ abmptoolsenv (vermouth + insane + gmx + tleap)
    └─ Martini 3 force field 4 ITP を ff/ に展開
[Step 2] 系の設定 (JSON)
    ├─ MembraneCGBuildConfig: peptide / lipid / window 仕様
    └─ MembraneCGBuilder.build() で全 input 生成 (~1 分)
[Step 3] MD 実行 (smoke / production を切替)
    └─ NT=8 bash run.sh   (~5 分 smoke / ~45 分 production)
        ├─ em → nvt → npt → pull
        ├─ make-windows (per-window start.gro 抽出)
        ├─ per-window mdrun (×13 smoke or ×31 production)
        └─ gmx wham → analysis/pmf.xvg
[Step 4] 結果確認
    ├─ analysis/pmf.xvg               PMF(z) 数値
    ├─ analysis/histo.xvg             per-window histogram (overlap チェック)
    ├─ analysis/gmx_wham.log          "Check your histograms" warning 数
    └─ docs/figures/cg_membrane_*.png 例の plot
[Step 5] 反復 / 本番化
    ├─ smoke 13 win × 1 ns × k=1000 → twin-peak + dip artifact
    └─ production 31 win × 5 ns × k=500 → smooth single-peak PMF
```

---

## §1. 環境準備

### 1.1 abmptoolsenv セットアップ (まだなら)

```bash
# micromamba 環境を作成 (一度きり)
micromamba create -n abmptoolsenv python=3.10 -c conda-forge
micromamba activate abmptoolsenv

# core 依存
pip install -e ~/path/to/abmptools[cg]    # vermouth + insane + pyyaml

# 外部ツール (subprocess で呼ぶ)
micromamba install -c conda-forge gromacs ambertools

# 確認
which insane martinize2 gmx tleap
# /home/.../envs/abmptoolsenv/bin/insane
# /home/.../envs/abmptoolsenv/bin/martinize2
# /home/.../envs/abmptoolsenv/bin.AVX2_256/gmx
# /home/.../envs/abmptoolsenv/bin/tleap
```

### 1.2 Martini 3 force field 4 ITP の取得

```bash
mkdir ff && cd ff
curl -L -o /tmp/m300.zip \
    https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip
unzip -o -j /tmp/m300.zip \
    "martini_v300/martini_v3.0.0.itp" \
    "martini_v300/martini_v3.0.0_solvents_v1.itp" \
    "martini_v300/martini_v3.0.0_ions_v1.itp" \
    "martini_v300/martini_v3.0.0_phospholipids_v1.itp" \
    -d .
cd ..
ls ff/
# martini_v3.0.0.itp
# martini_v3.0.0_solvents_v1.itp
# martini_v3.0.0_ions_v1.itp
# martini_v3.0.0_phospholipids_v1.itp
```

**Note**: cgmartini.nl 配布物は license 未明記なので abmptools には未同梱。
`.gitignore` に `ff/` を入れて誤 commit を防ぐ。

### 1.3 依存性検証

```bash
python -m abmptools.cg.membrane example > /tmp/probe.json
python -m abmptools.cg.membrane validate --config /tmp/probe.json --ff-dir ./ff
```

期待出力:

```
Configuration: OK (POPC bilayer, 128 lipids total, peptide=kgg, T=310.0 K, NaCl 0.15 M, 13 umbrella windows)

External tools:
  [OK]    insane      (/home/.../envs/abmptoolsenv/bin/insane)
  [OK]    martinize2  (/home/.../envs/abmptoolsenv/bin/martinize2)
  [OK]    gmx         (/home/.../envs/abmptoolsenv/bin.AVX2_256/gmx)
  [OK]    tleap       (/home/.../envs/abmptoolsenv/bin/tleap)

Martini 3 force field files (in ./ff):
  [OK]      martini_v3.0.0.itp
  [OK]      martini_v3.0.0_solvents_v1.itp
  [OK]      martini_v3.0.0_ions_v1.itp
  [OK]      martini_v3.0.0_phospholipids_v1.itp
```

全 `[OK]` でないと先に進めない。`[MISS]` ならツール不足、`[MISSING]` なら ITP 取得不足。

---

## §2. Smoke run (5-6 分で完走、pipeline 検証用)

### 2.1 Sample config を取得

```bash
cp /path/to/abmptools/sample/cg_membrane/kgg_popc_smoke.json kgg.json
```

または `python -m abmptools.cg.membrane example > kgg.json`。

主な設定 (smoke):

```json
{
  "peptide": {"sequence": "KGG", "initial_z_offset_nm": 3.0},
  "lipids": [{"resname": "POPC", "n_per_leaflet": 64}],
  "umbrella": {
    "window_spacing_nm": 0.25,        // smoke
    "force_constant_kj_mol_nm2": 1000, // smoke
    "window_nsteps": 50000             // 1 ns/window
  }
}
```

### 2.2 build (~1 分)

```bash
python -m abmptools.cg.membrane build --config kgg.json --ff-dir ./ff -o ./smoke_out
```

期待出力 (末尾):

```
=== Stage 7/7: run.sh ===
Build complete. Output: ./smoke_out
  final coordinates: ./smoke_out/system_ions.gro
  topology:          ./smoke_out/topol.top
  index:             ./smoke_out/index.ndx
  run script:        ./smoke_out/run.sh
  windows: 13, lipid: POPC
```

ディレクトリ構造:

```
smoke_out/
├── input/config.json           # config dump
├── molecules/kgg/{kgg_atomistic.pdb, kgg_cg.pdb, kgg.itp}
├── martini_v3.0.0*.itp (4 files、ff_dir からコピー)
├── topol.top
├── system_ions.gro
├── index.ndx
├── bilayer.gro                 # raw insane 出力 (audit)
├── insane_topol.top            # raw insane topology (audit)
├── mdp/{em,nvt,npt}.mdp
├── pull/pull.mdp
├── windows/win_000/window.mdp .. windows/win_012/window.mdp
└── run.sh
```

### 2.3 MD 実行 (~5-6 分 smoke、CPU 8-thread)

```bash
cd smoke_out
NT=8 bash run.sh
```

`NT=8` は OpenMP thread 数。CPU core 数に合わせて調整。`run.sh` は次の順:

1. em (~ < 1 秒)
2. nvt (0.5 ns ~ 8 秒)
3. npt (5 ns ~ 1.5 分)
4. pull (5 ns ~ 1.5 分)
5. make-windows (frame 抽出 ~ 数秒)
6. 13 windows × 1 ns ≈ 13 × 15 秒 = 3 分
7. wham (数秒)

**合計 ~5-6 分** (Martini 3 で ~5300 ns/day @ 8 thread)。

### 2.4 結果確認

```bash
ls analysis/
# gmx_wham.log  histo.xvg  pmf.xvg  pullx.dat  tpr.dat

# Check your histograms warning 数 (smoke は多い)
grep -c "Check your histograms" analysis/gmx_wham.log
# 19  ← smoke は overlap 不足のため大量 warning

# PMF 数値
head -5 analysis/pmf.xvg
```

PMF を簡易 plot:

```bash
python << 'PY'
import matplotlib.pyplot as plt
data=[]
for line in open('analysis/pmf.xvg'):
    if line.startswith(('@','#')): continue
    p=line.split()
    if len(p)>=2: data.append([float(p[0]), float(p[1])])
import numpy as np
arr = np.array(data)
plt.plot(arr[:,0], arr[:,1])
plt.xlabel('z (nm)'); plt.ylabel('PMF (kJ/mol)')
plt.title('Smoke PMF (KGG-POPC)')
plt.savefig('smoke_pmf.png')
PY
```

### 2.5 Smoke は何が見えれば成功か

- ✅ 13/13 windows 完走 (`ls windows/win_*/window.gro | wc -l` → 13)
- ✅ pull traj 範囲 +3 → -2 nm (peptide が bilayer を貫通)
- ✅ `analysis/pmf.xvg` 生成、200 点
- ⚠️ gmx wham warning 19 件 (default protocol の宿命、production で解消)
- ⚠️ PMF が twin-peak with z=-0.05 dip (sampling artifact、production で解消)

## §3. Smoke の限界 (sampling artifact の正体)

smoke の PMF が "見栄え悪い" 理由は **histogram overlap 不足**:

```
# Each window の z 標準偏差 (k=1000 + 310K で σ ≈ √(kT/k))
σ ≈ √(8.314e-3 × 310 / 1000) ≈ 0.05 nm

# Spacing 0.25 nm vs σ 0.05 nm
separation = 0.25 / 0.05 = 5σ

# 隣接 window の Gauss が ~5σ 離れている → overlap < 0.0001
# WHAM が PMF を重み付き連結できず、warning 連発
```

これは **Martini 3 固有ではなく** Umbrella Sampling default protocol の限界。
AA membrane (CHARMM Phase D) でも同じ症状を観察済み (memory
`project_membrane_charmm36_gotchas.md`)。

解決: **k を下げて σ を広げ、spacing を狭めて overlap を確保** = §4 production。

---

## §4. Production run (45 分、研究品質)

### 4.1 Production config

```bash
cp /path/to/abmptools/sample/cg_membrane/kgg_popc_production.json kgg_prod.json
```

主な変更点 (smoke 比):

| パラメータ | smoke | production | 効果 |
|---|---|---|---|
| `window_spacing_nm` | 0.25 | **0.10** | windows ×2.4 (13→31) |
| `force_constant_kj_mol_nm2` | 1000 | **500** | σ ×√2 ≈ 0.07 nm |
| `window_nsteps` | 50000 (1 ns) | **250000 (5 ns)** | per-window sampling ×5 |

理論的 overlap:

```
σ ≈ √(8.314e-3 × 310 / 500) ≈ 0.072 nm
separation = 0.10 / 0.072 = 1.32σ → overlap > 0.1 ✅
```

### 4.2 build + run

```bash
python -m abmptools.cg.membrane build \
    --config kgg_prod.json --ff-dir ./ff -o ./prod_out
cd prod_out
NT=8 bash run.sh
```

Wall time 内訳 (8-thread CPU):

| stage | time |
|---|---|
| em + nvt + npt + pull | ~3 分 |
| make-windows | ~10 秒 |
| 31 windows × 5 ns | ~42 分 |
| wham | ~10 秒 |
| **total** | **~45 分** |

途中経過の確認 (別 terminal):

```bash
grep -c "Performance:" run.log    # 完了 mdrun の数 (em+nvt+npt+pull で 4、+ windows 数)
ls windows/win_*/window.gro | wc -l   # 完走 windows 数 (max 31)
```

### 4.3 Production の結果確認

```bash
# warning が劇的に減るのが定石
grep -c "Check your histograms" analysis/gmx_wham.log
# 0   ← 完璧
```

per-window 標準偏差:

```bash
python << 'PY'
from pathlib import Path; import statistics
sigmas = []
for d in sorted(Path("windows").glob("win_*")):
    zs = []
    for line in open(d / "pullx.xvg"):
        if line.startswith(('@','#')): continue
        p = line.split()
        if len(p)>=2:
            try: zs.append(float(p[1]))
            except: pass
    if zs and len(zs)>1: sigmas.append(statistics.stdev(zs))
print(f"σ mean: {statistics.mean(sigmas):.3f} nm")
print(f"separation: {0.10/statistics.mean(sigmas):.2f}σ")
PY
# σ mean: 0.076 nm
# separation: 1.32σ
```

### 4.4 Smoke vs production 比較 plot

```bash
python << 'PY'
import matplotlib.pyplot as plt
import numpy as np

def parse(path):
    arr = []
    for line in open(path):
        if line.startswith(('@','#')): continue
        p = line.split()
        if len(p)>=2:
            try: arr.append([float(p[0]), float(p[1])])
            except: pass
    return np.array(arr)

smoke = parse('../smoke_out/analysis/pmf.xvg')
prod  = parse('analysis/pmf.xvg')

plt.figure(figsize=(8, 5))
plt.plot(smoke[:,0], smoke[:,1], 'r-', alpha=0.6, label='Smoke (13win × 1ns × k=1000 × Δz=0.25)')
plt.plot(prod[:,0], prod[:,1], 'b-', lw=2, label='Production (31win × 5ns × k=500 × Δz=0.10)')
plt.xlabel('z (nm)'); plt.ylabel('PMF (kJ/mol)')
plt.title('KGG (+1) PMF through POPC bilayer (Martini 3)')
plt.axhline(0, color='k', lw=0.5)
plt.axvline(0, color='gray', ls='--', alpha=0.5)
plt.legend(); plt.grid(alpha=0.3)
plt.savefig('compare_pmf.png', dpi=120)
PY
```

参考: 同じ plot は `abmptools/docs/figures/cg_membrane_kgg_popc_pmf_production.png`
にある (smoke vs production 比較 + per-window histograms)。

### 4.5 Production の代表値 (KGG +1 in POPC)

| 量 | 値 |
|---|---|
| Bilayer-center barrier peak | +68 kJ/mol at z = -0.49 nm |
| 水中 minimum (z = +1.7 nm 側) | -82 kJ/mol |
| 水中 minimum (z = -1.7 nm 側、anchor) | 0 kJ/mol |
| 左右非対称 (single-direction pull の宿命) | 82 kJ/mol |
| Barrier height (vs water、平均) | ~70-90 kJ/mol |

文献的に Lys (+1) の dehydration cost は ~30-80 kJ/mol (Vorobyov et al.
2014 ACS Chem Biol 等)。我々の ~70-90 kJ/mol は M3 の上端側で、対 +1 charge の
M3 parametrization は若干高めに出る系統傾向と一致。

---

## §5. Single-direction pull の左右非対称をどう扱うか

物理的には、bilayer は homogeneous なので PMF は **z=0 を中心に対称**であるべき。
しかし single-direction pull (peptide を +z → -z に 1 回 sweep) では、各 window
の sampling history が偏るため、左右で 30-80 kJ/mol の系統的 offset が出る。

対処法 (実装難度順):

1. **何もしない** (smoke / 早い検証): asymmetry は inherent と認識、PMF 形状の
   定性確認のみに使う
2. **両端を強制 align** (post-process): scipy で fit して midplane で symmetrize
3. **Bidirectional pulling + Bennett Acceptance** (ガチ): pull を +→- と -→+ の
   両方向で実行、各 window の forward/reverse force を Bennett 式で結合。
   `gmx wham` には bootstrap オプションあるが本格対応は別実装が必要
4. **Replica-exchange umbrella sampling**: 隣接 window 間で交換、kinetic
   trapping を緩和。[plumed](https://www.plumed.org/) 経由が定番

abmptools.cg.membrane 1.19.0 は (1) のみサポート。(2)〜(4) は v1.20+ または
ユーザーが post-process / 別ツールで実施。

---

## §6. 大きい系・別 lipid で試す

### 6.1 別 lipid (DOPC, POPE 等)

```json
{"lipids": [{"resname": "DOPC", "n_per_leaflet": 64}]}
```

Martini 3 phospholipids ITP (`martini_v3.0.0_phospholipids_v1.itp`) には
POPC / DOPC / POPE / DOPE / POPG / DPPC 等が含まれている。`gmx grompp` が
moleculetype を見つけられればそのまま動く。

### 6.2 大きい系 (lipid を増やす)

```json
{"lipids": [{"resname": "POPC", "n_per_leaflet": 128}]}
```

bilayer xy-extent が増えるので box サイズ自動調整 (insane が `-d` を満たすよう
xy を広げる)。wall time は ~2× (atom 数 ~2×、CG perf ~半減)。

### 6.3 大きいペプチド

```json
{"peptide": {"sequence": "RRRR", "initial_z_offset_nm": 3.0}}
```

Lys (+1) → Arg (+4 charge) で barrier が更に大きくなる。M3 charged residue
parametrization は Arg/Lys/Asp/Glu すべて支持。

複数残基や capping (ACE/NME) は cg.peptide が tleap で生成、CG 化は martinize2
が処理。`PeptideMembraneSpec.sequence` に 1-letter で書くだけで OK。

---

## §7. 失敗パターンと対処

### 7.1 grompp が "Pull group 1 is larger than half the box" で fatal

→ pbc-atom が window MDP に注入されていない。abmptools.cg.membrane 1.19.0 は
   builder で自動注入する設計。手動 build なら `find_pbc_center_atom` を呼んで
   `pull-group1-pbcatom` を render_pull_block の引数に渡す。

### 7.2 pull mdrun が step 0 で hang (99% CPU 無限)

→ pull MDP に pbcatom が指定されている。`direction-periodic` geometry は
   periodic image を内部処理するので pbcatom 不要。CG (dt=20 fs) で両者を
   併用すると integrator deadlock を起こす (AA dt=2 fs では発生しない)。
   abmptools.cg.membrane 1.19.0 は **windows のみ** に注入する設計で回避済み。

### 7.3 全 window の start.gro が同じ frame

→ `pull_rate_nm_per_ps` の符号が逆。`PeptideMembraneSpec.initial_z_offset_nm = +3.0`
   (peptide above bilayer) なら **rate は負** (peptide DOWN sweep)。
   1.19.0 のデフォルトは `-0.001`。

### 7.4 system_ions.gro の atom name で grompp が大量 warning

→ insane の `NA+` / `CL-` atom name と Martini 3 ions ITP の `NA` / `CL` の
   不一致。`topology_composer.normalize_ion_atom_names_gro` が builder の
   Stage 4 で自動処理する。`-maxwarn 5` でも warning が残るなら、builder の
   `_stage4_ions` が呼ばれていない可能性。

### 7.5 windows 一部失敗で wham が "no window.tpr / pullx.xvg pairs"

→ `set -euo pipefail` で run.sh が中断したか、特定 window で blow up。
   `windows/win_NNN/window.log` を tail で確認、`Step ... nan potential energy`
   なら start.gro の問題 (window mdp の pbcatom or pull init を確認)。

---

## §8. 参考リンク

- `docs/cg_membrane.md` — 設計詳解 + API + design rationale
- `abmptools/cg/membrane/README.md` — Quick Start
- `docs/membrane.md` — AA 版 (CHARMM36 / Lipid21)
- `docs/tutorial_membrane_us.md` — AA 版 tutorial (本 tutorial の sister)
- `sample/cg_membrane/{kgg_popc_smoke.json, kgg_popc_production.json}` — config
- `docs/figures/cg_membrane_kgg_popc_pmf_*.png` — example PMF plots
- `abmptools-sample` repo `sample/membrane_us/peptide-kgg_POPC64_cg_us_*/` —
  軽量 archive (~35-60 MB、xtc 除外)
- OneDrive `abmptools-dump/membrane-us/peptide-kgg_POPC64_cg_us_*/` — 完全
  archive (~102-405 MB、xtc 含む) for re-analysis
- Insane (GPL-2.0): https://github.com/Tsjerk/Insane
- Martini 3: https://cgmartini.nl/
- vermouth-martinize (Apache-2.0): https://github.com/marrink-lab/vermouth-martinize
- Souza et al. 2021, *Nat. Methods* 18:382-388 (Martini 3 reference paper)
