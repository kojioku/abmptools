# Tutorial: Peptide-Bilayer Umbrella Sampling (PMF) Pipeline

ペプチドの脂質膜透過 PMF を `abmptools.membrane` で計算するための step-by-step
チュートリアル。全コマンド実行例と実測データ付き。

**対象**: GROMACS と AmberTools が動く Linux / WSL2 環境を持っているユーザー。
基本的な MD 概念 (NPT、PME、constraints 等) は既知とする。

**所要時間** (poly-Ala 5-mer + POPC 64 lipid の最小例):

| ステップ | wall-time |
|---|---|
| 1. 環境準備 | 10-30 分 (初回のみ) |
| 2. 系の設定 + build | 1-2 分 |
| 3. MD (CPU 4-core) | ~7 時間 |
| 3'. MD (GPU RTX 4070 Ti) | ~30 分 (smoke) / ~3 時間 (本番) |
| 4. 結果確認 | 数分 |

最初は §1〜§5 を一通り通して動作確認、その後 §6 で本番品質に拡張する流れを推奨。

---

## §0. 全体像

```
[Step 1] 環境準備
    ├─ abmptoolsenv (tleap / packmol-memgen / parmed / abmptools)
    ├─ gmxcudaenv  (任意、GPU 利用時)
    └─ NumPy パッチ (packmol-memgen 互換)
[Step 2] 系の設定 (Python script)
    ├─ MembraneConfig: peptide / lipid / window 仕様
    └─ MembraneUSBuilder.build() で全入力ファイル生成
[Step 3] MD 実行
    └─ bash run.sh (CPU / GPU 切替)
[Step 4] 結果確認
    ├─ analysis/pmf.xvg          PMF(z)
    ├─ analysis/histo.xvg        per-window histogram (収束チェック)
    └─ pull/pullx.xvg            反応座標 trace
[Step 5] 反復 / 本番化
    └─ window 数 / 時間 / pull rate を調整
```

詳細リファレンスは [`membrane.md`](membrane.md) を参照。本チュートリアルは
**実行手順**に絞る。

---

## §1. 環境準備

### 1.1 abmptoolsenv (tleap + packmol-memgen + parmed + abmptools)

```bash
# 既存 abmptoolsenv があればこの手順は不要。
micromamba create -n abmptoolsenv -c conda-forge \
    python=3.10 \
    gromacs=2021 \
    'ambertools>=22' \
    parmed \
    openff-toolkit \
    numpy

# editable install of abmptools
cd ~/llm-project/fcews-workspace/abmptools
~/.local/share/mamba/envs/abmptoolsenv/bin/python -m pip install -e .
```

### 1.2 NumPy ≥ 1.24 互換性パッチ (packmol-memgen)

`packmol-memgen` 同梱の `pdbremix/v3numpy.py` は削除済 `np.float` を使う:

```bash
ENV=~/.local/share/mamba/envs/abmptoolsenv
LIB=$ENV/lib/python3.10/site-packages/packmol_memgen/lib/pdbremix
sed -i 's/np\.zeros(3, dtype=np\.float)/np.zeros(3, dtype=float)/' $LIB/v3numpy.py
sed -i 's/np\.array(args, dtype=np\.float, copy=True)/np.array(args, dtype=float, copy=True)/' $LIB/v3numpy.py
```

冪等。再実行しても問題なし。

### 1.3 (任意) GPU 加速用の gmxcudaenv

bioconda 版 `gromacs=2021.3` は OpenCL ビルドで、WSL2 では NVIDIA GPU が
見えない。CUDA ビルドの GROMACS だけを別 env に install:

```bash
micromamba create -n gmxcudaenv -c conda-forge 'gromacs[build=nompi_cuda*]' -y
~/.local/share/mamba/envs/gmxcudaenv/bin/gmx --version | grep "GPU support"
# → "GPU support: CUDA" が出れば OK
```

### 1.4 動作確認 (smoke shell test)

```bash
cd ~/llm-project/fcews-workspace/abmptools
bash tests/integration/run_membrane_us_smoke.sh
```

`=== smoke test PASSED ===` まで通れば環境 OK (~30 秒)。

---

## §2. 系の設定 + build

### 2.1 設定スクリプト

`run_my_pmf.py` のような Python スクリプトを書く:

```python
# run_my_pmf.py
import logging
logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")

from abmptools.membrane import (
    MembraneConfig, MembraneUSBuilder,
    LipidSpec, PeptideSpec, IonSpec,
    USProtocol, EquilibrationProtocol, PullingProtocol,
)

# --- US 設定: 13 windows, 5 ns each ---
us = USProtocol(
    z_min_nm=-1.5, z_max_nm=+1.5,
    window_spacing_nm=0.25,                  # → 13 windows
    window_nsteps=2_500_000, window_dt_ps=0.002,  # 5 ns/window
    window_nstxtcout=1000, window_nstenergy=500,
    force_constant_kj_mol_nm2=1000.0,
)

# --- 平衡化: 0.2 ns NVT + 5 ns NPT ---
eq = EquilibrationProtocol(
    em_steps=50_000, em_tol=1000.0,
    nvt_nsteps=100_000,                       # 0.2 ns
    npt_nsteps=2_500_000,                     # 5 ns
    dt_ps=0.002, temperature_K=310.0,
)

# --- 反応座標生成 (NVT pull) ---
pull = PullingProtocol(
    pull_rate_nm_per_ps=0.001,                # 10 ns で 10 nm 変位
    pull_force_constant=1000.0,
    nsteps=5_000_000,                         # 10 ns
    nstxout_compressed=500,
)

cfg = MembraneConfig(
    backend="amber",                          # or "charmm36"
    lipids=[LipidSpec(resname="POPC", n_per_leaflet=32)],
    peptide=PeptideSpec(
        name="aa5", sequence="AAAAA",
        cap_n="ACE", cap_c="NME",
    ),
    ions=IonSpec(salt_concentration_M=0.15),
    output_dir="./run01",
    seed=42,
    box_xy_nm=6.0,
    water_thickness_nm=2.5,
    distance_to_lipid_nm=1.5,
    equilibration=eq, pulling=pull, umbrella=us,
)

# --- 混合脂質の例 (v1.17.1+) ---
# POPC + CHL1 4:1 (lipid raft 模型) の場合は LipidSpec を並べるだけ:
#
# cfg = MembraneConfig(
#     backend="amber",
#     lipids=[
#         LipidSpec(resname="POPC", n_per_leaflet=80),
#         LipidSpec(resname="CHL1", n_per_leaflet=20),
#     ],
#     peptide=PeptideSpec(name="aa5", sequence="AAAAA"),
#     ...
# )
#
# 内部処理:
#   --lipids POPC:CHL1 --ratio 4:1   (n_per_leaflet を gcd 約分)
#   --distxy_fix sqrt(80*67 + 20*38) ≈ 78.2  (per-lipid APL)
# APL は LipidSpec.apl_angstrom2 を省略すると DEFAULT_LIPID_APL から
# 自動取得 (POPC=67, CHL1=38, POPE=56, DOPC=72, DPPC=63, ... 詳細は
# membrane.md の表を参照)。低温 gel 相など特殊な条件では explicit に:
#   LipidSpec(resname="DPPC", n_per_leaflet=64, apl_angstrom2=49.0)

result = MembraneUSBuilder(cfg).build()
print(f"\nbuild OK")
print(f"  output_dir : {result['output_dir']}")
print(f"  run.sh     : {result['run_script']}")
print(f"  windows    : {len(result['window_mdps'])}")
```

### 2.2 build 実行

```bash
ENV=~/.local/share/mamba/envs/abmptoolsenv
PATH=$ENV/bin:$PATH AMBERHOME=$ENV $ENV/bin/python3 run_my_pmf.py
```

実行例 (~1 分):

```
abmptools.membrane.builder | === Stage 1: bilayer construction ===
abmptools.membrane.bilayer | running tleap to build peptide
abmptools.membrane.bilayer | running packmol-memgen
abmptools.membrane.bilayer | bilayer PDB written: run01/build/bilayer_peptide.pdb
abmptools.membrane.builder | === Stage 2: force-field parameterisation (amber) ===
abmptools.membrane.parameterize_amber | running tleap
abmptools.membrane.parameterize_amber | parmed: saving run01/build/system.top (gromacs)
abmptools.membrane.builder | === Stage 3: equilibration MDPs ===
abmptools.membrane.builder | === Stage 4: pulling MDP ===
abmptools.membrane.builder | estimated initial pull-coord: -4.656 nm
abmptools.membrane.builder | pbc-atom for Bilayer: 6826
abmptools.membrane.builder | === Stage 5: umbrella window MDPs ===
abmptools.membrane.umbrella | wrote 13 window MDPs in run01/windows
build OK
  run.sh     : run01/run.sh
  windows    : 13
```

> ⚠ もし下記の警告が出たら **§6.1 を読む** (pull rate 不足):
> ```
> WARNING: Pull range may not cover all windows!
>   initial peptide z (rel. bilayer): -4.656 nm
>   z_max (last window):              +1.500 nm
>   required displacement:             6.156 nm
>   configured displacement:           5.000 nm
>   ...
> ```

### 2.3 ディレクトリレイアウト確認

```bash
tree run01 -L 2 -d
# run01/
# ├── analysis/
# ├── build/        ← system.top / system.gro / system.ndx 等
# ├── equil/        ← em.mdp / nvt.mdp / npt.mdp
# ├── input/        ← config.json (再現性)
# ├── pull/         ← pull.mdp
# └── windows/
#     ├── win_000/  ← window.mdp (rate=0, init=-1.5)
#     ├── win_001/
#     ├── ...
#     └── win_012/
```

---

## §3. MD 実行

### 3.1 CPU 実行 (4 cores、~7 時間)

```bash
cd run01
ENV=~/.local/share/mamba/envs/abmptoolsenv
PATH=$ENV/bin:$PATH AMBERHOME=$ENV \
    NT=4 \
    bash run.sh > run.log 2>&1 &
disown
```

`disown` で sshセッションを切っても継続。`tail -f run.log` でログ追跡。

### 3.2 GPU 実行 (RTX 4070 Ti、~3 時間)

```bash
cd run01
ENV=~/.local/share/mamba/envs/abmptoolsenv
CUDAENV=~/.local/share/mamba/envs/gmxcudaenv
PATH=$ENV/bin:$PATH AMBERHOME=$ENV \
    GMX=$CUDAENV/bin/gmx \
    MDRUN_OPTS="-nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu -pin on" \
    NT=4 \
    bash run.sh > run.log 2>&1 &
disown
```

ポイント:
- `GMX=$CUDAENV/bin/gmx` で **CUDA build** を指す
- `MDRUN_OPTS` で各 mdrun に GPU offload オプションを追加
- `PATH` は abmptoolsenv 優先 (python / tleap / packmol-memgen 用)
- run.sh 内で grompp と mdrun が同じ `$GMX` を使うので .tpr version mismatch しない

### 3.3 進行状況の確認

```bash
# 現在のステージ
grep -E "Stage|complete" run01/run.log | tail -5

# 各 stage の所要時間 (mtime)
stat -c "%y %n" run01/equil/*.gro run01/pull/*.gro run01/windows/win_*/window.gro 2>/dev/null

# GPU 使用率
nvidia-smi

# 実行中の gmx プロセス
pgrep -af gmx
```

---

## §4. 結果確認

### 4.1 PMF (analysis/pmf.xvg)

```bash
# z 範囲 / E 範囲
grep -v "^[#@]" run01/analysis/pmf.xvg | awk '
    BEGIN{mz=99;Mz=-99;me=99;Me=-99}
    {if($1<mz)mz=$1; if($1>Mz)Mz=$1; if($2<me)me=$2; if($2>Me)Me=$2}
    END{printf "z: %.3f to %.3f nm\nE: %.2f to %.2f kJ/mol\n", mz,Mz,me,Me}'

# プロット (xmgrace)
xmgrace run01/analysis/pmf.xvg

# プロット (matplotlib)
~/.local/share/mamba/envs/abmptoolsenv/bin/python3 - <<'PY'
import numpy as np, matplotlib.pyplot as plt
data = np.loadtxt("run01/analysis/pmf.xvg", comments=["#", "@"])
plt.plot(data[:,0], data[:,1])
plt.xlabel("z (nm)"); plt.ylabel("PMF (kJ/mol)"); plt.grid(True)
plt.savefig("run01/analysis/pmf.png", dpi=150)
print("saved pmf.png")
PY
```

### 4.2 Histogram overlap (収束チェック)

`analysis/histo.xvg` を可視化。各 window の histogram が **隣接 windowと
30% 以上 overlap** していれば WHAM の信頼性が高い:

```bash
xmgrace run01/analysis/histo.xvg
```

overlap が悪い領域 (gap がある) では window を細かくする (§6.3)。

### 4.3 反応座標 trace (pull/pullx.xvg)

pull stage で peptide が想定通りに traverse したか確認:

```bash
~/.local/share/mamba/envs/abmptoolsenv/bin/python3 - <<'PY'
import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt("run01/pull/pullx.xvg", comments=["#", "@"])
plt.plot(d[:,0]/1000, d[:,1])
plt.xlabel("time (ns)"); plt.ylabel("pull-coord z (nm)")
plt.axhline(-1.5, color="g", ls=":", label="z_min")
plt.axhline(+1.5, color="r", ls=":", label="z_max")
plt.legend(); plt.grid(True)
plt.savefig("run01/pull/pullx_trace.png", dpi=150)
print("saved pullx_trace.png")
PY
```

`z_min..z_max` の範囲を**完全にカバー**していなければ §6.1 で pull
パラメータを再調整。

### 4.4 期待される PMF の形

ペプチド + 脂質膜の場合:
- **z = ±2 nm 付近 (水相)**: PMF ≈ 0 (基準)
- **z = ±1.5 nm (頭部基)**: 浅い極小 (溶媒和シェル)
- **z = ±0.5 nm (尾部)**: 立ち上がり
- **z = 0 nm (中心)**: 最大 (free-energy barrier)
- **対称**: PMF(z) ≈ PMF(-z) で対称になることを確認 (非対称なら sampling
  不足、§6 参照)

---

## §5. 反復・本番化

### 5.1 Smoke から本番への拡張パス

| 段階 | windows | window 時間 | pull 時間 | 合計 MD | GPU walltime |
|---|---|---|---|---|---|
| tiny smoke | 7 (Δ=0.5) | 50 ps | 200 ps | ~1 ns | ~3 分 |
| production smoke | 7 (Δ=0.5) | 1 ns | 5 ns | ~13 ns | ~30 分 |
| **本番より** (このチュートリアル) | 13 (Δ=0.25) | 5 ns | 10 ns | ~81 ns | ~3 時間 |
| production-quality | 30 (Δ=0.1) | 50 ns | 30 ns | ~1500 ns | ~3 日 (1 GPU) |

### 5.2 収束チェック

PMF が安定したかを確認する基本手順:

```bash
# (a) 時間ブロック分割: 各 window を前半 / 後半 に分けて WHAM 比較
# 前半:
~/.local/share/mamba/envs/gmxcudaenv/bin/gmx wham \
    -it run01/analysis/tpr.dat -ix run01/analysis/pullx.dat \
    -e 2500 -o run01/analysis/pmf_first_half.xvg

# 後半:
~/.local/share/mamba/envs/gmxcudaenv/bin/gmx wham \
    -it run01/analysis/tpr.dat -ix run01/analysis/pullx.dat \
    -b 2500 -o run01/analysis/pmf_second_half.xvg

# (b) Bootstrap 誤差
python -m abmptools.membrane.pmf \
    --windows-dir run01/windows --analysis-dir run01/analysis \
    --config run01/input/config.json \
    --bootstrap-n 50 --gmx-path "$CUDAENV/bin/gmx"
```

### 5.3 よくある修正

| 症状 | 原因 | 修正 |
|---|---|---|
| PMF が単調増加 | pull rate × 時間 < window 範囲 | `pull_rate_nm_per_ps` を上げる、または `nsteps` を増やす |
| histogram overlap 悪い | window 間隔が粗い | `window_spacing_nm` を細かく (0.25 → 0.1) |
| PMF が時間で大きく変動 | window 時間不足 | `window_nsteps` を増やす (5 ns → 20-50 ns) |
| 膜中心の barrier が異常に高い | 膜変形 (peptide 側鎖が膜内で詰まる) | force constant を下げる、または窓中心付近の MD 時間を延長 |
| 配向が固定 | window 内で peptide が回転していない | REUS 検討 (GENESIS、別パッケージ) |

---

## §6. CHARMM36 backend (任意)

AMBER 路線が動いた後で CHARMM36 路線を試す:

```python
cfg = MembraneConfig(
    backend="charmm36",
    charmm_ff_dir="/abs/path/to/charmm36-jul2022.ff",  # Klauda port
    lipids=[LipidSpec(resname="POPC", n_per_leaflet=32)],
    peptide=PeptideSpec(name="aa5", sequence="AAAAA"),
    ions=IonSpec(cation="Na+", anion="Cl-"),  # 自動で SOD/CLA に翻訳
    output_dir="./run01_charmm",
    ...  # equil / pull / umbrella は同じ
)
```

CHARMM36 force-field の取得方法は [`membrane.md`](membrane.md#charmm36-gromacs-port-の取得) を参照。

商用利用 OK な手順 (CGenFF / CHARMM-GUI 経由は不可) を厳守。

---

## §7. 引用・ライセンス

論文発表時の引用:

| 機能 | 引用 |
|---|---|
| GROMACS | Abraham et al., SoftwareX, 2015 |
| OpenMP/CUDA mdrun | Páll et al., J. Chem. Phys., 2020 |
| AmberTools | Case et al., J. Comput. Chem., 2005 |
| ff19SB | Tian et al., J. Chem. Theory Comput., 2020 |
| Lipid21 | Skjevik et al., J. Phys. Chem. B, 2012 |
| TIP3P-JC ions | Joung & Cheatham, J. Phys. Chem. B, 2008 |
| packmol | Martínez et al., J. Comput. Chem., 2009 |
| packmol-memgen | Schott-Verdugo & Gohlke, J. Chem. Inf. Model., 2019 |
| WHAM | Hub et al., J. Chem. Theory Comput., 2010 |
| CHARMM36 (使用時) | Best 2012 / Klauda 2010 / Pastor 2011 |

商用利用 OK ルール: §1 と [`membrane.md`](membrane.md#ライセンス上のルール-商用利用前提) 参照。

---

## 付録 A: 実測パフォーマンス (RTX 4070 Ti + 4 CPU cores)

```
[poly-Ala 5-mer + POPC 32/leaflet, ~18k atoms]

Stage          MD time    wall time   ns/day
─────────────────────────────────────────────
em             —          ~30 sec     —
nvt (heating)  0.2 ns     ~30 sec     ~570
npt-semiiso    5.0 ns     ~13 min     ~570
pull (NVT)     10.0 ns    ~22 min     ~660
window × 13    5.0 ns/win ~12 min/win ~640
wham           —          ~30 sec     —
─────────────────────────────────────────────
Total          81 ns      ~3 hours
```

CPU 4-core 比 ~4-5×。1 core 比 ~15×。

## 付録 A.2: 実例 PMF (run02、2026-05-03)

上記 §2 の設定 (本番より、13 windows × 5 ns、pull rate=0.001 nm/ps × 10 ns、
合計 81 ns) を実際に実行した結果:

```
=== analysis/pmf.xvg サンプル ===
  z = -1.50 nm  →  E =  -3.0 kJ/mol   (下側 headgroup well)
  z = -1.00 nm  →  E =   6.5 kJ/mol
  z = -0.50 nm  →  E =  47.1 kJ/mol
  z =  0.00 nm  →  E =  83.7 kJ/mol   ← 膜中心 barrier
  z = +0.50 nm  →  E = 107.6 kJ/mol
  z = +1.00 nm  →  E = 122.1 kJ/mol
  z = +1.50 nm  →  E = 130.1 kJ/mol   (上側、sampling 不足で過大)

  PMF range:  z = -1.65..+1.65 nm  /  E = -3.8..+131.0 kJ/mol
  Pull traversed: z = -4.82..+5.43 nm (target 範囲を完全カバー)
  Histogram: 各 window 2476 samples
```

### 物理解釈

- ✅ **膜中心 z=0 で free-energy barrier の極大** (~84 kJ/mol)
- ✅ **z=-1.36 nm の headgroup 領域で well** (-3.8 kJ/mol、溶媒和シェル)
- ⚠ **対称性は不完全** (上下で 30 kJ/mol 差)
  - 5 ns/window では peptide 側鎖配向の thermalize が不十分
  - pull が一方向だったので後半 windows の初期構造が緩和不足
  - 本番 PMF は 20-50 ns/window + 双方向 pull + bootstrap 推奨 (§5.1)

### 完全データの保管場所

| Tier | サイズ | 場所 |
|---|---|---|
| **Light** (input + analysis) | ~13 MB | `abmptools-sample/sample/membrane_us/peptide-polyAla5_POPC32_us_20260503_13win5ns/` の analysis/ 部分 |
| **Medium** (上記 + tpr/gro/log) | ~85 MB | 同 sample dir 全体 (再 wham 解析可、軌跡解析不可) |
| **Full** (上記 + .xtc + .edr + .cpt) | ~3 GB | `OneDrive/abmptools-dump/membrane-us/peptide-polyAla5_POPC32_us_20260503_13win5ns/` |

軌跡解析 (RDF / 配向 / contact map 等) は OneDrive のフルデータから .xtc を取得して実施。

## 付録 B: 詳細リファレンス

- [`membrane.md`](membrane.md) — 全 API / config field / 設計判断
- [`amorphous.md`](amorphous.md) — 姉妹サブパッケージ (3D ランダム box)
- [`architecture.md`](architecture.md) — abmptools 全体構成
- [`licenses_third_party.md`](licenses_third_party.md) — third-party license 一覧
- パッケージ source: `abmptools/membrane/`
- 統合テスト: `tests/integration/run_membrane_us_smoke.sh` (AMBER) /
  `run_membrane_us_charmm_smoke.sh` (CHARMM36)
- 翻訳ユニットテスト: `tests/test_membrane_charmm_translate.py`
