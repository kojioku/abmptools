# AMBER Phase D — peptide-polyAla5 × POPC32 入力サンプル

`abmptools.membrane` の **AMBER backend** (ff19SB + Lipid21 + TIP3P + GAFF2) で
peptide-bilayer umbrella sampling PMF を実行する時の入力設定。
[`../charmm_phaseD/`](../charmm_phaseD/) (CHARMM36 backend) と対をなす
Phase D = L9 verification の AMBER ベースライン側。

CHARMM phaseD と完全に同じプロトコル (13 windows × 1 ns、`window_spacing_nm=0.25`)
で、力場のみ差し替え。AMBER vs CHARMM の PMF 比較 (memory
`project_membrane_charmm36_gotchas.md` 参照) で AMBER 側 reference として使用。

## ファイル

```
amber_phaseD/
└── input/
    └── config_phaseD.json    ← MembraneConfig (AMBER backend 用)
```

config の CHARMM 版との差分:

| Field | CHARMM phaseD | AMBER phaseD |
|---|---|---|
| `backend` | `"charmm36"` | **`"amber"`** |
| `charmm_ff_dir` | charmm36 ff dir のパス | `""` (空、未使用) |
| `amber_protein_ff` | `leaprc.protein.ff19SB` | 同左 |
| `amber_lipid_ff` | `leaprc.lipid21` | 同左 |
| `amber_water_ff` | `leaprc.water.tip3p` | 同左 |
| `output_dir` | `.../membrane_charmm_smoke_run12` | `.../membrane_amber_phaseD_run` |

その他 (系構成、equilibration / pulling / umbrella protocol、box / water /
distance パラメータ) は CHARMM 版と完全一致。両 backend の差分は **力場と
水モデル (AMBER TIP3P の H に LJ なし vs CHARMM-modified TIP3P)** のみ。

## 結果サマリ (memory 由来、Phase D AMBER reference)

| 指標 | 値 |
|---|---|
| Backend | AMBER (ff19SB + Lipid21 + GAFF2 + TIP3P) |
| 系 | poly-Ala 5 (ACE-AAAAA-NME) + POPC 32×2 + AMBER TIP3P + 0.15 M NaCl |
| プロトコル | em → nvt 0.1 ns → npt-semi 5 ns → pull 10 ns → 13 windows × 1 ns → WHAM |
| パイプライン | `packmol-memgen` (bilayer + water + ion + peptide) → `tleap`+`parmed` (prmtop+inpcrd → top+gro) → GROMACS pull/wham |
| APL | ~64-65 Å² (Lipid21 推奨範囲) |
| PMF (z=0 barrier) | **+86.7 kJ/mol** (CHARMM phaseD +97.9 と Δ-11.3 kJ/mol、典型 FF 差) |
| Wall time | 約 4 時間 (RTX 4070 Ti + 4 CPU、abmptoolsenv + gmxcudaenv) |

## 使用ツール (packmol-memgen 経路)

| Stage | ツール | 役割 |
|---|---|---|
| 1. 系構築 | **`packmol-memgen`** (AmberTools 同梱) | bilayer (POPC 32×2) + 水 + Na+/Cl- + peptide (ACE-AAAAA-NME) の初期配置を一発で生成。`MembraneConfig.packmol_memgen_path` で binary 指定 |
| 2. 力場割当 | `tleap` (AmberTools) | ff19SB + Lipid21 + GAFF2 + TIP3P を `loadAmberParams` + `loadpdb` で適用、prmtop + inpcrd を出力 |
| 3. GROMACS 変換 | `parmed` | AMBER prmtop+inpcrd → GROMACS top+gro |
| 4. EQ + Pull + US | `gmx grompp` / `gmx mdrun` | em → nvt → npt-semiisotropic → z-pulling → 13 window US |
| 5. PMF | `gmx wham` | umbrella histogram → free energy profile |

商用利用クリーン: **CGenFF Web server / CHARMM-GUI に依存しない** 設計
([`docs/membrane.md`](../../../docs/membrane.md) のライセンスルール節参照)。
小分子追加が必要な場合は AMBER `antechamber` + `parmchk2` で GAFF2 パラメータ生成。

## 実行例

```bash
ENV=~/.local/share/mamba/envs/abmptoolsenv
PATH=$ENV/bin:$PATH AMBERHOME=$ENV $ENV/bin/python3 - <<'PY'
from abmptools.membrane import MembraneConfig
from abmptools.membrane.builder import MembraneUSBuilder

config = MembraneConfig.from_json("input/config_phaseD.json")
# output_dir はその場で上書きしても良い
config.output_dir = "run_amber_phaseD"
builder = MembraneUSBuilder(config, output_dir=config.output_dir)
result = builder.build()
print(f"Build complete: {config.output_dir}/")
print("Next: cd run_amber_phaseD && bash run.sh")
PY
```

`run_amber_phaseD/run.sh` を実行すると em → nvt → npt-semiisotropic → pull
→ 13 window US → `gmx wham` まで一括実行 (GPU が見えれば自動利用)。
完了後 `analysis/profile.xvg` に PMF が出る。

## 重量級結果の保管先

- **軽量版** (input + analysis のみ、~13 MB):
  `abmptools-sample/sample/membrane_us/peptide-polyAla5_POPC32_us_20260503_13win5ns/` の `analysis/`
- **中量版** (上記 + tpr/gro/log、~85 MB、再 wham 解析可): 同 sample dir 全体
- **完全版** (xtc 含む、~1.4 GB): OneDrive `abmptools-dump/membrane-us/peptide-polyAla5_POPC32_us_20260503_13win5ns/`
- **PMF 比較プロット**: [`docs/figures/pmf_compare_amber_charmm.png`](../../../docs/figures/pmf_compare_amber_charmm.png)

## 関連 docs

- [`docs/membrane.md`](../../../docs/membrane.md) — Subsystem reference (AMBER backend が default、CHARMM36 は opt-in)
- [`docs/tutorial_membrane_us.md`](../../../docs/tutorial_membrane_us.md) — Step-by-step ops (AMBER 経路ベース)
- [`docs/dependencies.md`](../../../docs/dependencies.md) — AmberTools / packmol-memgen / parmed の install 手順
- [`docs/licenses_third_party.md`](../../../docs/licenses_third_party.md) — packmol-memgen / AmberTools / Lipid21 等の license posture
