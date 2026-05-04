# CHARMM36 Phase D — peptide-polyAla5 × POPC32 入力サンプル

`abmptools.membrane` v1.17.3 の **CHARMM36 backend** で AMBER vs CHARMM PMF 比較
(Phase D, L9 verification) を実行した時の入力設定。実行結果は別途保存:

- **軽量版** (~123 MB、xtc 等除く、再 wham 解析可): [`abmptools-sample/sample/membrane_us/peptide-polyAla5_POPC32_charmm_us_20260504_13win1ns/`](../../../../abmptools-sample/sample/membrane_us/peptide-polyAla5_POPC32_charmm_us_20260504_13win1ns/)
- **完全版** (~1.4 GB、xtc 含む): OneDrive `abmptools-dump/membrane-us/peptide-polyAla5_POPC32_charmm_us_20260504_13win1ns/`
- **PMF 比較プロット**: [`docs/figures/pmf_compare_amber_charmm.png`](../../../docs/figures/pmf_compare_amber_charmm.png)

## ファイル

```
charmm_phaseD/
└── input/
    └── config_phaseD.json    ← MembraneConfig (CHARMM36 backend 用)
```

## 結果サマリ

| 指標 | 値 |
|---|---|
| Backend | CHARMM36 (Klauda port、charmm36-feb2026_cgenff-5.0) |
| 系 | poly-Ala 5 (ACE-AAAAA-NME) + POPC 32×2 + CHARMM TIP3P + 0.15 M NaCl |
| プロトコル | em → nvt 0.1ns → npt-semi 5ns → pull 10ns → 13 windows × 1ns → WHAM |
| APL | 64.88 ± 0.82 Å² ⭐ (lit 65.0 完全一致) |
| PMF (z=0 barrier) | +97.9 kJ/mol (AMBER reference +86.7 と Δ+11.3 kJ/mol、典型 FF 差範囲内) |
| Wall time | 約 4 時間 (RTX 4070 Ti + 4 CPU、abmptoolsenv + gmxcudaenv) |

## 実行例

```bash
ENV=~/.local/share/mamba/envs/abmptoolsenv
PATH=$ENV/bin:$PATH $ENV/bin/python3 - <<'PY'
from abmptools.membrane import MembraneConfig
from abmptools.membrane.builder import MembraneUSBuilder

config = MembraneConfig.from_json("input/config_phaseD.json")
builder = MembraneUSBuilder(config, output_dir="run_charmm_phaseD")
result = builder.build()
print(f"Build complete: run_charmm_phaseD/")
PY

# 完成後、driver で全 stage 実行
GMX=~/.local/share/mamba/envs/gmxcudaenv/bin/gmx \
NT=4 bash run_charmm_phaseD/run.sh
```

`config_phaseD.json` の `charmm_ff_dir` (line 68) はローカルの ff dir パスに合わせて修正必要:

```json
"charmm_ff_dir": "/path/to/charmm36-feb2026_cgenff-5.0.ff"
```

DL: <https://mackerell.umaryland.edu/charmm_ff.shtml> から `charmm36-feb2026_cgenff-5.0.ff.tgz`。

## 関連 docs

- 操作手順: [`docs/tutorial_membrane_us.md`](../../../docs/tutorial_membrane_us.md) (§6 CHARMM backend)
- 検証階層 (L1-L10): [`docs/membrane.md`](../../../docs/membrane.md) Phase C+ セクション
- Klauda port gotchas: memory `project_membrane_charmm36_gotchas.md`

## ⚠ Sampling artifact について

本 config は `force_constant_kj_mol_nm2: 1000.0` + `window_spacing_nm: 0.25` だが、
これだと histogram overlap が 0.01-0.07 にとどまる (推奨 >0.1)。改善するには
以下のいずれか:

- `force_constant_kj_mol_nm2: 500.0` (k 半減で σ ≈ 0.07 nm に拡大)
- `window_spacing_nm: 0.10-0.15` (window 間隔細分、→ 21-31 windows)
- `window_nsteps: 2_500_000-5_000_000` (5-10 ns/window へ延長)

詳細: AMBER reference (`peptide-polyAla5_POPC32_us_20260503_13win5ns/`) も同じ
default で実行した結果、共通 artifact を確認している (CHARMM 固有ではない)。
