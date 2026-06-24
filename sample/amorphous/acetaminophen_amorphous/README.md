# acetaminophen amorphous — COOH を持たない薬物の hbond generic-mode 解析

公開薬物 **acetaminophen (paracetamol, C8H9NO2)** の amorphous MD を作り、
`abmptools.hbond` の **generic mode** で H-bond network を解析する end-to-end の例。

acetaminophen は **カルボキシル基を持たない**(2 級アミド N-H/C=O + フェノール
-OH のみ)。IMC 用の `--classify-mode imc`(COOH 4-species 前提)は使えないため
**`--classify-mode generic`** + element fallback で解析する。OpenFF SMIRNOFF の
per-atom unique atom type は element + bond-graph fallback(default ON)で
O/N/H を判定する(`reference_openff_smirnoff_atom_typing`)。

## 系

| 項目 | 値 |
|---|---|
| 分子 | acetaminophen (PubChem name 取得、3D SDF) |
| 規模 | 64 分子 × 20 atoms = 1280 atoms、box 2.37 nm、density 1.2 g/cm³ |
| FF / 電荷 | OpenFF Sage 2.1.0 / AM1-BCC |
| MD | GROMACS 5-stage anneal、T_high=400 K、production 500 ps(5 ps stride → 101 records) |
| 官能基 | 2 級アミド(amide_donor N-H + amide_O C=O)、フェノール(hydroxyl + hydroxyl_O) |

## 再現手順

```bash
cd <abmptools>/sample/amorphous/acetaminophen_amorphous
bash run_sample.sh                                          # build (~1 min, AM1-BCC)
export OMP_NUM_THREADS=4
cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 4' bash run_all.sh     # 5-stage MD
python wrap_pbc.py && python build_bdf.py                   # xtc -> UDF (101 rec)
cd ../
python -m abmptools.hbond md/05_npt_final.udf \
    --out-prefix output/apap_hbond \
    --classify-mode generic \
    --donor-groups amide_donor,hydroxyl \
    --acceptor-groups amide_O,hydroxyl_O \
    --dt 5.0 --no-colorize --no-copy-uncolored               # ~7 s
```

大規模な multi-record UDF(数百〜数千 record)を解析する場合は **`--record-stride N`** で
N record ごとに間引ける(例: `--record-stride 10 --record-end 1000` で 1001→100 frames)。
`run()` の record ループにのみ作用する hbond 専用オプション。lifetime/τ_HB を使うなら `--dt`
も同倍率でスケールする(間引くと sampled frame が離れ、連続性指標 `continuous` は意味が薄くなる)。

`--donor-groups`/`--acceptor-groups` を手で指定する代わりに **`--classify-mode auto`** も使える。
系の官能基を自動検出して generic 実行する(この系なら amide + hydroxyl を検出 → 上と同じ
`amide_donor,hydroxyl` / `amide_O,hydroxyl_O` 相当を自動選択。明示指定すればそちらが優先)。

本体 repo に追跡されるのは **軽量側**(input SDF / `run_sample.sh` /
`md/build_bdf.py` / `output/apap_hbond_*` の集計 CSV・PNG / 本 README)のみ。
`build/`・MD 生出力(xtc/trr/edr/tpr/gro/log/udf, ~13 MB)は `run_sample.sh`+
MD で再生成できるため追跡しない。

## 結果(101 records、Luzar-Chandler default、`--dt 5.0`)

### pair ごとの平均 H-bond 数 / record

| donor → acceptor | 平均 /rec | d(D...A) mean [Å] | peak [Å] |
|---|---|---|---|
| **hydroxyl → amide_O**(phenol OH → amide C=O) | **39.3** | **2.75** | **2.67** |
| amide_donor → amide_O(amide N-H → C=O) | 26.0 | 3.02 | 2.92 |
| amide_donor → hydroxyl_O | 22.3 | 3.08 | 2.98 |
| hydroxyl → hydroxyl_O | 16.3 | 2.91 | 2.83 |

→ **phenol OH → amide carbonyl が最多かつ最短**(最強の H-bond)。amorphous
acetaminophen は OH と amide の相互供与で密な network を作る。

### lifetime / autocorrelation

- unique pairs: 185、mean occupancy 0.561、longest continuous run 101 frames
- **τ_HB(∫C(t)dt)= 212 ps**(dt=5 ps、Luzar-Chandler unbiased 自己相関）

## 出力ファイル(`output/apap_hbond_*`)

`pair_stats.csv` / `pairs.csv` / `count.png` / `distance_stats.csv` /
`distance_hist.{csv,png}` / `distance_by_class.png` / `distance_angle_2d.png` /
`lifetime.csv` / `autocorr.{csv,png}`。

## 備考

- この sample 作成中に **generic mode で lifetime/autocorr が空になるバグ**を発見・
  修正(`_write_lifetime_outputs` が imc 専用の `hbonds_cc/ca` のみ参照していた →
  generic mode では `hbonds_by_pair_type` を使うよう修正、回帰テスト追加）。
- 関連: `docs/hbond.md`(generic mode / element fallback)、`project_abmptools_hbond`。
