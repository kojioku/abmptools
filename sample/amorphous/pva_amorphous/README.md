# PVA amorphous (for abmptools.hbond generic-mode demo)

ポリビニルアルコール (PVA) 10-mer atactic を amorphous box にして MD を流し、
`abmptools.hbond` の **generic mode** (donor-type × acceptor-type pair stats)
の入力にする sample。COOH を持たない系で `--classify-mode generic
--donor-groups hydroxyl --acceptor-groups hydroxyl_O` が動作することの demo。

## System

| 項目 | 値 |
|---|---|
| 分子 | PVA 10-mer atactic (`CH3-(CH(OH)-CH2)9-CH(OH)-CH3`) |
| SMILES | `CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)C` |
| atoms / mol | 75 (incl. H), 31 heavy |
| MW / mol | 456.57 g/mol |
| 分子数 | 30 (= 約 2,250 atoms) |
| 目標密度 | 1.2 g/cm³ (PVA バルク ~1.25) |
| 力場 | OpenFF Sage 2.1.0 (`openff_unconstrained-2.1.0.offxml`) |
| 電荷 | AM1-BCC |
| T_high | **400 K** (memory `feedback_openff_thigh_for_small_organic` 通り、600 K だと super-critical 化リスク) |
| T_prod | 300 K |

## Build + MD

```bash
cd sample/amorphous/pva_amorphous
bash run_sample.sh
cd md && bash run_all.sh && python wrap_pbc.py
```

5-stage MD (~10-20 min on CPU 8 cores): EM → NVT (T_high) → NPT (T_high)
→ annealing → NPT (T_prod)。出力: `md/test_05_output_rec*.bdf`。

## H-bond 解析 (generic mode)

xtc → multi-record UDF へ変換してから hbond CLI を呼ぶ:

```bash
# 1) xtc -> UDF (101 records)
cd md
python build_bdf.py

# 2) hbond generic-mode 解析 (antechamber patch 不要)
cd ..
python -m abmptools.hbond md/05_npt_final.udf \
    --out-prefix output/pva_hbond \
    --classify-mode generic \
    --donor-groups hydroxyl \
    --acceptor-groups hydroxyl_O \
    --colorize-mode both
```

`--donor-groups`/`--acceptor-groups` を省きたい場合は **`--classify-mode auto`** が使える。
系の官能基を自動検出して generic 実行する(PVA は hydroxyl のみ → `hydroxyl` / `hydroxyl_O`
を自動選択。明示指定すればそちらが優先)。大規模 multi-record は `--record-stride N` で間引ける。

**Element + bond-graph fallback について**: OpenFF Sage は UDF の
``Atom_Type_Name`` に per-atom unique な ``MOL0_X`` を書く (SMIRNOFF が atom
type の概念を持たないため)。``abmptools.hbond`` v1.28+ では fallback (default
ON) が自動で element + bond graph から ``hydroxyl_O`` / ``hydroxyl_H`` 等を
推定するので、**antechamber で GAFF type を別途生成して UDF に patch する
必要はない**。strict mode が必要なら ``--no-element-fallback`` で無効化。

出力:
- `output/pva_hbond_pair_stats.csv`: per-record の `hydroxyl -> hydroxyl_O` 集計
  (n_hbonds, n_uniq_donors, n_uniq_acceptors, ratio_*_busy)
- `output/pva_hbond_pairs.csv`: H-bond ペア一覧
- `output/pva_hbond_count.png`: ペア種別ごとの本数 vs record(`hydroxyl→hydroxyl_O`)
- `output/pva_hbond_diagram.{png,svg}`: 2D 構造式に H-bond サイトを色分け
  (`abmptools 2.4.0`+ で解析時に自動生成、`2.4.1`+ で both=magenta)。**PVA は全 OH が
  donor かつ acceptor → 全 OH がマゼンタ**(下記の色規約と同じ)
- `output/pva_hbond.bdf`: Mol_Name 維持 + Attributes に
  `hbond=Donor/Acceptor/Both` を append (J-OCTA Attribute フィルタ用)
- `output/pva_hbond_action.bdf` + `_show.act`: gourmet autorun
- `output/pva_hbond_show.py`: J-OCTA Python パネル用

可視化:
- OCTA gourmet: `gourmet output/pva_hbond_action.bdf`
- J-OCTA: `output/pva_hbond.bdf` を開いて Python パネルで `output/pva_hbond_show.py`
  を Load → Run

色: Donor=red / Acceptor=cyan / Both=magenta (= 同 atom が donor + acceptor 兼任、
hydroxyl O は OH-H として donate しつつ別 H-bond の acceptor になり得る) /
Candidate=描画 skip。
