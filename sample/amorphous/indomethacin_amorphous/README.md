# Indomethacin (IMC) amorphous — 多官能基薬物の非晶質系

**インドメタシン** (1-(4-クロロベンゾイル)-5-メトキシ-2-メチルインドール-3-酢酸) の
amorphous 系を組む例。カルボン酸・3 級アミド・メトキシエーテルを併せ持つ分子で、
**官能基が多い薬物でも OpenFF + AM1-BCC がそのまま通る**ことの確認になる。

## 系

| 項目 | 値 |
|---|---|
| 分子 | indomethacin (C19H16ClNO4、MW 357.79)、PubChem `indometacin` から 3D SDF 取得 |
| 力場 / 電荷 | OpenFF Sage 2.1.0 (`openff_unconstrained-2.1.0.offxml`) + AM1-BCC |
| 規模 | 48 分子 × 41 atoms = 1968 atoms、box ~2.80 nm、density 1.3 g/cm³ |
| MD | 5-stage (EM → NVT 500 K → NPT 500 K → anneal → NPT 300 K prod 500 ps)、101 frames @ 5 ps |

## 再現手順

```bash
cd <abmptools>/sample/amorphous/indomethacin_amorphous
bash run_sample.sh                                          # build (OpenFF + AM1-BCC + packmol)
cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 8' bash run_all.sh     # 5-stage MD
python wrap_pbc.py                                          # 分子を whole 化 (05_npt_final_pbc.xtc)
```

## 備考

- この系は水素結合解析の題材としても使っている。解析側のワークフローは
  本パッケージには含まれない (README の "Beyond this package" を参照)。
- 塩素を含む分子なので、力場割当が通ることの確認用にも使える。
