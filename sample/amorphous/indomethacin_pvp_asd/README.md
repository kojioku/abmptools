# Indomethacin + PVP ASD — 2 成分混合系の構築

**インドメタシン (薬物) + PVP (ポリマー)** を共非晶質化する例。
非晶質固体分散体 (ASD = amorphous solid dispersion) を模した **薬物とポリマーの
混合系**で、`sample/amorphous/mixture_json` の低分子どうしの混合に対し、
**分子サイズが大きく異なる 2 成分**を詰めるケースになっている。

## 系

| 項目 | 値 |
|---|---|
| 成分 | indomethacin (41 atoms) × 24 + PVP 5-mer (90 atoms) × 6 |
| 力場 / 電荷 | OpenFF Sage 2.1.0 + AM1-BCC (両成分) |
| 規模 | 30 分子 / 1704 atoms、box ~2.49 nm、density 1.3 g/cm³ |
| MD | 5-stage (EM → NVT/NPT 500 K → anneal → NPT 300 K prod 500 ps)、101 frames @ 5 ps |

PVP オリゴマーは `input/build_pvp_oligomer.py` が RDKit で組む。

## 再現手順

```bash
cd <abmptools>/sample/amorphous/indomethacin_pvp_asd
bash run_sample.sh                                          # 2成分 build (OpenFF + AM1-BCC + packmol)
cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 8' bash run_all.sh     # 5-stage MD
python wrap_pbc.py
```

## 備考

- この系は種間水素結合の解析にも使っている。解析側のワークフローは
  本パッケージには含まれない (README の "Beyond this package" を参照)。
- 単一成分の対照は `indomethacin_amorphous/` と `pvp_amorphous/`。
