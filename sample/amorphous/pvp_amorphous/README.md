# PVP amorphous — RDKit で組んだオリゴマーの非晶質系

**ポリビニルピロリドン (PVP)** 5-mer の amorphous 系を組む例。モノマーが
市販の SDF/SMILES で手に入らないので、**`input/build_pvp_oligomer.py` が RDKit で
オリゴマーを構築**してから packmol に渡す。SDF 入力でも SMILES 入力でもない
第 3 の経路の実例。

## 系

| 項目 | 値 |
|---|---|
| 分子 | PVP 5-mer (C31H49N5O5、ピロリドン 5 個/鎖)。`input/build_pvp_oligomer.py` で RDKit 構築 |
| 力場 / 電荷 | OpenFF Sage 2.1.0 + AM1-BCC |
| 規模 | 20 鎖 × 90 atoms = 1800 atoms、box ~2.51 nm、density 1.2 g/cm³ |
| MD | 5-stage (EM → NVT/NPT 500 K → anneal → NPT 300 K prod 500 ps)、101 frames @ 5 ps |

## 再現手順

```bash
cd <abmptools>/sample/amorphous/pvp_amorphous
bash run_sample.sh                                          # oligomer 構築 + OpenFF + AM1-BCC + packmol
cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 8' bash run_all.sh     # 5-stage MD
python wrap_pbc.py                                          # whole 化
```

追跡されるのは**軽量側** (`run_sample.sh` / `input/build_pvp_oligomer.py` / 本 README)
のみ。`build/`・MD 生出力・`md/` 生成スクリプトは再生成できるため追跡しない。

## 備考

- この系は水素結合解析の題材としても使っている (PVP は N-H を持たない 3 級アミドで、
  donor を持たない系の例になる)。解析側のワークフローは本パッケージには
  含まれない (README の "Beyond this package" を参照)。
- 混合系の例は `indomethacin_pvp_asd/`。
