# ATOMセクション
■  MP2D
    - MP2D + NBO + ESP  --mp2d,  --nbo,  --esp
    - MP2D + NBO        --mp2d, --nbo
    - MP2D + ESP        --mp2d, --esp
    - MP2D              --mp2d

■ MP2D以外
    - NBO + ESP     --hf, --esp
    - NBO           --hf, --nbo
    - ESP           --hf, --esp
    - 計算のみ      --hf

-> 計8種類
--mp2d,  --nbo,  --esp
--mp2d, --nbo
--mp2d, --esp
--mp2d
--hf, --nbo --esp
--hf, --nbo
--hf, --esp
--hf

# DIPOLE部分
■ MP2D
■それ以外

上記ものをつかう
--mp2d
--hf

# Monomerセクション
■ MP2
■ HF, CIS
■ MP3
■ LRD

--mp2
--hf
--mp3
--lrd

# Dimer セクション
■ PIEDA on
    ■ solv on
        ■ BSSE ON
            - CIS or HF
            - MP2
            - MP3
            - LRD
        ■ BSSE OFF
            - CIS or HF
            - MP2
            - MP3
            - LRD
    ■ solv OFF
        ■ BSSE ON
            - CIS or HF
            - MP2
            - MP3
            - LRD
        ■ BSSE OFF
            - CIS or HF
            - MP2
            - MP3
            - LRD
■ PIEDA OFF
    ■ solv on
        ■ BSSE ON
            - CIS or HF
            - MP2
            - MP3
            - LRD
        ■ BSSE OFF
            - CIS or HF
            - MP2
            - MP3
            - LRD
    ■ solv OFF
        ■ BSSE ON
            - CIS or HF
            - MP2
            - MP3
            - LRD
        ■ BSSE OFF
            - CIS or HF
            - MP2
            - MP3
            - LRD

for pieda
    for solve
        for bsse
            for ['HF', 'MP2', 'MP3', 'LRD']
