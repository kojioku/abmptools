## ABMPTools 簡易マニュアル
2020.04 -  
Author: Koji Okuwaki

実行サンプルはsampleディレクトリにあり  
run.sh, batchrun.sh などでDL, 実行可

### getifiepieda.py
log ファイルから各種IFIE情報を取得するスクリプト  
現状の機能、オプションは以下の通り

#### fragids  
N:1、それに付随する1:1が出る

#### ffmatrix  
frag-fragのマトリックス  
`-ffmatrix i-j k-l`

e.g.) `python $abmptdir/getifiepieda.py --ffmatrix 1-100 101-200 -i ../6lu7-multi-fmolog/6lu7orig_md040j8_163neu-100-hopt-ps-mod_forabmp_192n-2p-24t.log`

#### frag_dist  
あるフラグメントを基準に距離でスクリーニング

#### multi_fragid_dist  
時系列ifie: あるフラグメントを基準に距離でスクリーニング

#### fraginmol  
分子aのx番目のフラグメントと、分子名STRのy番目のフラグメント間の相互作用

#### selectmolid  
分子 i からの距離r以内の相互作用

#### selectmolids  
分子 i-j からの距離r以内の相互作用

#### multi-fragid-dist  
frag i からみた距離r以内の時系列の相互作用

#### multi-fragid-fragid  
frag i とjのIFIEの時系列変化

#### multi-fragid-molname  
frag i と分子AのIFIEの時系列変化

#### multi_fragids-tfmatrix  

`--tfmatrix i-j k-l
-t 開始時間 終了時間 間隔  
-exclude 省く残基  
-np 読み込みの並列数`

fragment i-j, k-l 間の時系列の相互作用の出力  

### pdb2fmo.py
pdbにフラグメント情報を割り当てて、FMO実行可能なajf, pdbセットを作成するスクリプト  
中分子分散体やポリマーのpdb構造に、フラグメントを割り当てる際に主に使用。  
1分子ずつのフラグメント情報を事前に準備する必要有

### udf2fmo.py
udfにフラグメント情報を割り当てて、FMO実行可能なajfを作成するスクリプト

### ajf2config.py
フラグメント情報を読み込んで、pythonの辞書型に保存するスクリプト  
各種フラグメント作成で使用

usage: `python ajf2config.py -i xxx.ajf yyy.ajf`

### ajfserial.py
雛形ajfから連番を作成する  
e.g.) `python ~/fmodpd/abmptools/ajfserial.py -i 1eo8-ff14sb-xxxps-renamed-HF-STO-3G-nbo.ajf -t 100000 200000 5000 -str xxx`

### generateajf.py
指定したオプションのajfファイルを作成する  

e.g.) `python ~/fmodpd/abmptools/generateajf.py -c 1eo8-ff14sb-xxxps-renamed.pdb -cmm -mem 6000 -np 1 -lc NA 1 -rs Na 0.0 -basis STO-3G --method HF`

### pdbmodify.py
PDBの情報を編集するスクリプト  
残基名の変更、残基番号、原子名のつけなおし、など  

e.g.) `python ~/fmodpd/abmptools/pdbmodify.py -mode rename -str ' NA' 'NA ' 'CYX' 'CYS' ' CL' 'CL ' -i \*pdb`

e.g.) `python $abmptdir/pdbmodify.py -mode resnum -aresname -reatm -i cyc.pdb.1`

### addsolvfrag.py
雛形ajfのフラグメント情報に、読み込んだPDBの追加溶媒情報を追加して、新たにajfを作成するスクリプト  
手動分割で溶媒の数が異なっている複数構造の調整の際に使用
