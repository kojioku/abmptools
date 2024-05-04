## ABMPTools 簡易マニュアル
2020.04 -  
Author: Koji Okuwaki
実行サンプルはsampleディレクトリにあり  
run.sh, batchrun.sh などでDL, 実行可


### 各種IFIE情報の取得(getifiepieda)
log ファイルから各種IFIE情報を取得するモジュール

    - 一点構造
        - 指定した"フラグメント"や"分子"を基準に、"距離"や"対象フラグメント, 分子名"を指定したIFIEテーブル出力
        - 指定フラグメント種間のIFIEを取得(同一分子種の分散体などを対象）(fraginmol)
        - 指定フラグメント群間の相互作用マトリックス出力 (ffmatrix)

    - MDサンプリング構造
        - 対象フラグメント間のIFIEの時系列変化 
        - 対象フラグメント-指定した残基種間 IFIE和の時系列出力(水との相互作用を一括計算など)
        - 対象フラグメント-フラグメント群間の時間-IFIEの相互作用マトリックス出力 (1:1, 1:N, N:1, N:Nの時間変化の一括出力, tfmatrix）

### オプション一覧
  -f [FRAG ...], (--frag)
    基準とするフラグメントidを指定。(int) 単体でも可。
  -m MOL, (--mol)
    基準とする分子ID 。フラグメントと分子の対応は、フラグメント間距離から自動認識される
  -d DIST, (--dist)
    基準フラグメントや分子からの距離。基準は -fや-mなどで事前に指定。
  -mn MOLNAME, (--molname)
    PDBに入力された分子名(非タンパク)
  -fi FRAGINMOL FRAGINMOL FRAGINMOL FRAGINMOL, (--fraginmol)
    分子ID, 分子内の対象フラグメントid, 分子名, 分子2内の対象フラグメントID
  -ff FFMATRIX FFMATRIX, (--ffmatrix)
    fragmentマトリックスを出力するモード。 引数1は範囲1 範囲2, 範囲はハイフン区切り
  -mul [MULTI ...], (--multi)
    多サンプルのIFIEを出力するモード
  -z ZP, (-zp, --zp)
    連番指定の際のzeropadding数
  -tf TFMATRIX TFMATRIX, (--tfmatrix)
    時間-フラグメントのmatrix出力
  -t TIME TIME TIME, (--time)
    連番指定の際の時間(ファイル番号)指定(start end interval)
  -np PYNP, (--pynp)
    並列数
  -ex [EXCLUDE ...], (--exclude)
    出力対象にしないフラグメント
  -i INPUT, (--input)
    インプットファイル名
  -ix INPUTX, (--inputx) 
    連番用複数 Input log name (e.g., file-xxx-bbb.log)
  -dimeres, (--dimeres) 
    Dimer-ES範囲内の情報のみを取得
  -nof90, (--nof90so)
    Fortranライブラリを使用しないフラグ
  -nores, (--noresinfo)
    残基名情報を出力結果に記載しない（処理時間短縮)
  -dimene DIMENE DIMENE, (--dimene)
    指定したダイマーのダイマーエネルギーを出力
  -momene MOMENE, (--momene)
    指定したモノマーのモノマーエネルギーを出力
  -imd, (--is_momdim)
    モノマー, ダイマーエネルギーの出力フラグ

### 複数フラグメント間のIFIEマトリックス情報取得(ffmatrix)
- 概要
  指定した残基のi,jマトリックスを出力可能です。
  タンパク質界面や、解析対象のフラグメントが1つでない場合に有効です。
- 入力
  --ffmatrix フラグメント群1 フラグメント群2 
  -i ログファイル名

  logとajf, pdbを同じフォルダに配置し
  python -m abmptools.getifiepieda --ffmatrix 1-100 101-200 -i ログファイル名
  という形で実行可能です。

  - 範囲選択
    --ffmatrix の後に引数が2つ必要で、それぞれ範囲指定が可能です。

    1. 一括の範囲指定（ハイフン)
        例えば --ffmatrix 1-100 101-200 で、1～100, 101～200の間のIFIEや最近接距離情報が出力できます。

    2. 任意の残基指定(カンマ）
        カンマ区切りで個別に指定することも可能です。
        例えば、　--ffmatrix 1,4,5 7,10,11　であれば、　1,4,5と7,10,11 の総当たりのデータを出力できます。


- 出力
    1ファイルにPIEDA情報を全て出せないので、別々のファイルが出力されます。
    Distance-ffmatrix.csvが距離情報のマトリックスで、
    その他PIEDAの各成分やMP2まで含めたIFIEなどがそれぞれファイルごとに出力されています。
    - CT-ffmatrix.csv　CT項
    - Distance-ffmatrix.csv　距離
    - ES-ffmatrix.csv　ES項
    - EX-ffmatrix.csv　EX項
    - HF-ffmatrix.csv　HF項
    - MP2corr-ffmatrix.csv　MP2相関のみ
    - MP2total-ffmatrix.csv　HF + MP2Corr
    - PRMP2corr-ffmatrix.csv PR-MP2 (MP2の過大評価補正をした値）
    - PRMP2total-ffmatrix.csv HF + PRMP2Total

    [img]

- (参考) PPIの議論での利用イメージ

界面の残基を選定-> 近距離で接している相手残基と水を選定 -> IFIEを見る
という手続きかと思いますので、事前に構造解析で対象の残基群を特定いただくか、
ffmatrixモードを使用されるのであれば以下で完了できるかもしれません。
（ざっくりイメージですので、足りない部分がありましたら申し訳ありません)

1. タンパク1　タンパク2間の距離マトリックス出力 → PPI界面残基のID選出
2. タンパク1 PPI界面残基と水の距離マトリックス算出 -> タンパク1 PPI界面付近の水ID選出
3. タンパク2 PPI界面残基と水の距離マトリックス算出 -> タンパク2 PPI界面付近の水ID選出
4. 1-3ででてきた残基番号で OR をとりIFIEマトリックス出力 


-  サンプル
getifiepieda のサンプル入出力一覧を以下にアップいたしました。
sample/getifiepieda/ffmatrix-6lu7/csv　内が結果の出力例となっており、以下のような形で実行した結果となります。
python -m abmptools.getifiepieda --ffmatrix 1-100 101-200 -i xxx.log

e.g.) `python -m abmptools.getifiepieda --ffmatrix 1-100 101-200 -i xxx.log`


#### フラグメント群間IFIE (fragids)
あるフラグメント群(i1-i2, j1-j2)間のIFIEを取得

- 入力
python -m abmptools.getifiepieda --frag i1-i2 j1-j2 -i xxx.log
--frag 基準フラグメント(i1-i2) 相手フラグメント(j1-j2)。 単体も可。範囲指定は-。
-i 入力ログファイル名

- 出力
i1-i2に対するj1-j2 のIFIEがiの各フラグメントごとに1ファイルで出力


#### 基準フラグメント-指定距離間IFIE(frag_dist)
指定フラグメント(i)を基準に、指定した距離(r)内のフラグメントとのIFIEを出力
- 入力
python -m abmptools.getifiepieda --frag i(int) -d r(float) -i xxx.log

- 出力
frag i から指定距離rに含まれるフラグメントとのIFIEを1ファイルで出力


#### 基準分子-指定距離間IFIE(mol-dist)
分子 i からの距離r以内の相互作用を出力
- 入力
python -m getifiepieda --mol i(int) -d r(float) -i xxx.log
--mol 基準分子(i)
-d 距離(r)
-i 入力ログファイル

- 出力


#### 指定フラグメント種間のIFIE (fraginmol)
分子iのj番目のフラグメントと、分子名STRのk番目のフラグメント間の相互作用
python -m abmptools.getifiepieda --fraginmol i j WAT k -i xxx.log


### 多サンプル: 基準フラグメントから指定距離以内の時系列相互作用(multi-fragid-dist)
- 入力
python -m abmptools.getifiepieda --multi 10 -d 8.0 -t 100 3100 1000 -i '["aaa", "-bbb.log"]' --exclude 102 -np 4
--multi 対象フラグメント(int)
-d 距離(r)
-t 開始番号　終了番号 間隔

- 出力


#### selectmolids  
分子 i-j からの距離r以内の相互作用



#### multi-fragid-fragid  
frag i とjのIFIEの時系列変化

    python -m abmptools.getifiepieda --multi 1-100 101-200 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", -hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4
    python -m abmptools.getifiepieda --multi 1-100 101-200 -t 100 3100 1000 -ix 6lu7orig_md040j8_163neuxxx-hopt-ps-mod_forabmp_192n-2p-24t.log --exclude 102 -np 4

    python -m abmptools.getifiepieda --multi 20 --molname WAT -d 8.0 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", "-hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4

#### multi-fragid-molname  
frag i と分子AのIFIEの時系列変化


#### multi_fragids-tfmatrix  

`--tfmatrix i-j k-l
-t 開始時間 終了時間 間隔  
-exclude 省く残基  
-np 読み込みの並列数`

fragment i-j, k-l 間の時系列の相互作用の出力  

#### 構文と引数

'''
usage: e.g)
    # python -m abmptools.getifiepieda --frag 1-10 101-200 -i xxx.log
    # python -m abmptools.getifiepieda --frag 10 -d 0.8 -i xxx.log
    # python -m abmptools.getifiepieda --frag 10 --molname WAT -i xxx.log
    # python -m abmptools.getifiepieda --frag 10 --molname WAT -d 0.8 -i xxx.log
    # python -m abmptools.getifiepieda --mol 1-10 -i xxx.log
    # python -m abmptools.getifiepieda --fraginmol 1-10 2 000 1 -i xxx.log
    # python -m abmptools.getifiepieda --ffmatrix 1-100 101-200 -i xxx.log
    # python -m abmptools.getifiepieda --tfmatrix 1-100 101-200 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", "-hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4

Analysis script for ABINIT-MP log


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


'''
usage: python generateajf.py -i xxx.pdb

generate ABNITMP input file (ajf)

options:
  -h, --help            show this help message and exit
  -i INCOORD, --incoord INCOORD
                        coordinate file (pdb)
  -pb, --solvation      parameter file
  -th THICNV, --thicnv THICNV
                        threshold
  -arad ATMRAD, --atmrad ATMRAD
                        atmrad
  -ajfv AJFVERSION, --ajfversion AJFVERSION
                        ajf version
  -np NPRO, --npro NPRO
                        ajf version
  -nopieda, --nopieda   pieda
  -cmm, --cmm           cmm
  -nocpf, --nocpf       cpfflag
  -cpfv CPFVER, --cpfver CPFVER
                        cpf version
  -basis BASISSET, --basisset BASISSET
                        basis
  -m METHOD, --method METHOD
                        method
  -ml, --mldat          mldat flag
  -mll MLLIMIT, --mllimit MLLIMIT
                        mldat fraglimit
  -disp, --disp         flag disp
  -dg, --dgemm          dgemm
  -rp, --resp           resp
  -nonbo, --nonbo       nonbo
  -mem MEMORY, --memory MEMORY
                        memory
  -lc LIGANDCHARGE LIGANDCHARGE, --ligandcharge LIGANDCHARGE LIGANDCHARGE
                        ligand charge
  -rs RSOLV RSOLV, --rsolv RSOLV RSOLV
                        rsolv
  -ma MANUAL, --manual MANUAL
                        manual table
  -bsse, --bsse         bsse
'''

### pdbmodify.py
PDBの情報を編集するスクリプト  
残基名の変更、残基番号、原子名のつけなおし、など  

e.g.) `python ~/fmodpd/abmptools/pdbmodify.py -mode rename -str ' NA' 'NA ' 'CYX' 'CYS' ' CL' 'CL ' -i \*pdb`

e.g.) `python $abmptdir/pdbmodify.py -mode resnum -aresname -reatm -i cyc.pdb.1`

usage: e.g.)
                python pdbmodify.py -addc 307 312 C -i *.pdb
                python pdbmodify.py -move -p 10 10 10 -aresname -reid -reatm -addc 307 312 C -i *pdb
                python pdbmodify.py -move -mol 3 --into -reres -reatm -addc 307 312 C -i *pdb
                python pdbmodify.py -mode rename -str 001 MRT 002 CD7 003 NA 004 WAT -i *pdb

description

options:
  -h, --help            show this help message and exit
  -i [INPUT ...], --input [INPUT ...]
                        input pdb info
  -move, --move         move
  -p POS POS POS, --pos POS POS POS
                        move position
  -addc ADDCHAIN ADDCHAIN ADDCHAIN, --addchain ADDCHAIN ADDCHAIN ADDCHAIN
                        add chain
  -mol MOL, --mol MOL   move mol
  -into, --into         move
  -aresname, --assignresname
                        assignresname
  -reid, --refreshresid
                        refreshresid
  -reatm, --refreshatmtype
                        refreshatmtype
  -mode MODE, --mode MODE
                        [rfile], [resnum], [rename]
  -str [STRING ...], --string [STRING ...]
                        rename data
  -s SORT SORT, --sort SORT SORT
                        sort

end

### addsolvfrag.py
雛形ajfのフラグメント情報に、読み込んだPDBの追加溶媒情報を追加して、新たにajfを作成するスクリプト  
手動分割で溶媒の数が異なっている複数構造の調整の際に使用

usage: e.g.)
                python addsolvfrag.py -temp 6lu7-covneu-nowat-hinagata0514.ajf -solv HOH WAT NA -i *.pdb

description

options:
  -h, --help            show this help message and exit
  -i [INCOORD ...], --incoord [INCOORD ...]
                        input frag info
  -temp TEMPLATE, --template TEMPLATE
                        template file
  -solv [SOLVMOL ...], --solvmol [SOLVMOL ...]
                        output config file
  -pb, --solvation      parameter file
  -th THICNV, --thicnv THICNV
                        threshold
  -arad ATMRAD, --atmrad ATMRAD
                        atmrad
  -ajfv AJFVERSION, --ajfversion AJFVERSION
                        ajf version
  -np NPRO, --npro NPRO
                        ajf version
  -nopieda, --nopieda   pieda
  -cmm, --cmm           cmm
  -nocpf, --nocpf       cpfflag
  -cpfv CPFVER, --cpfver CPFVER
                        cpf version
  -basis BASISSET, --basisset BASISSET
                        basis
  -m METHOD, --method METHOD
                        method
  -ml, --mldat          mldat flag
  -mll MLLIMIT, --mllimit MLLIMIT
                        mldat fraglimit
  -disp, --disp         flag disp
  -dg, --dgemm          dgemm
  -rp, --resp           resp
  -nonbo, --nonbo       nonbo
  -mem MEMORY, --memory MEMORY
                        memory
  -lc LIGANDCHARGE LIGANDCHARGE, --ligandcharge LIGANDCHARGE LIGANDCHARGE
                        ligand charge
  -rs RSOLV RSOLV, --rsolv RSOLV RSOLV
                        rsolv
  -ma, --manual         manual table
  -bsse, --bsse         bsse

end

