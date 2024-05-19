## ABMPTools マニュアル

2020.04 -
Author: Koji Okuwaki

### Overview
フラグメント分子軌道計算プログラムABINIT-MPの解析支援ツールです。

- IFIE解析支援(getifiepieda)
    フラグメント間相互作用エネルギー(IFIE, PIEDA)の解析機能

- CPF parser(cpfmanager)
    - cpf残基切り出し、バージョンコンバート(convertcpf)
    - 動的平均cpf作成機能(generate_difie)

- logからのcpf生成(log2cpf)

- FMO実行支援
    - PDB読み込み(pdbread)
    - 雛形ajf生成(generateajf)
    - 非タンパク系のPDB, フラグメント情報の処理
        - 一分子のフラグメント情報の全分子への割り当て: 分子集合体向け(ajf2config)
        - 任意のサイズでの分子系の切り出し＋フラグメント割り当て(pdb2ajf)
        - OCTA(COGNAC) MD結果のpdb化 (udf2pdb)

- 各種ファイル変換
    - 結晶構造(cif)ファイルからFMO用ファイルへの変換
    - OCTA udf ファイルからFMO用ファイルへの変換(udf2ajf, udf2pdb)
    - 古典MD後からFMO用ファイルへの変換


### 各種IFIE情報の取得(getifiepieda)
log ファイルから各種IFIE情報を取得するモジュール
主に以下の機能を提供しています。

- 一点構造
    - 指定した"フラグメント"や"分子"を基準に、"距離"や"対象フラグメント, 分子名"を指定したIFIEテーブル出力
    - 指定フラグメント種間のIFIEを取得(同一分子種の分散体などを対象）(fraginmol)
    - 指定フラグメント群間の相互作用マトリックス出力 (ffmatrix)

- MDサンプリング構造
    - 対象フラグメント間のIFIEの時系列変化
    - 対象フラグメント-指定した残基種間 IFIE和の時系列出力(水との相互作用を一括計算など)
    - 対象フラグメント-フラグメント群間の時間-IFIEの相互作用マトリックス出力 (1:1, 1:N, N:1, N:Nの時間変化の一括出力, tfmatrix）


#### オプション一覧

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

#### 複数フラグメント間のIFIEマトリックス情報取得(ffmatrix)
- 概要
  - 指定した残基のi,jマトリックスを出力可能です。
  - タンパク質界面や、解析対象のフラグメントが1つでない場合に有効です。
- 入力
  - `python -m abmptools.getifiepieda --ffmatrix i1-i2 j1-j2 ログファイル名`
  - logとajf, pdbを同じフォルダに配置して実施
  - --ffmatrix 範囲選択
  - 引数が2つ必要で、それぞれ単体での指定のほか、範囲指定が可能です。

    1. 一括の範囲指定（ハイフン)
        - 例えば --ffmatrix 1-100 101-200 で、1～100, 101～200の間のIFIEや最近接距離情報が出力できます。

    2. 任意の残基指定(カンマ）
        - カンマ区切りで個別に指定することも可能です。
        - 例えば、　--ffmatrix 1,4,5 7,10,11　であれば、　1,4,5と7,10,11 の総当たりのデータを出力できます。
  - -i ログファイル名

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


- 出力例

|-------|SER1(1)|GLY2(2)|PHE3(3)|ARG4(4)|LYS5(5)|MET6(6)|ALA7(7)|PHE8(8)|PRO9(9)|SER10(10)|
|-------|-------|-------|-------|-------|-------|-------|-------|-------|---------|-------|
|SER1(1)|0.0000|0.0000|-0.6369|33.3948|17.3142|-0.2083|0.8635|-1.4288|0.5616|0.3865|
|GLY2(2)|0.0000|0.0000|0.0000|-0.1801|0.2554|0.1650|-0.1017|0.0314|0.0596|-0.0364|
|PHE3(3)|-0.6369|0.0000|0.0000|0.0000|2.6525|-0.0665|0.1443|0.2930|-0.1663|-0.0082|
|ARG4(4)|33.3948|-0.1801|0.0000|0.0000|0.0000|-2.1988|3.1808|-1.2092|0.2981|0.6150|
|LYS5(5)|17.3142|0.2554|2.6525|0.0000|0.0000|0.0000|-3.2486|2.4843|-0.4618|0.2761|
|MET6(6)|-0.2083|0.1650|-0.0665|-2.1988|0.0000|0.0000|0.0000|-3.9709|0.4826|0.0157|
|ALA7(7)|0.8635|-0.1017|0.1443|3.1808|-3.2486|0.0000|0.0000|0.0000|-3.9012|0.4255|
|PHE8(8)|-1.4288|0.0314|0.2930|-1.2092|2.4843|-3.9709|0.0000|0.0000|0.0000|-3.6327|
|PRO9(9)|0.5616|0.0596|-0.1663|0.2981|-0.4618|0.4826|-3.9012|0.0000|0.0000|0.0000|
|SER10(10)|0.3865|-0.0364|-0.0082|0.6150|0.2761|0.0157|0.4255|-3.6327|0.0000|0.0000|


- (参考) PPIの議論での利用イメージ

界面の残基を選定-> 近距離で接している相手残基と水を選定 -> IFIEを見る
という手続きかと思いますので、
- 事前に構造解析で対象の残基群を特定いただくか、
- ffmatrixモードを使用されるのであれば以下で完了できる

1. タンパク1　タンパク2間の距離マトリックス出力 → PPI界面残基のID選出
2. タンパク1 PPI界面残基と水の距離マトリックス算出 -> タンパク1 PPI界面付近の水ID選出
3. タンパク2 PPI界面残基と水の距離マトリックス算出 -> タンパク2 PPI界面付近の水ID選出
4. 1-3ででてきた残基番号で OR をとりIFIEマトリックス出力

-  サンプル
sample/getifiepieda/ffmatrix に配置
run.shで実行可能です。


#### フラグメント群間IFIE (fragids)
あるフラグメント群(i1-i2, j1-j2)間のIFIEを取得

- 入力
    `python -m abmptools.getifiepieda --frag i1-i2 j1-j2 -i xxx.log`
    - --frag 基準フラグメント(i1-i2) 相手フラグメント(j1-j2)。 単体も可。範囲指定は-。
    - -i 入力ログファイル名

- 出力
    - log-fragi-fragj1-j2-ifie.csv
        - i1-i2に対するj1-j2 のIFIEがiの各フラグメントごとに1ファイルで出力

| |I|J|DIST|DIMER-ES|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|02J307(307)|SER1(1)|40.376438|1|-0.070281|0.0|0.0|0.0|0.0|0.0|||||
|1|02J307(307)|GLY2(2)|37.224025|1|-0.001255|0.0|0.0|0.0|0.0|0.0|||||
|2|02J307(307)|PHE3(3)|29.580733|1|-0.003137|0.0|0.0|0.0|0.0|0.0|||||
|3|02J307(307)|ARG4(4)|29.78684|1|-0.1123242|0.0|0.0|0.0|0.0|0.0|||||
|4|02J307(307)|LYS5(5)|21.19564|1|-0.2660640|0.0|0.0|0.0|0.0|0.0|||||
|5|02J307(307)|MET6(6)|26.868553|1|-0.002510|0.0|0.0|0.0|0.0|0.0|||||

    - log-fragi1-i2-fragj1-j2n-1sum-ifie.csv
        - i1-i2に対するj1-j2 のN:1 IFIEをまとめて出力

|   |I|J|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|
|---|---|---|-------|--------|--------|------|----|----|--|--|------|-------|---|
|0|02J307(307)-010312(311)|SER1(1)|-0.033885|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|
|1|02J307(307)-010312(311)|GLY2(2)|0.0018825|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|
|2|02J307(307)-010312(311)|PHE3(3)|0.0006275|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|
|3|02J307(307)-010312(311)|ARG4(4)|-0.080321|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|
|4|02J307(307)-010312(311)|LYS5(5)|-0.299949|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|
|5|02J307(307)-010312(311)|MET6(6)|0.0075301|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0|


#### 基準フラグメント-指定距離間IFIE(frag_dist)
指定フラグメント(i)を基準に、指定した距離(r)内のフラグメントとのIFIEを出力
- 入力
    `python -m abmptools.getifiepieda --frag i(int) -d r(float) -i xxx.log`

- 出力
    - frag i から指定距離rに含まれるフラグメントとのIFIEを1ファイルで出力
    - log-i-ifie_distr.csv

- サンプル
    - sample/getifiepieda/frag_dist

| |I|J|DIST|DIMER-ES|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|
|-|-|-|----|--------|-------|--------|--------|------|----|----|--|--|------|-------|-------|
|0|PJE311(310)|LEU27(27)|3.661634|0|0.125502|-0.171938|-0.166918|-0.125502|-0.102284|-0.163780|0.187172|0.004672|-0.066598|-0.171789|0.000221|
|1|PJE311(310)|HIS41(41)|2.839272|0|1.250626|-1.506650|-1.433859|-1.254392|-1.127635|-1.125125|0.472950|1.778051|-1.000187|-1.506963|0.006989|
|2|PJE311(310)|PHE140(140)|2.675780|0|-0.088479|-2.649345|-2.532001|-2.158005|-1.912022|-2.073919|-1.104702|2.332716|-1.316674|-2.649320|-0.012656|
|3|PJE311(310)|LEU141(141)|2.057548|0|-8.572408|-3.326428|-3.231674|-2.605420|-2.245229|-2.815008|-11.569676|5.224549|-2.227418|-3.326582|0.023936|
|4|PJE311(310)|ASN142(142)|2.491730|0|3.468245|-5.085337|-4.885162|-4.090107|-3.592492|-4.086342|1.724414|3.380922|-1.637174|-5.085337|0.013705|
|5|PJE311(310)|GLY143(143)|3.065182|0|0.512048|-0.560366|-0.542168|-0.436119|-0.374623|-0.478790|0.224691|0.051412|0.235970|-0.560092|0.000268|


#### 基準分子-指定距離間IFIE(mol-dist)
分子 i からの距離r以内にあるフラグメントの相互作用を出力
- 入力
    `python -m getifiepieda --mol i(int) -d r(float) -i xxx.log`
    - --mol 基準分子(i)  'i1-i2' の形で、複数分子i1-i2に対して一括で出力することも可能
    - -d 距離(r)
    - -i 入力ログファイル

- 出力
- xxx_ifie-fragmol-molidi-idistr.csv

||I|J|DIST|DIMER-ES|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|
|-|-|-|----|--------|-------|--------|--------|------|----|----|--|--|------|-------|-------|
|0|NA4(17)|CD71(1)|17.064399|1|-13.371600|0.000000|0.000000|0.000000|0.000000|0.000000|-13.371664|0.000000|0.000000|0.000000|0.000000|
|1|NA4(17)|CD71(2)|13.493569|1|-1.123870|0.000000|0.000000|0.000000|0.000000|0.000000|-1.123905|0.000000|0.000000|0.000000|0.000000|
|2|NA4(17)|CD71(3)|9.684535|1|-1.576931|0.000000|0.000000|0.000000|0.000000|0.000000|-1.576847|0.000000|0.000000|0.000000|0.000000|
|3|NA4(17)|CD71(4)|2.331214|0|-119.741367|-4.147210|-3.917542|-3.661518|-3.418044|-2.671936|-118.224746|3.668100|-5.184485|-4.147445|0.161143|
|4|NA4(17)|CD71(5)|4.620488|0|-3.423064|-1.019075|-0.963855|-0.884161|-0.816390|-0.687750|-2.011678|-0.003013|-1.408399|-1.018930|0.026662|
|5|NA4(17)|CD71(6)|7.365636|1|2.518823|0.000000|0.000000|0.000000|0.000000|0.000000|2.518628|0.000000|0.000000|0.000000|0.000000|

- xxx_ifiemol-mol-molidi-idistr.csv

||I|J|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|PR_TYPE1|
|-|-|-------|--------|--------|------|----|----|--|--|------|-------|-------|--------|-|
|0|[17]|"[1,2,3,4,5,6,7,8,9,10,11,12,13,14]"|-262.205481|-5.166286||-4.545679|-4.234434|-3.359686|-259.276920|3.665087|-6.592884|-5.166375|0.187805|-4.881397|
|1|[18]|"[1,2,3,4,5,6,7,8,9,10,11,12,13,14]"|-249.266235|24.695010||26.160244|26.892548|7.053207|-242.924904|0.485723|-6.827464|24.694368|0.174430|24.529975|
|2|[22]|"[1,2,3,4,5,6,7,8,9,10,11,12,13,14]"|-15.416654|-2.046936||-1.573794|-1.337223|-1.791540|-20.037975|7.003602|-2.381964|-2.046723|0.022800|-1.986068|
|3|[25]|"[1,2,3,4,5,6,7,8,9,10,11,12,13,14]"|-2.684486|-3.071659||-2.620480|-2.395204|-2.163653|-2.181071|2.422642|-2.925730|-3.071805|-0.026609|-2.951177|
|4|[27]|"[1,2,3,4,5,6,7,8,9,10,11,12,13,14]"|-12.932971|-2.766689||-2.151103|-1.842995|-2.373869|-22.762547|13.364711|-3.535093|-2.766748|0.025231|-2.690761|


- xxx_ifiesummol-mol-molidi-idistr.csv

| |HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|
|-|-|----|--------|-------|--------|--------|------|----|----|--|--|
|0|-1111.631862|-309.681596|-299.146339|-241.282431|-207.073122|-264.691674|-1527.736941|746.607169|-330.500542|-309.686615|1.024979|


#### 指定フラグメント種間のIFIE (fraginmol)
- 分子iのj番目のフラグメントと、分子名STRのk番目のフラグメント間の相互作用
- `python -m abmptools.getifiepieda --fraginmol i j WAT k -i xxx.log`
- 'ii-i2' の形で、複数分子i1-i2に対して一括で出力することも可能



#### 多サンプル: 基準フラグメントから指定距離以内の時系列相互作用(multi-fragid-dist)
- 入力
    `python -m abmptools.getifiepieda --multi 10 -d 8.0 -t 100 3100 1000 -i '["aaa", "-bbb.log"]' --exclude 102 -np 4`
    - --multi 対象フラグメント(int)
    - -d 距離(r)
    - -t 開始番号　終了番号 間隔

- 出力
    - fragi-distr-ifiedt.csv 対象のフラグメント単体ペアのIFIE (時刻ごと)

|-|I|J|DIST|DIMER-ES|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|TIMES|
|-|-|-|----|--------|-------|--------|--------|------|----|----|--|--|------|-------|-------|-----|
|0|PJE311(310)|LEU27(27)|3.661634|0|0.125502|-0.171938|-0.166918|-0.125502|-0.102284|-0.163780|0.187172|0.004672|-0.066598|-0.171789|0.000221|100|
|1|PJE311(310)|HIS41(41)|2.839272|0|1.250626|-1.506650|-1.433859|-1.254392|-1.127635|-1.125125|0.472950|1.778051|-1.000187|-1.506963|0.006989|100|
|2|PJE311(310)|PHE140(140)|2.675780|0|-0.088479|-2.649345|-2.532001|-2.158005|-1.912022|-2.073919|-1.104702|2.332716|-1.316674|-2.649320|-0.012656|100|
|3|PJE311(310)|LEU141(141)|2.057548|0|-8.572408|-3.326428|-3.231674|-2.605420|-2.245229|-2.815008|-11.569676|5.224549|-2.227418|-3.326582|0.023936|100|
|4|PJE311(310)|ASN142(142)|2.491730|0|3.468245|-5.085337|-4.885162|-4.090107|-3.592492|-4.086342|1.724414|3.380922|-1.637174|-5.085337|0.013705|100|
|5|PJE311(310)|GLY143(143)|3.065182|0|0.512048|-0.560366|-0.542168|-0.436119|-0.374623|-0.478790|0.224691|0.051412|0.235970|-0.560092|0.000268|100|
|6|PJE311(310)|SER144(144)|2.572440|0|0.332580|-2.026856|-1.941514|-1.658508|-1.474647|-1.571284|0.859822|0.611152|-1.138180|-2.026809|0.005500|100|
|7|PJE311(310)|CYS145(145)|0.000000|0|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000||||||100|
|8|PJE311(310)|GLY146(146)|3.515848|0|1.091239|-0.495733|-0.473142|-0.402861|-0.356425|-0.390938|0.954076|0.027657|0.109643|-0.496003|0.001026|100|
|9|PJE311(310)|HIS163(163)|2.103530|0|-9.427076|-4.167918|-3.992843|-3.307603|-2.877131|-3.440635|-12.660462|5.725049|-2.491843|-4.168072|-0.025862|100|
|10|PJE311(310)|HIS164(164)|3.957990|0|-1.525476|-0.472515|-0.449297|-0.378388|-0.330698|-0.384036|-1.165180|-0.000111|-0.359885|-0.472670|0.001293|100|

    - fragi-distr-ifiesum.csv 対象のフラグメント単体ペアのIFIE (時刻ごとの合計)

|-|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|DimerEnergy(2-1)|HF-BSSE|MP2-BSSE|MonomerEnergy(1)|
|-|-------|--------|--------|------|----|----|--|--|------|-------|-------|----------------|-------|--------|---------------|
|100|-51.020917|-42.938594|-41.391783|-34.376227|-30.094102|-34.829916|-79.133890|59.399802|-31.287967|-42.939307|0.086640||0.000000|0.000000|
|1100|-45.781212|-40.052050|-38.429938|-32.244577|-28.342095|-32.117192|-57.886142|38.425953|-26.321731|-40.050398|0.030896||0.000000|0.000000|
|2100|-35.330669|-41.331542|-39.777829|-32.725249|-28.423044|-34.262647|-73.418502|75.503551|-37.416453|-41.331105|0.039566||0.000000|0.000000|
|3100|-51.484019|-32.736544|-31.472740|-26.209818|-22.944258|-26.556203|-61.744878|34.014306|-23.752203|-32.737827|0.052560||0.000000|0.000000|

#### multi-fragid-fragid
frag i とj(単体指定)のIFIEの時系列変化

    python -m abmptools.getifiepieda --multi 1 101 -t 100 3100 1000 -i '["1l2y-", "-hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4
    python -m abmptools.getifiepieda --multi 1 101 -t 100 3100 1000 -ix 1l2y-xxx-hopt-ps-mod_forabmp_192n-2p-24t.log --exclude 102 -np 4

|-|I|J|DIST|DIMER-ES|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|
|-|-|-|----|--------|-------|--------|--------|------|----|----|--|--|------|-------|-------|
|100|PJE311(310)|HIS41(41)|2.839272|0|1.250626|-1.506650|-1.433859|-1.254392|-1.127635|-1.125125|0.472950|1.778051|-1.000187|-1.506963|0.006989|
|1100|PJE311(310)|HIS41(41)|2.872523|0|0.801330|-1.345380|-1.282629|-1.108809|-0.990838|-1.027233|0.434798|1.292914|-0.926495|-1.345674|0.005414|
|2100|PJE311(310)|HIS41(41)|3.541479|0|-0.078439|-1.024096|-0.973267|-0.834588|-0.739834|-0.800075|0.241795|0.289265|-0.609326|-1.023895|0.004512|
|3100|PJE311(310)|HIS41(41)|3.020180|0|1.555596|-1.914532|-1.783382|-1.585717|-1.421309|-1.445154|0.997147|2.091160|-1.532707|-1.914223|0.010777|
|4100|PJE311(310)|HIS41(41)|3.286906|0|-0.886671|-1.080571|-1.024723|-0.870356|-0.764934|-0.865963|-0.678428|0.419987|-0.628261|-1.080443|0.004824|
|5100|PJE311(310)|HIS41(41)|2.743530|0|2.285390|-2.110942|-2.067016|-1.737574|-1.551204|-1.615210|0.618316|3.574962|-1.907958|-2.110877|0.009481|


#### multi-fragid-molname
frag i と分子AのIFIEの時系列変化

    python -m abmptools.getifiepieda --multi 20 --molname WAT -d 8.0 -t 100 3100 1000 -i '["6lu7orig_md040j8_163neu", "-hopt-ps-mod_forabmp_192n-2p-24t.log"]' --exclude 102 -np 4

|-|HF-IFIE|MP2-IFIE|PR-TYPE1|GRIMME|JUNG|HILL|ES|EX|CT-mix|DI(MP2)|q(I=>J)|
|-|-------|--------|--------|------|----|----|--|--|------|-------|-------|
|1700|-60.841441|-35.270428|-34.260137|-27.705173|-23.923173|-29.686848|-130.209592|102.751008|-33.378411|-35.270256|0.123733
|2200|-84.277039|-31.714959|-30.613680|-25.431706|-22.287256|-25.647569|-114.673609|56.683863|-26.290034|-31.716089|0.099634
|2700|-105.633070|-29.253239|-28.389158|-22.798676|-19.573277|-24.987429|-170.539890|98.662858|-33.756243|-29.253696|0.138081
|3200|-116.114989|-28.104269|-27.287252|-22.173676|-19.210576|-23.452541|-162.763167|75.261902|-28.611254|-28.103367|0.086623
|3700|-120.810643|-23.484544|-22.793656|-18.474508|-15.969490|-19.717604|-159.616780|64.032750|-25.227879|-23.485234|0.110387
|4200|-139.041677|-30.624347|-29.770934|-23.596868|-20.088462|-26.703668|-206.489264|107.279775|-39.835033|-30.622384|0.109662


#### multi_fragids-tfmatrix
    `--tfmatrix i-j k-l
    -t 開始時間 終了時間 間隔
    -exclude 省く残基
    -np 読み込みの並列数`

- 出力
    - fragment iとk-l 間の時系列の相互作用の出力(i-jごとに1ファイル)

|---|SER1(1)|GLY2(2)|PHE3(3)|ARG4(4)|LYS5(5)|MET6(6)|ALA7(7)|PHE8(8)|PRO9(9)|
|---|-------|-------|-------|-------|-------|-------|-------|-------|-------|
|100|-0.033885|0.001882|0.000627|-0.080321|-0.299949|0.007530|0.008785|-0.060240|0.026982|
|1100|-0.249748|0.004392|0.026982|-0.333835|-0.655119|0.065260|0.010040|-0.082203|0.060240|


### generateajf
#### 概要
- 指定したオプションのajfファイルを作成する

e.g.) `python -m abmptools.generateajf -c 1eo8-ff14sb-xxxps-renamed.pdb -cmm -mem 6000 -np 1 -lc NA 1 -rs Na 0.0 -basis STO-3G --method HF`


    usage: python -m abmptools.generateajf -i xxx.pdb

    generate ABNITMP input file (ajf)

    options:
      -h, --help            show this help message and exit
      -i INCOORD, --incoord INCOORD
                            入力PDBファイル
      -pb, --solvation      parameter file
                            溶媒効果計算の実施フラグ
      -th THICNV, --thicnv THICNV
                            溶媒効果計算時の収束閾値
      -arad ATMRAD, --atmrad ATMRAD
                            溶媒効果計算時の原子半径
      -ajfv AJFVERSION, --ajfversion AJFVERSION
                            ajf version
      -np NPRO, --npro NPRO
                            フラグメントあたりに割り当てるMPI数
      -nopieda, --nopieda   pieda PIEDA計を無効にするフラグ
      -cmm, --cmm           cmm
                            多重極展開による遠距離フラグメント計算時ダイマーES高速化を有効にするフラグ
      -nocpf, --nocpf       cpfflag
                            CPFを出力しないフラグ
      -cpfv CPFVER, --cpfver CPFVER
                            cpf version
      -basis BASISSET, --basisset BASISSET
                            基底関数
      -m METHOD, --method METHOD
                            method
      -ml, --mldat          mldat flag
                            mldatファイルの出力
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
                            BSSE計算フラグ

#### サンプル
- sample/generateajf
- run.sh にて実行可
- ABINIT-MPのinputオプションは[Web版マニュアル](https://fmodd.jp/member_contents/manual_ABINIT-MP/)をご参照ください。

### cpfmanager
#### 概要
- cpfをpythonコマンドラインで読み込み情報を抽出するモジュール

#### 実行
    python
    # 環境ロード
    import abmptools
    cpf = abmptools.CPFManager()

    # ファイル読み込み
    cpf.parse("xxxin.cpf")

    # 主なクラス変数
    cpfver バージョン
    atominfo 原子の情報
    fraginfo フラグメント情報
    condition 計算条件
    static_data 簡易統計所法
    mominfo モノマー結果
    diminfo ダイマー結果
    labels cpfラベル情報

    >>> print(cpf.atominfo)
        alabels elems elemtypes resnames  resids  fragids  xcoords  ycoords  zcoords chainids optflags    MUL-HF
    0         1     N       N        GLY       1        1    0.162   -0.202    0.000                 1 -0.397224
    1         2     C       CA       GLY       1        1    1.612   -0.031    0.000                 1 -0.038224
    2         3     C       C        GLY       1        2    1.985    1.432    0.000                 1  0.302354
    3         4     O       O        GLY       1        2    1.137    2.319    0.000                 1 -0.310874
    4         5     H      1H        GLY       1        1   -0.293    0.722    0.000                 1  0.177669
    5         6     H      2H        GLY       1        1   -0.122   -0.724   -0.841                 1  0.155283
    6         7     H      1HA       GLY       1        1    2.055   -0.519    0.888                 1  0.072199
    7         8     H      2HA       GLY       1        1    2.055   -0.519   -0.888                 1  0.056519
    8         9     N       N        GLY       2        2    3.277    1.700    0.000                 1 -0.355555
    9        10     C       CA       GLY       2        2    3.784    3.070    0.000                 1 -0.017052
    10       11     C       C        GLY       2        3    3.412    3.784   -1.278                 1  0.295201
    11       12     O       O        GLY       2        3    3.055    4.958   -1.286                 1 -0.277340

    >>> print(cpf.static_data)
    {'nuclear_repulsion_energy': 1889.49334720085, 'total_electronic_energy': -2985.05965937836, 'total_energy': -1095.56631217751, 'natom': 38, 'nfrag': 5, 'ndimer': 10, 'ntrimer': '0', 'ntetramer': '0'}

    >>> print(cpf.labels)
    {'charge': ['MUL-HF'], 'DPM': ['DPM-HF-X', 'DPM-HF-Y', 'DPM-HF-Z'], 'monomer': ['NR', 'HF'], 'dimer': ['NR', 'HF', 'ES', 'EX', 'CT', 'DQ']}

    >>> print(cpf.mominfo)
       fragi  DPM-HF-X  DPM-HF-Y  DPM-HF-Z          NR          HF
    0      1  0.785971  2.812967 -1.326046   32.173322 -111.217901
    1      2  1.450489 -3.348173 -1.614646  105.581624 -294.954214
    2      3 -0.524705 -3.573679  1.565234  105.580761 -294.944053
    3      4  2.456669 -1.604684  2.801527  105.592376 -294.954885
    4      5  3.686567 -4.789664 -0.248676  302.972089 -692.375396

    >>> print(cpf.diminfo)
       fragi  fragj  min-dist          NR          HF         ES        EX        CT         DQ
    0      1      2  0.000000   89.835821 -104.597691 -14.464324 -0.170246 -0.127299  10.081019
    1      1      3  6.434985   48.115142  -48.112201   0.002926  0.000020 -0.000004   0.000026
    2      2      3  0.000000  151.596219 -166.354359 -14.460761 -0.171096 -0.126284  10.021917
    3      1      4  6.423184   47.234065  -47.232722   0.001358  0.000038 -0.000053   0.000029
    4      2      4  5.209663  103.751130 -103.749010   0.002019  0.000735 -0.000634   0.001864
    5      3      4  0.000000  151.609836 -166.363653 -14.456629 -0.170776 -0.126411  10.018572
    6      1      5  6.693533   69.358996  -69.356242   0.002788  0.000003 -0.000038  -0.000007
    7      2      5  3.979030  159.493576 -159.501166  -0.007162  0.004484 -0.004912   0.014454
    8      3      5  5.209633  181.521410 -181.515319   0.005989  0.000747 -0.000645   0.001845
    9      4      5  0.000000  235.076980 -249.830849 -14.456895 -0.171821 -0.125153  10.010322

    # 変数情報出力
    >>> print(vars(cpf))
    >>> print(dirs(cpf))

    # ファイル書き込み
    cpf.write('Title', 'xxxout.cpf')


### リスト出力 cpf2ifielist
- cpfからIFIEのリストを出力するモジュール

#### 実行
    python -m abmptools.cpf2ifielist -i hoge.cpf -f 1 834 --and
    -i 入力cpf名                            # 入力cpf名の指定
    -f 開始フラグメント　終了フラグメント	# フラグメント範囲指定
    --and   						        # フラグメント指定を＆にする
        - フラグメントペアの両方が指定範囲の際に出力する(and)
        - フラグメントペアの片方が指定範囲に含まれていれば出力する(or, デフォルト)

    入力ファイル名: hoge.cpf
    出力ファイル名:
    フラグメント番号範囲: [1, 834]
    andflag: True
    start read hoge.cpf
    cpfver 10
    atom section
    fragment section
    dimer distance section
    dipole moment section
    condition and static data section
    monomer section
    dimer section
    Read CPF file finished.
    Filtering...
    Writing output file...
    Write output file finished: hoge_filtered.csv

#### 出力

|fragi|fragj|min-dist|ES|DI|PR-MP2|SCS-MP2(Grimme)|EX|CT|DQ|MP2-Total
|---|---|---|---|---|---|---|---|---|---|---|
0|2|1|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000
1|3|1|3.830886|53.492467|-0.318749|-0.343800|-0.253076|-0.006184|-0.283089|-0.000217|52.884446
2|3|2|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000
3|4|1|6.872217|3.587648|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|3.587648
4|4|2|3.391023|-1.807980|-1.127326|-1.069723|-0.938376|0.007501|-1.557434|0.003896|-4.485238
5|4|3|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000|0.000000


### log2cpf
#### 概要
- logからcpfを生成するモジュール
- HFとMP2のみ対応

#### 実行
    python -m abmptools.log2cpf -i 入力log名 -o cpf名(任意)

#### サンプル
    sample/log2cpf内にサンプルデータあり
    bash run.sh で実行


### CPFの整形(コンパクト化) convertcpf
- 指定した範囲だけのcpfに変換する機能
- 溶媒が多すぎてBiostation Viewerで開けない際などに使用
- cpfのバージョン変換も可能

#### 実行
    python -m abmptools.convertcpf -i hoge.cpf -f 1-834

    Arguments received:
    Input CPF: hoge.cpf
    Target fragments: 1-834
    Output Version: 23
    start read hoge.cpf
    cpfver 10
    atom section
    fragment section
    dimer distance section
    dipole moment section
    condition and static data section
    monomer section
    dimer section
    hoge.cpf is created.


### 動的IFIE (Dynamical IFIE) 生成 generatedifie
#### 概要
- 指定した複数のcpfから、”平均、標準偏差をひとまとめにしたcpfファイル”を生成
- 3D座標は指定した1構造を代表して表示
- 最新のBiostation Viewerで読み込み可能 (テスト中)

#### サンプルデータ
    abmptools/sample内に2例(TrpCage, CS4)  run.shで実行
    cd sample/
    cd generate_difie/
    cd TrpCage/
    bash run.sh

#### 実行
    python -m abmptools.generate_difie  オプション

    - 入力オプション
        -i 入力cpf名  ※指定必須 (file-xxx-bbb.cpfの形、数字部分をxxxで表記）
        -t 開始番号 終了番号　読込間隔 　※指定必須
        -z ゼロ埋めの桁数 ※デフォルト: "1" (ゼロ埋めなし)
        -s 代表構造の番号: 平均cpfの表示構造番号の指定　※デフォルト: 最初の構造(＝開始番号)
        -f 対象フラグメントの指定:  周囲の水は入れ替わるため、タンパクとリガンド残基番号を手動指定　※デフォルト:  削除なし　半角ハイフンで範囲指定, 1-のみ可能
        -v 出力バージョンの指定：現状はrev23のみ ※デフォルト"23"
        -np 並列数

    - 出力
        - DIFIEcpf本体 ("入力名-DIFIE.cpf")


### 実行例
    `python -m abmptools.generate_difie -i egfr-HYZ_pr_xxx_fmo_mask-renamed-MP2-6-31Gd-nbo.cpf.gz -t 2 6 1 -z 1 -s 0 -f 1-323 -v 23 -np 5`
    `python -m abmptools.generate_difie -i TrpCage-xxx.cpf -t 1 5 1 -z 1 -s 0 -f 1-20 -v 23 -np 5`
    `python -m abmptools.generate_difie -i CS4_ligX_md1_xxxns_mod_mp2_631gd.cpf.gz -t 12 20 2 -z 1 -s 0 -f 1-299 -v 23 -np 5`


### 出力書式
    M- (Mean), S- (STD)の接頭文字でヘッダー出力

    CPF Open1.0 rev23 DIFIE (Generated by ABMPTools 2024-01-16 23:48:59.894215)
           304        20
    M-MUL-HF M-MUL-MP2 M-NPA-HF M-NPA-MP2 M-ESP-HF M-ESP-MP2 S-MUL-HF S-MUL-MP2 S-NPA-HF S-NPA-MP2 S-ESP-HF S-ESP-MP2
    DPM-HF-X DPM-HF-Y DPM-HF-Z DPM-MP2-X DPM-MP2-Y DPM-MP2-Z
    NR HF MP2 MP3
    M-Total M-NR M-HF M-ES M-MP2 M-SCS-MP2(Grimme) M-MP3 M-SCS-MP3(MP2.5) M-HF-BSSE M-MP2-BSSE M-SCS-MP2-BSSE M-MP3-BSSE M-SCS-MP3-BSSE M-EX M-CT M-DQ S-Total S-NR S-HF S-ES S-MP2 S-SCS-MP2(Grimme) S-MP3 S-SCS
    -MP3(MP2.5) S-HF-BSSE S-MP2-BSSE S-SCS-MP2-BSSE S-MP3-BSSE S-SCS-MP3-BSSE S-EX S-CT S-DQ

![DIFIE](img/difie.png)


### pdb2fmo
- pdbにフラグメント情報を割り当てて、FMO実行可能なajf, pdbセットを作成するモジュール
- 中分子分散体やポリマーのpdb構造に、フラグメントを割り当てる際に主に使用。
- 1分子ずつのフラグメント情報を事前に準備する必要有

### udf2fmo
- udfにフラグメント情報を割り当てて、FMO実行可能なajfを作成するモジュール

### ajf2config
- フラグメント情報を読み込んで、pythonの辞書型に保存するモジュール
- 各種フラグメント作成で使用
- usage: `python -m abmptools.ajf2config -i xxx.ajf yyy.ajf`

### ajfserial
- 雛形ajfから連番を作成する
- e.g.) `python -m abmptools.ajfserial -i 1eo8-ff14sb-xxxps-renamed-HF-STO-3G-nbo.ajf -t 100000 200000 5000 -str xxx`

### pdbmodify
PDBの情報を編集するモジュール
残基名の変更、残基番号、原子名のつけなおし、など

e.g.) `python -m abmptools.pdbmodify -mode rename -str ' NA' 'NA ' 'CYX' 'CYS' ' CL' 'CL ' -i \*pdb`
e.g.) `python -m abmptools.pdbmodify -mode resnum -aresname -reatm -i cyc.pdb.1`

    usage: e.g.)
                    python -m abmptools.pdbmodify -addc 307 312 C -i *.pdb
                    python -m abmptools.pdbmodify -move -p 10 10 10 -aresname -reid -reatm -addc 307 312 C -i *pdb
                    python -m abmptools.pdbmodify -move -mol 3 --into -reres -reatm -addc 307 312 C -i *pdb
                    python -m abmptools.pdbmodify -mode rename -str 001 MRT 002 CD7 003 NA 004 WAT -i *pdb

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


### addsolvfrag
    雛形ajfのフラグメント情報に、読み込んだPDBの追加溶媒情報を追加して、新たにajfを作成するモジュール
    手動分割で溶媒の数が異なっている複数構造の調整の際に使用


    usage: e.g.)
                    python -m abmptools.addsolvfrag -temp xxx-temp.ajf -solv HOH WAT NA -i *.pdb
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
