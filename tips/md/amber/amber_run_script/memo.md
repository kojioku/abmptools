# frcmod取得のcmd
reduce sustiva.pdb > sustiva_h.pdb
antechamber -i sustiva_new.pdb -fi pdb -o sustiva.mol2 -fo mol2 -c bcc -s 2
parmchk2 -i sustiva.mol2 -f mol2 -o sustiva.frcmod

antechamberはqsubでいれないと時間内におわんない

流れとして
パラメータない1分子のpdbを準備
antechamber でamber名称のmol2にする
parmchk2でmol2からfrcmodをつくる


# 大事なこと

まず
1. MOEでbulkモデルpdbをつくる
2. そこから1分子をとりだして、1分子のantechamberを確実に行
3. parmchk2 -i sustiva.mol2 -f mol2 -o sustiva.frcmodでfrcmodをつくる
3. mol2, frcmodを元に、leap_getlib.sh xxx.mol2 xxx.frcmod でlibファイルをつくる
出力するlibの名称は引数にしてないからファイル編集する必要あり
ここでエラーが出る場合は、connectテーブルの誤りの可能性が高いのでopenbabelのコネクト情報を拝借

3. 1のbulkモデルpdbに、2のmol2の原子ラベル情報を入れ込む(原子ラベル番号と、残基名情報)
タンパク側の読み込み問題がある場合はpdb4amberも打つ必要があって少々面倒

4. 3のpdbファイルをleap_pcpdb.sh pdb xxx.lib で実行
loadoffするのがlibファイル
frcmod変えたり、solvationどうするかとかはべた書きで編集する必要あり
