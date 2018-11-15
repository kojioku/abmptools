## 必要なファイル
* segmentdata.dat
    - 特殊な分割がある場合はここで指定
* pdbファイル(繰り返し構造を定義した場合、繰り返し一単位をconnectつきのpdbで配置する)
* xyzファイル

python ~/usr_program/py/udf_to_pdb.py で構造を出しておくとよい
python rev_md_fmo.py udf(bdf)名 出力ファイルhead
