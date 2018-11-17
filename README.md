FMORmap (FMO - Reversemap)
====

Overview

[FMORmap](http://www.cenav.org/fcews_ver1_rev2/)
COGNACのMD結果をFMO計算にかけるシステム
## Description


## Requirement

- numpy
- openbabel
- [FCEWS](http://www.cenav.org/fcews_ver1_rev2/)
- [FCEWSMB](http://www.cenav.org/fcews_ver1_rev2/)
- [COGNAC](http://www.cenav.org/fcews_ver1_rev2/)


## Usage

* segmentdata.dat
    - 特殊な分割がある場合はここで指定
* pdbファイル(繰り返し構造を定義した場合、繰り返し一単位をconnectつきのpdbで配置する)
* xyzファイル

python ~/usr_program/py/udf_to_pdb.py で構造を出しておくとよい
python rev_md_fmo.py udf(bdf)名 出力ファイルhead


詳細はdocs/User_Manual-ja.pdf を参照ください。

## Install
python setup.py install

## Contribution

## Licence

## Author

[Koji Okuwaki](okuwaki@rikkyo.ac.jp)

## main 内ファイル構成

## example
python set_monomer -all
python getreverse.py dpd_xxx.bdf 7.1


