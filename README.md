ABMPTools (ABINIT-MP Tools)
====

Overview
abinitmpの解析支援ツール

## Description

- 解析支援
    - 一点構造
        - 対象フラグメントから距離を指定したifieテーブル出力
        - 対象分子から一定距離以内にある分子とのifieテーブル出力
        - 分子内フラグメントIDによる指定(特定のユニット間の相互作用の解析）
        - 対象のフラグメント群-フラグメント群間の相互作用マトリックス出力(svd用)

    - MDサンプリング構造】
        - 対象フラグメント間のIFIEの時系列変化
        - 対象フラグメント-指定した残基種間 IFIE和の時系列出力(水との相互作用を一括計算など)
        - 対象フラグメント-フラグメント群間の時間-IFIEの相互作用マトリックス出力 (1:1, 1:N, N:1, N:Nの時間変化の一括出力）

    - SVD実行スクリプト


- FMO実行支援
    - PDB読み込み　-> 雛形ajf生成
    - 非タンパク系のPDB, フラグメント情報の処理（分子集合体向け）
        - 一分子のフラグメント情報の全分子への割り当て: 分子集合体向け
        - 任意のサイズでの分子系の切り出し＋フラグメント割り当て

## Requirement


## Usage


## Install
python setup.py install

## Contribution


## Licence


## Author

[Koji Okuwaki](okuwaki@rikkyo.ac.jp)

## example


## update
v1.0.1 create repository
v1.2.0 reform directory
v1.4.0 add svd and get matrix for mp3
v1.4.1 add single N:N mode

