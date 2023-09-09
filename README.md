ABMPTools (ABINIT-MP Tools)
====

Overview
フラグメント分子軌道計算プログラムABINIT-MPの解析支援ツールです。
※ 2022/09現在、意図部を除くサンプルファイルは非公開リポジトリとさせて頂いています。

## Description

- 解析支援(getifiepieda)
    - 一点構造
        - 対象フラグメントから距離を指定したIFIEテーブル出力
        - 対象分子から一定距離以内にある分子とのIFIEテーブル出力
        - 分子内フラグメントIDによる指定(特定のユニット間の相互作用の解析）
        - 対象のフラグメント群-フラグメント群間の相互作用マトリックス出力(SVD用)

    - 古典MDサンプリング構造
        - 対象フラグメント間のIFIEの時系列変化
        - 対象フラグメント-指定した残基種間 IFIE和の時系列出力(水との相互作用を一括計算など)
        - 対象フラグメント-フラグメント群間の時間-IFIEの相互作用マトリックス出力 (1:1, 1:N, N:1, N:Nの時間変化の一括出力）

    - SVD実行スクリプト

- CPF parser(cpfmanager)
    - cpf残基切り出し、バージョンコンバート(convertcpf)
    - 動的平均cpf作成機能(generate_difie)

- FMO実行支援
    - PDB読み込み(pdbread)　-> 雛形ajf生成(generateajf)
    - 非タンパク系のPDB, フラグメント情報の処理（分子集合体向け）
        - 一分子のフラグメント情報の全分子への割り当て: 分子集合体向け
        - 任意のサイズでの分子系の切り出し＋フラグメント割り当て

- 各種ファイル変換
    - 結晶構造(cif)ファイルからFMO用ファイルへの変換
    - OCTA udf ファイルからFMO用ファイルへの変換
    - 古典MD後からFMO用ファイルへの変換

## Requirement

`numpy`
`pandas`

## Usage

 `python -m abmptools.addsolvfrag`
 `python -m abmptools.generateajf`
 `python -m abmptools.getifiepieda`
 `python -m abmptools.pdb2fmo`
 `python -m abmptools.pdbmodify`
 `python -m abmptools.pdbread`
 `python -m abmptools.setupserialpdb`
 `python -m abmptools.udf2fmo`

## Install
`python setup.py install`

## Contribution

## Licence

## Author

[Koji Okuwaki](koujioku81@gmail.com)

## example

