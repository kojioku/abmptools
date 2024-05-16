import abmptools
import pandas as pd
import argparse
import numpy as np


def getargs():

    parser = argparse.ArgumentParser(description='引数の例')

    # 入力ファイル名の引数を追加
    parser.add_argument('-i', '--input', type=str,
                        required=True, help='入力ファイル名')

    # 出力ファイル名の引数を追加
    parser.add_argument('-o', '--output', type=str,
                        default='', help='出力ファイル名')

    # フラグメント番号範囲の引数を追加（2つの引数を受け取る）
    parser.add_argument('-f', '--fragment', type=int,
                        nargs=2, required=True, help='フラグメント番号範囲 (2つの引数)')

    # オプションのandflagの引数を追加
    parser.add_argument('--and', dest='andflag',
                        action='store_true', help='andflagをonにする')

    # 引数を解析
    args = parser.parse_args()

    # 引数を表示（デバッグ用）
    print(f"入力ファイル名: {args.input}")
    print(f"出力ファイル名: {args.output}")
    print(f"フラグメント番号範囲: {args.fragment}")
    print(f"andflag: {args.andflag}")

    return args


if __name__ == '__main__':
    args = getargs()
    cpf = abmptools.CPFManager()
    cpf.parse(args.input)
    print("Read CPF file finished.")

    # print(cpf.diminfo.head())
    # filter section
    print("Filtering...")
    if args.andflag:
        filtered = cpf.diminfo[
            (cpf.diminfo['fragi'].between(args.fragment[0], args.fragment[1])) &
            (cpf.diminfo['fragj'].between(args.fragment[0], args.fragment[1]))]
    else:
        filtered = cpf.diminfo[
            (cpf.diminfo['fragi'].between(args.fragment[0], args.fragment[1])) |
            (cpf.diminfo['fragj'].between(args.fragment[0], args.fragment[1]))]

    filtered_copy = filtered.copy()
    filtered_copy.loc[filtered_copy['min-dist'] == 0.0,
                      filtered_copy.select_dtypes(include=[np.float64]).columns] = 0.0

    filtered_copy = filtered_copy.loc[:, (filtered_copy != 0.0).any(axis=0)]
    filtered_copy['MP2-Total'] = filtered_copy['NR'].fillna(0) + filtered_copy['HF'].fillna(0) \
        + filtered_copy['MP2'].fillna(0)
    filtered_copy = filtered_copy.drop(columns=['NR', 'HF'])
    filtered_copy = filtered_copy.rename(columns={'MP2': 'DI'})

    cols_to_exclude = ['min-dist', 'DQ']

    # float型の列を選択し、除外する列を除く
    cols_to_multiply = filtered_copy.select_dtypes(
        include=[float]).columns.difference(cols_to_exclude)

    # 選択した列を627.5095倍する
    filtered_copy[cols_to_multiply] = filtered_copy[cols_to_multiply] * 627.5095
    filtered_copy['min-dist'] = filtered_copy['min-dist'] * 0.529177

    print("Writing output file...")
    if args.output == '':
        output = args.input.replace('.cpf', '_filtered.csv')
    else:
        output = args.output
    filtered_copy.to_csv(output, float_format='%.6f')

    print("Write output file finished: " + output)
