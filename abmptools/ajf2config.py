"""ABINIT-MPフラグメント情報ファイル（AJF）からsegment_data設定ファイルを生成するCLIツール。"""
import sys
import math
import os
import abmptools as ampt
import argparse

def get_args():
    """コマンドライン引数を解析する。"""
    parser = argparse.ArgumentParser(
                prog='ajf2config.py', # program name
                usage='python ajf2config.py -i xxx.frag yyy.frag', # program usage
                description='generate frag config py file from ABNITMP fragment info (segment_data) file',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--input',
                        help='input frag info',
                        nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-o', '--output',
                        help='output config file',
                        default='segment_data.dat')

    # get args
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    print('input =', args.input)
    print('output =', args.output)


    obj = ampt.abinit_io()
    obj.getfragdict(args.input, args.output)





