"""テンプレートAJFファイル中の文字列を連番に置換し、複数のAJFファイルを一括生成するCLIツール。"""
import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools as ampt
import argparse

# Koji Okuwaki: Update 2020/03/22
# moveintocell and assignmolname

def get_args():
    """コマンドライン引数を解析する。"""
    parser = argparse.ArgumentParser(
                prog='pdbmodify.py', # program name
                usage='''e.g.)
                python ajfserial.py -i template.ajf -t 1 100 1 -str xxx
                ''', # program usage
                description='description',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--input',
                        help='input template ajf',
                        # action='append',
                        required=True)

    parser.add_argument('-t', '--time',
                        help='time info',
                        type=int,
                        nargs=3,
                        )

    parser.add_argument('-str', '--string',
                        help='rename data',
                        default='xxx',
                        )

    args = parser.parse_args()
    return args


def main():
    args = get_args()
    print('input =', args.input)
    print('time =', args.time)
    print('replaceword =', args.string)

    ## -- user setting --
    # read info

    infile = args.input
    lines = open(infile, 'r').readlines()
    replaceword = args.string
    fname = infile.split(replaceword)
    start, end, intv = args.time

    for i in range(start, end+1, intv):
        oname = fname[0] + str(i) + fname[1]
        with open(oname, 'w') as outf:
            for line in lines:
                line = line.replace(replaceword, str(i))
                print(line[:-1], file=outf)
            print('finished', oname)


if __name__ == "__main__":
    main()

