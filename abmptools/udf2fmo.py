"""UDFファイル（COGNAC形式）とsegment_data設定からABINIT-MP入力ファイルを生成するCLIツール。"""
import ast
import numpy as np
from UDFManager import *
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools as ampt
import argparse

def get_args():
    """コマンドライン引数を解析する。"""
    parser = argparse.ArgumentParser(
                prog='udf2fmo.py', # program name
                usage='python udf2fmo.py -i xxx.udf -o yyy', # program usage
                description='generate ABNITMP input (ajf,pdb set) from udf(cognac) and segment_data file',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--incoord',
                        help='coordinate file (udf)',
                        # nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-p', '--parameter',
                        help='parameter file',
                        default='input_param')

    parser.add_argument('-o', '--output',
                        help='output file',
                        default=None)

    parser.add_argument('-s', '--solutes',
                        help='solute id',
                        nargs = '*',
                        type=int,
                        default=None,)

    parser.add_argument('-r', '--record',
                        help='record id',
                        type=int,
                        default=None,)

    # get args
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # main
    args = get_args()

    print('coord(udf) =', args.incoord)
    print('parameter = ', args.parameter)

    is_solutes = False
    if args.solutes is not None:
        is_solutes = True

    _udf_ = UDFManager(args.incoord)

    aobj = ampt.setfmo()
    totalMol, totalRec = aobj.gettotalmol_rec(_udf_)
    path = ['.', '.']

    if args.record is None:
        tgtrec = totalRec-1
    else:
        tgtrec = args.record

    # print(totalRec)
    print('tgtrec:',tgtrec)
    # sys.exit()

    param_read = {}
    with open(args.parameter, 'r') as _f:
        _tree = ast.parse(_f.read())
    for _node in _tree.body:
        if isinstance(_node, ast.Assign) and len(_node.targets) == 1 and isinstance(_node.targets[0], ast.Name):
            param_read[_node.targets[0].id] = ast.literal_eval(_node.value)
    param_rfmo = param_read['param']
    aobj.setrfmoparam(param_rfmo)

    if is_solutes:
        aobj.solutes = args.solutes

    print('solutes', aobj.solutes)

    if args.output is None:
        oname= os.path.splitext(args.incoord)[0].split('/')[-1] + '-' + aobj.cutmode + '-' + 'rec' + str(tgtrec)

        if aobj.cutmode == 'around':
            oname += '-solu' + str(aobj.solutes[0]) + '-' + str(aobj.solutes[-1])
    else:
        oname = args.output
    print('out(pdb, ajf) = ', oname)

    # sys.exit()

    aobj.mainpath = '.'
    aobj.getcontact_rmapfmo(
        tgtrec, _udf_, totalMol, totalMol, path, oname)
