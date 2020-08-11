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

if __name__ == "__main__":
    # main
    # argvs = sys.argv
    # fname = str(argvs[1])
    # oname = str(argvs[2])

    # create parser
    parser = argparse.ArgumentParser(
                prog='generate ABNITMP input (ajf,pdb set) from udf(cognac) and segment_data file', # program name
                usage='Demonstration of argparser', # program usage
                description='description',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-c', '--coord',
                        help='coordinate file (udf)',
                        # nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-p', '--parameter',
                        help='parameter file',
                        default='input_param')

    parser.add_argument('-o', '--output',
                        help='output file',
                        default='out-forfmo')

    # get args
    args = parser.parse_args()

    print('coord(udf) =', args.coord)
    print('parameter = ', args.parameter)
    print('out(pdb, ajf) = ', args.output)

    _udf_ = UDFManager(args.coord)

    aobj = ampt.setfmo()
    totalMol, totalRec = aobj.gettotalmol_rec(_udf_)
    path = ['.', '.']

    param_read = {}
    exec(open(args.parameter, 'r').read(), param_read)
    param_rfmo = param_read['param']
    aobj.setrfmoparam(param_rfmo)

    aobj.mainpath = '.'
    aobj.getcontact_rmapfmo(
        totalRec-1, _udf_, totalMol, totalMol, path, args.output)
