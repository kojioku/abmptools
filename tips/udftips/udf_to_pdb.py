import numpy as np
from UDFManager import *
import sys
import os
import math
import subprocess
import re
from multiprocessing import Pool
import abmptools as ampt

if __name__ == "__main__":
    # main
    argvs = sys.argv
    argc = len(argvs)
    if (argc == 1):
        print(" python udf_to_pdf.py filename")
        quit()

    fff = ampt.udfrm_io()

    ## get tgtrec and move
    tgtflag = False
    moveflag = True
    for ival in range(1, argc):
        if argvs[ival] == "--rec":
            tgtrec = argvs[ival + 1]
            tgtflag = True
        if argvs[ival] == "--nomove":
            moveflag = False
        if argvs[ival] == "--mol":
            tgtmol = argvs[ival + 1]
            fff.molflag = True

    if tgtflag == False:
        tgtrec = -1

    if fff.molflag == False:
        tgtmol = -1

    recskip = False
    molskip = False
    calc_args = []
    for ival in range(1, argc):
        filename = argvs[ival]
        if recskip == True:
            recskip = False
            continue
        if molskip == True:
            molskip = False
            continue
        if argvs[ival] == "--rec":
            recskip = True
            continue
        if argvs[ival] == "--mol":
            molskip = True
            continue
        if argvs[ival] == "--nomove":
            continue
        calc_args.append([filename, int(tgtrec), int(tgtmol), moveflag])
        # main(filename)

    print(calc_args)

    p = Pool()
    p.map(fff.run_convert, calc_args)
