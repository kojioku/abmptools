import numpy as np
from UDFManager import *
import sys
import os
import math
import subprocess
import re
from multiprocessing import Pool
import rmdpd.udfrm_io as udio

if __name__ == "__main__":
    # main
    argvs = sys.argv
    argc = len(argvs)
    if (argc == 1):
        print(" python udf_to_pdf.py filename")
        quit()

    calc_args = []
    for ifile in range(1, argc):
        filename = argvs[ifile]
        calc_args.append(filename)
        # main(filename)

    fff = udio.udfrm_io()

    p = Pool()
    p.map(fff.run_convert, calc_args)
