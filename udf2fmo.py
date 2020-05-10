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
# Matrix operation

if __name__ == "__main__":
    # main
    argvs = sys.argv
    fname = str(argvs[1])
    oname = str(argvs[2])
    _udf_ = UDFManager(fname)

    aobj = ampt.setfmo()
    totalMol, totalRec = aobj.gettotalmol_rec(_udf_)
    path = ['.', '.']

    param_read = {}
    exec(open("input_param", 'r').read(), param_read)
    param_rfmo = param_read['param']
    aobj.setrfmoparam(param_rfmo)

    # tgtpos = param_rfmo['tgtpos']
    # criteria = param_rfmo['criteria']
    # molname = param_rfmo['molname']

    aobj.mainpath = '.'
    aobj.getcontact_rmapfmo(
        totalRec-1, _udf_, totalMol, totalMol, path, oname)
