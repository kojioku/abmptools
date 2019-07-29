import numpy as np
from UDFManager import *
import sys
import os
import math
import subprocess
import re
import time
import copy
import fcewsmb.udf_io as mbu
import fcews.abinit_io as fab
import fcewsmb.udfcreate as ufc
import rmdpd.udfrm_io as rud
import fmor.rev_md_fmo as frev
# Matrix operation

if __name__ == "__main__":
    # main
    argvs = sys.argv
    fname = str(argvs[1])
    oname = str(argvs[2])
    _udf_ = UDFManager(fname)

    obj = frev.rmap_fmo()
    totalMol, totalRec = obj.gettotalmol_rec(_udf_)
    path = ['.', '.']

    param_read = {}
    exec(open("input_param", 'r').read(), param_read)
    param_rfmo = param_read['param']
    obj.setrfmoparam(param_rfmo)

    # tgtpos = param_rfmo['tgtpos']
    # criteria = param_rfmo['criteria']
    # molname = param_rfmo['molname']

    obj.getcontact_rmapfmo(
        totalRec-1, _udf_, totalMol, totalMol, path, oname)
