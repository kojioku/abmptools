import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import fmor.rev_md_fmo as fr
# Matrix operation


if __name__ == "__main__":
    # main
    argvs = sys.argv
    fname = str(argvs[1])
    oname, ext = os.path.splitext(fname)
    oname = oname.split()[-1] + '-for_abmp'

    obj = fr.rmap_fmo()
    path = ['.', '.']

    param_read = {}
    exec(open("input_param", 'r').read(), param_read)
    param_rfmo = param_read['param']
    obj.setrfmoparam(param_rfmo)

    tgtpos = param_rfmo['tgtpos']
    criteria = param_rfmo['criteria']
    molname = param_rfmo['molname']
    print(molname)
    obj.getcontact_rmapfmopdb(
        path, molname, fname, oname, tgtpos, criteria)
