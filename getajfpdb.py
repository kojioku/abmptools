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

    obj = fr.rmap_fmo()
    path = ['.', '.']

    param_read = {}
    exec(open("input_param", 'r').read(), param_read)
    param_rfmo = param_read['param']
    obj.setrfmoparam(param_rfmo)

    if obj.cutmode != 'cube':
        oname = oname.split()[-1] + '-' + obj.cutmode + '-' + str(obj.criteria) + '-for_abmp'

    else:
        oname = oname.split()[-1] + '-' + obj.cutmode + '-' + str(obj.criteria[0]) + '-' + str(obj.criteria[1]) + '-' + str(obj.criteria[2]) + '-for_abmp'

    print(obj.molname)
    print(obj.solutes)
    obj.getcontact_rmapfmopdb(path, fname, oname)
