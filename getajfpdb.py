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
    # fname = argvs

    for i in range(len(argvs)):
        if i == 0:
            continue
        fname = argvs[i]
        oname, ext = os.path.splitext(fname)

        obj = fr.rmap_fmo()
        path = ['.', 'for_abmp']

        param_read = {}
        exec(open("input_param", 'r').read(), param_read)
        param_rfmo = param_read['param']
        obj.setrfmoparam(param_rfmo)

        if obj.cutmode == 'sphere' or obj.cutmode == 'around':
            oname = oname + '-' + obj.cutmode + '-' + str(obj.criteria) + '-for_abmp'

        elif obj.cutmode == 'cube':
            oname = oname + '-' + obj.cutmode + '-' + str(obj.criteria[0]) + '-' + str(obj.criteria[1]) + '-' + str(obj.criteria[2]) + '-for_abmp'

        if obj.cutmode == 'none':
            oname = oname +  '-for_abmp'

        print('molset', obj.molname)
        print('solute', obj.solutes)
        if len(obj.ionname) != 0:
            print('ion mode: ', obj.ionmode)
            print('ion name: ', obj.ionname)
        obj.getcontact_rmapfmopdb(path, fname, oname)
