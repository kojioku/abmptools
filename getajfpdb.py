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
        path = 'for_abmp'

        param_read = {}
        exec(open("input_param", 'r').read(), param_read)
        param_rfmo = param_read['param']
        obj.setrfmoparam(param_rfmo)

        if obj.cutmode == 'sphere' or obj.cutmode == 'around':
            oname = oname + '-' + obj.cutmode + 'p' + str(obj.tgtpos[0]) + '_' + str(obj.tgtpos[1]) + '_' + str(obj.tgtpos[2]) + '_ar' + str(obj.criteria)

        elif obj.cutmode == 'cube':
            oname = oname + '-' + obj.cutmode + 'p' + str(obj.tgtpos[0]) + '_' + str(obj.tgtpos[1]) + '_' + str(obj.tgtpos[2]) + '_x' + str(obj.criteria[0]) + '_y' + str(obj.criteria[1]) + '_z' + str(obj.criteria[2])

        if obj.cutmode == 'none':
            oname = oname +  '-for_abmp'

        if obj.cutmode == 'around':
            print('molset', obj.molname)
            print('solute', obj.solutes)

        print('piedaflag', obj.piedaflag)

        if len(obj.ionname) != 0:
            print('ion mode: ', obj.ionmode)
            print('ion name: ', obj.ionname)
        obj.getcontact_rmapfmopdb(path, fname, oname)
