import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools as abmp
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

        aobj = abmp.setfmo()
        path = 'for_abmp'

        param_read = {}
        exec(open("input_param", 'r').read(), param_read)
        param_rfmo = param_read['param']
        aobj.setrfmoparam(param_rfmo)

        if aobj.cutmode == 'sphere':
            oname = oname + '-' + aobj.cutmode + 'p' + str(aobj.tgtpos[0]) + '_' + str(aobj.tgtpos[1]) + '_' + str(aobj.tgtpos[2]) + '_ar' + str(aobj.criteria)

        elif aobj.cutmode == 'around':
            oname = oname + '-' + aobj.cutmode + '_ar' + str(aobj.criteria)

        elif aobj.cutmode == 'cube':
            oname = oname + '-' + aobj.cutmode + 'p' + str(aobj.tgtpos[0]) + '_' + str(aobj.tgtpos[1]) + '_' + str(aobj.tgtpos[2]) + '_x' + str(aobj.criteria[0]) + '_y' + str(aobj.criteria[1]) + '_z' + str(aobj.criteria[2])

        if aobj.cutmode == 'none':
            oname = oname +  '-for_abmp'

        if aobj.cutmode == 'around':
            print('molset', aobj.molname)
            print('solute', aobj.solutes)

        print('piedaflag', aobj.piedaflag)

        if len(aobj.ionname) != 0:
            print('ion mode: ', aobj.ionmode)
            print('ion name: ', aobj.ionname)
        aobj.getcontact_rmapfmopdb(path, fname, oname)
