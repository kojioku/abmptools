import numpy as np
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
    # fname = argvs

    for i in range(len(argvs)):
        if i == 0:
            continue
        fname = argvs[i]
        oname, ext = os.path.splitext(fname)

        obj = ampt.setfmo()
        path = ['.', 'for_abmp']

#         param_read = {}
#         exec(open("input_param", 'r').read(), param_read)
#         param_rfmo = param_read['param']
#         obj.setrfmoparam(param_rfmo)

        # if obj.cutmode == 'sphere' or obj.cutmode == 'around':
        #     oname = oname + '-' + obj.cutmode + '-' + str(obj.criteria) + '-for_abmp'

        # elif obj.cutmode == 'cube':
        #     oname = oname + '-' + obj.cutmode + '-' + str(obj.criteria[0]) + '-' + str(obj.criteria[1]) + '-' + str(obj.criteria[2]) + '-for_abmp'

        # if obj.cutmode == 'none':
        #     oname = oname +  '-for_abmp'

        # if obj.cutmode == 'around':
        #     print('molset', obj.molname)
        #     print('solute', obj.solutes)

        # print('piedaflag', obj.piedaflag)

        # if len(obj.ionname) != 0:
        #     print('ion mode: ', obj.ionmode)
        #     print('ion name: ', obj.ionname)
        # obj.getcontact_rmapfmopdb(path, fname, oname)
        # totalMol, atomnameMol_orig, molnamelist_orig, posMol_orig, heads_orig, labs_orig, chains_orig ,resnums_orig ,codes_orig ,occs_orig ,temps_orig ,amarks_orig ,charges_orig = obj.getpdbinfo(fname)
        obj.readpdb(fname)
        # print(obj.posRes)
        posMol_orig = obj.posRes
        xs = []
        ys = []
        zs = []
        for poss in posMol_orig:
            for pos in poss:
                xs.append(pos[0])
                ys.append(pos[1])
                zs.append(pos[2])
        print('atomnum', len(xs), len(ys), len(zs))
        print('x max',  max(xs), 'min', min(xs))
        print('y max', max(ys), 'min', min(ys))
        print('z max', max(zs), 'min', min(zs))
        print('center =', (max(xs) + min(xs))/2.0, (max(ys) + min(ys))/2.0, (max(zs) + min(zs))/2.0)
