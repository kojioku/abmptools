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
    ## -- user setting --
    tgtmol = 2
    moveflag = False
    ## -- setting end --

    # main
    argvs = sys.argv
    fname = str(argvs[1])
    oname, ext = os.path.splitext(fname)
    oname = oname.split()[-1] + '-moved'

    obj = fr.rmap_fmo()

    print('infile:', fname)
    print('centered-molid:', tgtmol - 1)

    # get pdbinfo
    totalMol, atomnameMol, molnames, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(fname)
    mollist = [i for i in range(totalMol)]
    cellsize = obj.getpdbcell(fname)

    if moveflag == True:
        # get center of solute
        coctgt = obj.getCenter(posMol[tgtmol-1])
        transVec = np.array(-coctgt, dtype='float64')
        # print(transVec)

        # move
        posmoveMol = []
        for i in range(totalMol):
            posmove = obj.movemoltranspdb(posMol[i], transVec)
            posmoveMol.append(posmove)

        #
        posintoMol = obj.moveintocellpdb(posmoveMol, totalMol, cellsize)

    else:
        posintoMol = obj.moveintocellpdb(posMol, totalMol, cellsize)

    # write
    obj.amarkflag = True
    obj.exportardpdbfull(oname, mollist, posintoMol, atomnameMol, molnames, heads, labs, chains, resnums, codes, occs, temps, amarks, charges)


