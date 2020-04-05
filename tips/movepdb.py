import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import fmor.rev_md_fmo as fr

# Koji Okuwaki: Update 2020/03/22
# moveintocell and assignmolname

if __name__ == "__main__":
    ## -- user setting --
    # read info
    mode = 'resnum' #rfile, resnum
    assignmolname = False
    refreshatmtype = False

    # move info
    moveflag = False
    movemode = 'mol' # pos or mol

    # --- mol mode
    tgtmol = 2
    # --- pos mode
    tgtpos = [10.0, 10.0, 10.0]
    intoflag = False
    ## -- setting end --

    # main
    argvs = sys.argv
    # fname = str(argvs[1])

    for arg in argvs:
        if arg == '--move':
            moveflag == True
        if arg == '--nomove':
            moveflag == False

    for i in range(len(argvs)):
        if i == 0:
            continue
        if argvs[i][0:2] == '--':
            continue
        fname = argvs[i]
        oname, ext = os.path.splitext(fname)
        if ext != '.pdb':
            oname = oname.split('.pdb')[0] + ext.split('.')[1] + '-moved'
        else:
            oname = oname + '-moved'

        obj = fr.rmap_fmo()
        obj.getmode = mode
        obj.assignmolname = assignmolname
        obj.refreshatmtype = refreshatmtype

        print('infile:', fname)
        print('oname:', oname)
        print('centered-molid:', tgtmol - 1)

        # get pdbinfo
        totalMol, atomnameMol, molnames, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(fname)
        mollist = [i for i in range(totalMol)]
        cellsize = obj.getpdbcell(fname)
        obj.cellsize = cellsize

        print('totalMol:', totalMol)

        if len(obj.cellsize) == 0:
            obj.cellsize = 0
            print('cellinfo: None')
        else:
            print('cellsize:', obj.cellsize)


        if moveflag == True:
            # get center of solute
            if movemode == 'mol':
                coctgt = obj.getCenter(posMol[tgtmol-1])
            elif movemode == 'pos':
                coctgt = tgtpos
            transVec = np.array(-coctgt, dtype='float64')
            # print(transVec)

            # move
            posmoveMol = []
            for i in range(totalMol):
                posmove = obj.movemoltranspdb(posMol[i], transVec)
                posmoveMol.append(posmove)

            if intoflag == True:
                posintoMol = obj.moveintocellpdb(posmoveMol, totalMol, cellsize)

            else:
                posintoMol = copy.deepcopy(posmoveMol)

        else:
            if intoflag == True:
                posintoMol = obj.moveintocellpdb(posMol, totalMol, cellsize)

            else:
                posintoMol = copy.deepcopy(posMol)

        # write
        obj.exportardpdbfull(oname, mollist, posintoMol, atomnameMol, molnames, heads, labs, chains, resnums, codes, occs, temps, amarks, charges)


