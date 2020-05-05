import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import ampt.cutset_fmo as cutf

# Koji Okuwaki: Update 2020/03/22
# moveintocell and assignmolname

if __name__ == "__main__":
    ## -- user setting --
    # read info
    mode = 'resnum' #rfile, resnum
    assignmolname = False
    refreshatmtype = False
    refreshresid = False

    # move info
    moveflag = False
    movemode = 'mol' # pos or mol

    addchain = True
    addres_start = 308
    addres_end = 312
    chainlab = 'C'


    addres = [i for i in range(addres_start, addres_end+1)]

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
            oname = oname.split('.pdb')[0] + ext.split('.')[1] + '-mod'
        else:
            oname = oname + '-mod'

        obj = cutf.cutset_fmo()
        obj.getmode = mode
        obj.assignmolname = assignmolname
        obj.refreshatmtype = refreshatmtype
        obj.refresresid = refreshresid

        print('infile:', fname)
        print('oname:', oname)
        print('centered-molid:', tgtmol - 1)

        # get pdbinfo
        totalMol, atomnameMol, molnames, nlabmols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(fname)

        mollist = [i for i in range(totalMol)]
        cellsize = obj.getpdbcell(fname)
        obj.cellsize = cellsize

        print('totalMol:', totalMol)

        # print(resnums)
        print(addres)
        alreadys = []
        for i in range(len(resnums)):
            # print(resnums[i][0])
            tgt = int(resnums[i][0])
            if tgt in alreadys:
                continue
            if tgt in addres:
                print(resnums[i])
                # print(nlabmols[i])
                for j in range(len(chains[i])):
                    chains[i][j] = chainlab
                alreadys.append(tgt)

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


