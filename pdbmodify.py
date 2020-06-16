import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools as ampt

# Koji Okuwaki: Update 2020/03/22
# moveintocell and assignmolname

if __name__ == "__main__":
    ## -- user setting --
    # read info
    mode = 'resnum' #rfile, resnum
    assignresname = False
    refreshatmtype = False
    refreshresid = False

    # move info
    moveflag = False
    movemode = 'mol' # pos or mol

    addchain = True
    addres_start = 307
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

        aobj = ampt.setfmo()
        aobj.getmode = mode
        aobj.assignresname = assignresname
        aobj.refreshatmtype = refreshatmtype
        aobj.refreshresid = refreshresid

        print('infile:', fname)
        print('oname:', oname)
        print('centered-molid:', tgtmol - 1)

        # get pdbinfo
        # totalMol, atomnameMol, molnames, nlabmols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(fname)

        aobj = aobj.getpdbinfowrap(fname)
        totalMol = aobj.totalRes
        atomnameMol = aobj.atmtypeRes
        molnames = aobj.resnames
        nlabmols = aobj.gatmlabRes
        posMol = aobj.posRes
        heads = aobj.headRes
        labs = aobj.labRes
        chains = aobj.chainRes
        resnums = aobj.resnumRes
        codes = aobj.codeRes
        occs = aobj.occRes
        temps = aobj.tempRes
        amarks = aobj.amarkRes
        charges = aobj.chargeRes

        mollist = [i for i in range(totalMol)]
#         cellsize = aobj.getpdbcell(fname)
#         aobj.cellsize = cellsize

        print('totalMol:', totalMol)

        # print(resnums)
        if addchain == True:
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
                        aobj.chainRes[i][j] = chainlab
                    alreadys.append(tgt)

#         if len(aobj.cellsize) == 0:
#             aobj.cellsize = 0
#             print('cellinfo: None')
#         else:
#             print('cellsize:', aobj.cellsize)


        if moveflag == True:
            # get center of solute
            if movemode == 'mol':
                coctgt = aobj.getCenter(posMol[tgtmol-1])
            elif movemode == 'pos':
                coctgt = tgtpos
            transVec = np.array(-coctgt, dtype='float64')
            # print(transVec)

            # move
            posmoveMol = []
            for i in range(totalMol):
                posmove = aobj.movemoltranspdb(posMol[i], transVec)
                posmoveMol.append(posmove)

            if intoflag == True:
                posintoMol = aobj.moveintocellpdb(posmoveMol, totalMol, cellsize)

            else:
                posintoMol = copy.deepcopy(posmoveMol)

        else:
            if intoflag == True:
                posintoMol = aobj.moveintocellpdb(posMol, totalMol, cellsize)

            else:
                posintoMol = copy.deepcopy(posMol)

        # write
        aobj.posRes = posintoMol
        aobj.exportardpdbfull(oname, mollist)



