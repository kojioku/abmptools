import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import abmptools as ampt
import collections
# Matrix operation


if __name__ == "__main__":
    ## -- user setting --
    # read info
    mode = 'resnum' #rfile, resnum
    assignmolname = False
    refreshatmtype = False
    ## -- setting end --

    aobj = ampt.setfmo()
    aobj.getmode = mode
    aobj.assignmolname = assignmolname
    aobj.refreshatmtype = refreshatmtype

    # main
    argvs = sys.argv
    # fname = str(argvs[1])

    pdbname = sys.argv[1]
    ajfname = sys.argv[2]

    print('infile:', pdbname)

    # get pdbinfo
    # totalMol, atomnameMol, molnames, anummols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfowrap(pdbname)
    aobj = aobj.getpdbinfowrap(pdbname)

    totalMol = aobj.totalRes
    atomnameMol = aobj.atmtypeRes
    molnames = aobj.resnames
    anummols = aobj.gatmlabRes
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

    # get ajfinfo
    aobj = aobj.readajf(ajfname)

    fatomnums = aobj.fatomnums
    fchgs = aobj.fchgs
    fbaas = aobj.fbaas
    fatminfos = aobj.fatminfos
    connects = aobj.connects

    # print(fatomnums, fchgs, fbaas, fatminfos, connects)

    # -- devide asn and nag --
    print('# === check nag === ')
    nagatoms = []
    nagbdas = []
    doubles = []
    nagmolids = []
    for i in range(len(molnames)):
        if molnames[i] == 'NAG':
            nagatoms.append(list(map(int, anummols[i])))
            nagmolids.append(i)
            # print(i, anummols[i][0])
            # print(molnames[i], anummols[i])
            print(atomnameMol[i])
            bdaidx = atomnameMol[i].index(' C1 ')
            nagbdas.append(int(anummols[i][bdaidx]))
            if len(anummols[i]) == 27:
                print('mol', i, 'and', i+1, 'are connected')
                doubles.append([i, i+1])
    print(nagatoms)
    print('nagbdas', nagbdas)
    print('doubles', doubles)

    # -- check CYS bridge --
    print('# === check cys bridge === ')
    cysatoms = []
    cysbdas = []
    cysmolid = []
    for i in range(len(molnames)):
        if molnames[i] == 'CYS' or molnames[i] == 'CYX':
            cysatoms.append(list(map(int, anummols[i])))
            cysmolid.append(i)

    bridgeds = []
    for i in range(len(cysatoms)):
        for j in range(len(fatminfos)):
            if cysmolid[i] == j:
                continue
#             if molnames[j] == 'CYX' or molnames[j] == 'CYS':
#                 print(j, 'is', molnames[j])
            if cysatoms[i][0] in fatminfos[j]: # search (asn + nag)
                print('CYS', cysmolid[i], 'info is in frag', j, '(bridged)')
                bridgeds.append(cysmolid[i])
                break

    fatomnums, fchgs, fbaas, connects, fatminfos = aobj.modifyfragparam(totalMol, atomnameMol, molnames, anummols, posMol, heads, labs, chains
                                                                   ,resnums ,codes ,occs ,temps ,amarks ,charges, fatomnums, fchgs, fbaas, fatminfos, connects, bridgeds, doubles, nagatoms, nagmolids, nagbdas)

    aobj.ajf_method = "MP2"
    aobj.ajf_basis_set = "6-31G*"
    aobj.abinit_ver = 'rev10'
    aobj.abinitmp_path = 'abinitmp'
    aobj.pbmolrad = 'vdw'
    aobj.pbcnv = 1.0
    aobj.piedaflag = True
    aobj.readgeom = pdbname

    aobj.saveajf(os.path.splitext(pdbname)[0] + '-nagfragmod.ajf')

