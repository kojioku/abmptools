from UDFManager import *
import sys
import os
import rmdpd.udfrm_io as uio
import numpy as np
import math
import subprocess
import re
import time
import copy
import fmor.rev_md_fmo as fr


if __name__ == "__main__":
    fhead = 'cbzrich-2-for_abmp.pdb'
    tstart = 1
    tend = 4
    interval = 1

    # tgtname = "Cholesterol"
    tgtname = "CBZ"

    tgtmolid = 1
    print('target:', tgtname, "mol", tgtmolid)

    ## -- user setting --
    # read info
    mode = 'resnum' #rfile, resnum
    assignmolname = False
    ## -- setting end --

    # main
    # fname = str(argvs[1])

    cobj = uio.udfrm_io()
    tgtmolid += 1
    totalrec = int((tend - tstart + 1)/interval)
    oname =  os.path.splitext(fhead)[0] + '-' + str(tstart) + '-' + str(tend) + '.log'

    for i in range(tstart, tend + 1, interval):
        fname = fhead + '.' + str(i)
        obj = fr.rmap_fmo()
        obj.getmode = mode
        obj.assignmolname = assignmolname

        print('infile:', fname)

        # get pdbinfo
        totalMol, atomnameMol, molnames, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfo(fname)
        mollist = [i for i in range(totalMol)]

    # get molname list
    # print(totalmol, totalrec)
    # print(molnames)
    print(molnames)
    totalmol = len(molnames)

    # get target mol index
    targets = []
    for i in range(totalmol):
        if molnames[i] == tgtname:
            targets.append(i)
    print('target id', targets)
# 
# 
    # get target posmol for all rec
    poss = []
    mol = targets[tgtmolid]
    print("--- target mol id:", mol, "---")

    for record in range(totalrec):
        # print(posMol[mol])
        poss.append(posMol[mol])

    print(poss)
# 
    # centerOfMol: com of each rec
    cocs = []
    for i in range(len(poss)):
        cocs.append(cobj.getCenter(poss[i]).tolist())
    # print('cocs', cocs)

    # getdist
    dists = []
    for i in range(totalrec-1):
        dists.append(cobj.getdist_list(cocs[i], cocs[i+1]))
    print('dists', dists)

    f = open(oname, 'w')
    rec = tstart
    print('rec', 'dist', file=f)
    for i in range(totalrec-1):
        print (rec, dists[i], file=f)
        rec += interval

