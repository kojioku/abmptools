from UDFManager import *
import sys
import math
# import fcewsmb.udf_io as uio
import abmptools
import numpy as np

if __name__ == "__main__":
    fname = sys.argv[1]
    mode = 'particle' # mol, particle
    tgtname = "Cholesterol"
    tgtpart = "cho1"

    centerpos = [5.0, 5.0, 5.0]
    record = 101

    step = 1.0
    # mol = 0
    # print (fname, record, mol)
    print('target:', tgtname)
    uobj = UDFManager(fname)
    cobj = abmptools.udfrm_io()

    if mode == 'mol':
        ofile= 'dist-' + tgtname + 'from' + str(centerpos[0]) + '-' + str(centerpos[1]) + '-' + str(centerpos[2]) + 'rec' + str(record) + '.txt'
        o2file= 'hist-' + tgtname + 'from' + str(centerpos[0]) + '-' + str(centerpos[1]) + '-' + str(centerpos[2]) + 'rec' + str(record) + '.txt'

    if mode == 'particle':
        ofile= 'dist-' + tgtname + '-' + tgtpart + 'from' + str(centerpos[0]) + '-' + str(centerpos[1]) + '-' + str(centerpos[2]) + 'rec' + str(record) + '.txt'
        o2file= 'hist-' + tgtname + '-' + tgtpart + 'from' + str(centerpos[0]) + '-' + str(centerpos[1]) + '-' + str(centerpos[2]) + 'rec' + str(record) + '.txt'


    totalmol, totalrec, totalatm, cell = cobj.getudfinfowrap(uobj)

    uobj.jump(record)
    cobj.moveintocell_rec(uobj, record, totalmol)

    # get molname list
    molnames = []
    for i in range(totalmol):
        # if i % 10 == 0:
        #     print("read mol name", i)
        molnames.append(cobj.getmolname(i, uobj))
    # print('totalmol', totalmol)
    print('totalrec', totalrec)
    # print(molnames)

    # get target mol index
    targets = []
    for i in range(totalmol):
        if molnames[i] == tgtname:
            targets.append(i)
    print('target id', targets)

    # get target posmol for all rec
    # poss = []
    # mol = targets[tgtmolid]
    # print("--- target mol id:", mol, "---")

    poss = []
    for tgt in targets:
        poss.append(cobj.getposmol(uobj, tgt))
    # print('poss0', poss[0])
    print('num of tgt', len(poss))

    if mode == 'mol':
#     # centerOfMol: com of each molecules
        cocs = []
        for i in range(len(poss)):
            cocs.append(cobj.getCenter(poss[i]).tolist())
        # print('cocs', cocs)
    #
        # getdist
        dists = []
        for i in range(len(cocs)):
            dists.append(cobj.getdist_list(cocs[i], centerpos))
        print('dists', dists)


    if mode == 'particle':
        atoms = []
        for tgt in targets:
            atoms.append(cobj.getnameAtom(uobj, tgt))
        # print(atoms)

        tgtposs = []
        for i in range(len(poss)):
            for j in range(len(atoms[i])):
                if atoms[i][j] == tgtpart:
                    tgtposs.append(poss[i][j])

        # getdist
        dists = []
        for i in range(len(tgtposs)):
            dists.append(cobj.getdist_list(tgtposs[i], centerpos))
        print('dists', dists)


    f = open(ofile, 'w')
    for dist in dists:
        print('{:8.3f}'.format(dist), file=f)
    f.close()

    # check distribution
    print('--- histogram ---')
    ddatas = []
    for i in range(0, math.ceil(cell[0]/step)):
        dcount = 0
        for dist in dists:
            if (dist > step * i - (step/2.0)) and (dist < (step * i + (step/2.0))):
                dcount += 1
        print(step*i, dcount)
        ddatas.append([step*i, dcount])

    f = open(o2file, 'w')
    for ddata in ddatas:
        print('{0[0]:8.3f} {0[1]:>10}'.format(ddata), file=f)
    f.close()

