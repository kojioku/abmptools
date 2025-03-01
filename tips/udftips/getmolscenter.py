from UDFManager import *
import sys
import math
# import fcewsmb.udf_io as uio
import abmptools as ampt
import numpy as np

if __name__ == "__main__":
    fname = sys.argv[1]
    mode = 'particle' # mol, particle
    tgtnames = ["Cholesterol", "DOPC"]

    record = 101

    step = 1.0
    # mol = 0
    # print (fname, record, mol)
    print('target:', tgtnames)
    uobj = UDFManager(fname)
    cobj = ampt.udfrm_io()


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
        if molnames[i] in tgtnames:
            targets.append(i)
    print('target id', targets)

    poss = []
    for tgt in targets:
        poss.append(cobj.getposmol(uobj, tgt))
    # print('poss0', poss[0])
    print('num of tgt', len(poss))

    cocs = []
    for i in range(len(poss)):
        cocs.append(cobj.getCenter(poss[i]).tolist())
    # print('cocs', cocs)
#
    print('tgtcenter:', cobj.getCenter(cocs))
