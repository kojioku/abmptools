from UDFManager import *
import sys
import abmptools
import numpy as np

if __name__ == "__main__":
    fname = sys.argv[1]
    # tgtname = "Cholesterol"
    tgtname = "POPC"

    tgtmolid = 10
    contactdist = 5.0
    # tgtname = "DOPC"
    # record = 1
    # mol = 0
    # print (fname, record, mol)
    print('target:', tgtname, "mol", tgtmolid)
    uobj = UDFManager(fname)
    cobj = abmptools.udfrm_io()

    totalmol, totalrec, totalatm, cell = cobj.getudfinfowrap(uobj)

    # get molname list
    molnames = []
    for i in range(totalmol):
        # if i % 10 == 0:
        #     print("read mol name", i)
        molnames.append(cobj.getmolname(i, uobj))
    print(totalmol, totalrec)
    # print(molnames)

    # get target mol index
    targets = []
    for i in range(totalmol):
        if molnames[i] == tgtname:
            targets.append(i)
    print('target id', targets)

    # get target posmol for all rec
    poss = []
    mol = targets[tgtmolid]
    print("--- target mol id:", mol, "---")

    for record in range(totalrec):
        poss.append(cobj.getposmolrec(uobj, mol, record))
    # print('poss0', poss[0])

    # centerOfMol: com of each molecules
    cocs = []
    for i in range(len(poss)):
        cocs.append(cobj.getCenter(poss[i]).tolist())
    # print('cocs', cocs)

    # getdist
    dists = []
    for i in range(totalrec-1):
        dists.append(cobj.getdist_list(cocs[i], cocs[i+1]))
    print('dists', dists)

    # get cholesterol mol pos at one record

    countss = []
    for record in range(totalrec):
        if record % 10 == 0:
            print('rec', record)
        tgtposs = []
        for mol in targets:
            tgtposs.append(cobj.getposmolrec(uobj, mol, record))
        # print('tgtposs0', tgtposs[0])

        # centerOfMol: com of each molecules
        tgtcocs = []
        for i in range(len(tgtposs)):
            tgtcocs.append(cobj.getCenter(tgtposs[i]).tolist())
        # print('tgtcocs', tgtcocs)

        for k in range(3):
            # print('take', k)
            for i in range(len(tgtcocs)):
                for j in range(3):
                    if tgtcocs[i][j] > cell[0]:
                        # print("big", tgtcocs[i][j])
                        tgtcocs[i][j] -= cell[0]
                    if tgtcocs[i][j] < 0:
                        # print("small", tgtcocs[i][j])
                        tgtcocs[i][j] += cell[0]

        # # check contact
        counts = []
        for i in range(len(tgtposs)):
            count = 0
            for j in range(len(tgtposs)):
                if i == j:
                    continue
                dist = cobj.getdist_list(tgtcocs[i], tgtcocs[j])
                # print(i, j, dist)
                if dist <= contactdist:
                    # print(i, j, "contact")
                    count += 1
            counts.append(count)
        countss.append(counts)

    onecount = []
    for i in range(len(countss)):
        onecount.append(countss[i][tgtmolid])
    print("number_mol_around", onecount)
