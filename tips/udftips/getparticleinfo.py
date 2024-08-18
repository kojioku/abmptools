from UDFManager import *
import sys
import copy
import abmptools
import numpy as np

if __name__ == "__main__":
    fname = sys.argv[1]
    tgtname = "Cholesterol"

    # particle_id = [3, 4]
    particle_id = [3, 4, 5]
    record = 101
    tgtpos = [5.0, 5.0, 5.0]
    area = 5.0

    print (fname, record)
    print('target:', tgtname)
    uobj = UDFManager(fname)
    cobj = abmptools.udfrm_io()
    if len(particle_id) == 2:
        ofile= 'dist-' + tgtname + 'part' + str(particle_id[0]) + '-' + str(particle_id[1]) + 'rec' + str(record) + '.txt'
    elif len(particle_id) == 3:
        ofile= 'angle-' + tgtname + 'part' + str(particle_id[0]) + '-' + str(particle_id[1]) + '-' + str(particle_id[2]) + 'rec' + str(record) + '.txt'

    totalmol, totalrec, totalatm, cell = cobj.getudfinfowrap(uobj)

    uobj.jump(record)
    cobj.moveintocell_rec(uobj, record, totalmol)

    # get molname list
    molnames = []
    for i in range(totalmol):
        # if i % 10 == 0:
        #     print("read mol name", i)
        molnames.append(cobj.getmolname(i, uobj))
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
    print('poss0', poss[0])
    print('num of tgt', len(poss))

#     # centerOfMol: com of each molecules
    # limit area mode
    if len(tgtpos) != 0:
        if len(particle_id) == 2:
            ofile= 'dist-' + tgtname + 'part' + str(particle_id[0]) + '-' + str(particle_id[1]) + 'rec' + str(record) + 'pos' + str(tgtpos[0]) + str(tgtpos[1]) + str(tgtpos[1]) + '-area' + str(area) + '.txt'
        elif len(particle_id) == 3:
            ofile= 'angle-' + tgtname + 'part' + str(particle_id[0]) + '-' + str(particle_id[1]) + '-' + str(particle_id[2]) + 'rec' + str(record) + 'pos' + str(tgtpos[0]) + str(tgtpos[1]) + str(tgtpos[1]) + '-area' + str(area) + '.txt'

        print('get center of coordinate')
        print('tgtpos', tgtpos)
        print('area', area)
        cocs = []
        for i in range(len(poss)):
            cocs.append(cobj.getCenter(poss[i]).tolist())
        # print('cocs', cocs)

        tdists = []
        targets_new = []
        poss_new = []
        for i in range(len(targets)):
            tdist = cobj.getdist_list(cocs[i], tgtpos)
            tdists.append(tdist)
            if tdist <= area:
                # print (targets[i], 'is in the area')
                targets_new.append(targets[i])
                poss_new.append(poss[i])
        # print('tdists', tdists)

        targets = []
        poss = []
        targets = copy.deepcopy(targets_new)
        poss = copy.deepcopy(poss_new)
        print('area target =', targets)
    if len(particle_id) == 2:
        dists = []
        for i in range(len(targets)):
            p1 = particle_id[0]
            p2 = particle_id[1]
            dists.append(cobj.getdist_list(poss[i][p1], poss[i][p2]))
        print('dists', dists)

        f = open(ofile, 'w')
        for dist in dists:
            print(dist, file=f)


    elif len(particle_id) == 3:
        angles = []
        for i in range(len(targets)):
            p1 = particle_id[0]
            p2 = particle_id[1]
            p3 = particle_id[2]
            l1 = cobj.getdist_list(poss[i][p2], poss[i][p1])
            l2 = cobj.getdist_list(poss[i][p2], poss[i][p3])
            angles.append(cobj.getangle(poss[i][p1], poss[i][p2], poss[i][p3], [l1, l2]))
        print('angles', angles)

        f = open(ofile, 'w')
        for angle in angles:
            print(angle, file=f)


