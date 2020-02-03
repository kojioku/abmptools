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
    # read info
    mode = 'rfile' #rfile, resnum
    calcmode = 'molpair' #molpair, id
    step = 1.0
    assignmolname = False
    refreshatmtype = False

    # tgtmol = 2
    tgt1 = '000'
    tgt2 = '000'
    ## -- setting end --

    centermolid = 20
    tgtmol = '000'
    # main
    argvs = sys.argv
    # fname = str(argvs[1])

    if calcmode == 'molpair':
        ofile= 'dist-' + tgt1 + '-' + tgt2 + '.txt'
        o2file= 'hist-' + tgt1 + '-' + tgt2 + '.txt'
        o3file= 'rdf-' + tgt1 + '-' + tgt2 + '.txt'


    if calcmode == 'id':
        ofile= 'dist-mol' + str(centermolid) + '-' + tgtmol + '.txt'
        o2file= 'hist-mol' + str(centermolid) + '-' + tgtmol + '.txt'
        o3file= 'rdf-mol' + str(centermolid) + '-' + tgtmol + '.txt'



#     if mode == 'particle':
#         ofile= 'dist-' + tgtname + '-' + tgtpart + 'from' + str(centerpos[0]) + '-' + str(centerpos[1]) + '-' + str(centerpos[2]) + 'rec' + str(record) + '.txt'
#         o2file= 'hist-' + tgtname + '-' + tgtpart + 'from' + str(centerpos[0]) + '-' + str(centerpos[1]) + '-' + str(centerpos[2]) + 'rec' + str(record) + '.txt'

    for arg in argvs:
        if arg == '--move':
            moveflag == True
        if arg == '--nomove':
            moveflag == False

    for ii in range(len(argvs)):
        if ii == 0:
            continue
        if argvs[ii][0:2] == '--':
            continue
        fname = argvs[ii]
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
        # print('centered-molid:', tgtmol - 1)

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

        if calcmode == 'molpair':
            targets = []
            for i in range(totalMol):
                # print(molnames[i])
                if molnames[i] == tgt1:
                    targets.append(i)
            print('target1 id', targets)

            targets2 = []
            for i in range(totalMol):
                # print(molnames[i])
                if molnames[i] == tgt2:
                    targets2.append(i)
            print('target2 id', targets2)

        # centerOfMol: com of each molecules
            cocs = []
            for i in targets:
                cocs.append(obj.getCenter(posMol[i]).tolist())
            # print('cocs', cocs)

            cocs2 = []
            for i in targets2:
                cocs2.append(obj.getCenter(posMol[i]).tolist())
            # print('cocs', cocs)

            # getdist
            dists = []
            for i in range(len(targets)):
                for j in range(i, len(targets2)):
                    if targets[i] == targets[j]:
                        continue
                    dists.append(obj.getdist_list(cocs[i], cocs2[j]))
            print('dists', dists)

            print(len(cocs), len(cocs2), len(dists))

        if calcmode == 'id':
            targets = []
            for i in range(totalMol):
                # print(molnames[i])
                if molnames[i] == tgtmol:
                    targets.append(i)
            print('tgtmol id', targets)
            print('coc', posMol[centermolid])


        # centerOfMol: com of each molecules
            cocs = []
            for i in targets:
                cocs.append(obj.getCenter(posMol[i]).tolist())
            # print('cocs', cocs)


            # getdist
            dists = []
            for i in range(len(targets)):
                dists.append(obj.getdist_list(cocs[i], cocs[centermolid]))
            print('dists', dists)

            print(len(dists))

#
#
# #         if calcmode == 'particle':
# #             atoms = []
# #             for tgt in targets:
# #                 atoms.append(obj.getnameAtom(uobj, tgt))
# #             # print(atoms)
# #
# #             tgtposMol = []
# #             for i in range(len(posMol)):
# #                 for j in range(len(atoms[i])):
# #                     if atoms[i][j] == tgtpart:
# #                         tgtposMol.append(posMol[i][j])
# #
# #             # getdist
# #             dists = []
# #             for i in range(len(tgtposMol)):
# #                 dists.append(obj.getdist_list(tgtposMol[i], centerpos))
# #             print('dists', dists)
#
#
        f = open(ofile, 'w')
        for dist in dists:
            print('{:8.3f}'.format(dist), file=f)
        f.close()

        # check distribution
        print('--- histogram ---')
        ddatas = []
        for i in range(0, math.ceil(100/step)):
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

        # get areavol
        areas = []
        rdfs = []
        astart = 0
        for i in range(0, math.ceil(100/step)):
            rend = step * i + (step/2.0)
            aend = 4.0/3.0 * math.pi * rend ** 3
            areas.append(aend - astart)
            rdfs.append([ddatas[i][0], ddatas[i][1]/(aend-astart)])
            astart = copy.deepcopy(aend)

        print(areas)

        print('--- rdf ---')
        print(rdfs)
        f = open(o3file, 'w')
        for rdf in rdfs:
            print('{0[0]:8.3f} {0[1]:>10}'.format(rdf), file=f)
        f.close()

        print('out:', ofile, o2file, o3file)

