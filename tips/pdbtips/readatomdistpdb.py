import sys
import os
import numpy as np
import math
import statistics
import copy
import argparse
import glob
import shutil
import abmptools as ampt


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                prog='readcif.py', # program name
                usage='python readcif.py -c xxx.cif -odir dir', # program usage
                description='readcif script',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-c', '--coord',
                        help='coordinate file (pdb)',
                        nargs='*',
                        action='append',
                        required=True)

    parser.add_argument('-calcd', '--calcdist',
                        help='calcflag',
                        action='store_true',
                        )

    parser.add_argument('-intra', '--intra',
                        help='nointra',
                        action='store_true',
                        )

    parser.add_argument('-d', '--dist',
                        help='dist',
                        type=float,
                        default=2.0,
                        )

    parser.add_argument('-center', '--center',
                        help='center',
                        type=int,
                        nargs='*',
                        default=[],
                        )

    parser.add_argument('-a', '--tgtatom',
                        help='tgtatm',
                        nargs=2,
                        action='append',
                        default=[])

    parser.add_argument('-od', '--odir',
                        help='outdir',
                        default='cifout'
                        )

    parser.add_argument('-an', '--atomnum',
                        help='atom num',
                        nargs='*',
                        type=int,
                        action='append',
                        required=True)

    parser.add_argument('-l', '--layer',
                        help='layer',
                        type=int,
                        default=1
                        )

    parser.add_argument('-nopdb', '--nopdb',
                        help='flag nopdb',
                        action='store_false')

    parser.add_argument('-noout', '--noout',
                        help='flag nopdb',
                        action='store_true')

    parser.add_argument('-min', '--min',
                        help='min',
                        action='store_true')

    # get args
    args = parser.parse_args()

    print('coord(cif) =', args.coord)
    print('odir = ', args.odir)
    print('atomnum = ', args.atomnum)
    print('layer =', args.layer)
    print('pdbflag =', args.nopdb)
    print('calcdist', args.calcdist)
    print('intra', args.intra)
    print('dist', args.dist)
    print('tgtatom', args.tgtatom)
    print('out', args.noout)
    print('min', args.min)
    print('center', args.center)

    ## -- user setiing
    calc_dist = args.calcdist
    tgtdist = args.dist

    if args.tgtatom:
        tgtatoms = args.tgtatom[0]

    if args.intra:
        nointra = False
    else:
        nointra = True
    maxnum = 10

    ## -- user setting end
    # argvs = sys.argv
    # print(argvs)

    # infile = argvs[1]
    infiles = args.coord
    # odir = argvs[2]
    odir = args.odir

    anum_inmol = args.atomnum[0]
    image = args.layer
    pdbflag = args.nopdb

    print('atom num mol:', anum_inmol)
    print('########## Read Start #########')

    pdb = ampt.setfmo()
    pdb.solutes = args.center

    minvals = []
    mindatas = []
    for infile in infiles[0]:

        out, ext = os.path.splitext(infile)
        out = out.split('/')[-1]

        print('\n##  Start Read', infile)

        pdb.readpdb(infile)
        # print ('position\n', pdb.posRes)

        # 4. -- screening dist from solute --
        # i: residue(mol) index

        distdatas = []
        distvals = []

        for i in range(len(pdb.posRes)):
            # check ion
            if i % 100 == 0:
                print('check mol', i)
            # solute: skip
            if i+1 in pdb.solutes:
               continue
            # no solute: check dist between solute and no solute
            # j: atom loop in mol i
            for j in range(len(pdb.posRes[i])):
                # k: solute mol id
                # print(pdb.atmtypeRes[i][j])
                if pdb.atmtypeRes[i][j].strip() != tgtatoms[1]:
                    continue

                for solid in pdb.solutes:
                    k = solid -1
                    # l: solute atom id in mol k
                    for l in range(len(pdb.posRes[k])):
                        if pdb.atmtypeRes[k][l].strip() != tgtatoms[0]:
                            continue
                        dist = pdb.getdist(np.array(pdb.posRes[k][l]), np.array(pdb.posRes[i][j]))
                        if dist < tgtdist:
                            # print(pdb.atmtypeRes[i][j], pdb.atmtypeRes[k][l], dist)
                            print('mol', i+1, 'atom', j+1, pdb.atmtypeRes[i][j], '- mol', k+1, 'atom', l+1, pdb.atmtypeRes[k][l], "{:6.3f}".format(dist))
                            distdatas.append([i+1, j+1, pdb.atmtypeRes[i][j], k+1, l+1, pdb.atmtypeRes[k][l]])
                            distvals.append(dist)


        if args.min:
            minidx = distvals.index(min(distvals))
            minval = min(distvals)
            mindata = distdatas[minidx]
            print('minval', minval)
            print('distdata', mindata)
            minvals.append(minval)
            mindatas.append(['file', infile, mindata])
        # print('as_list', as_list)
        # print('len as_list', len(as_list))

        ## write for file

    if args.min:
        ofile = 'mindist-' + tgtatoms[0] + '-' + tgtatoms[1] + '.log'
        fomin = open(ofile, 'w')
        for i in range(len(minvals)):
            print(mindatas[i], minvals[i], file=fomin)
        print('generated' + ofile)
# 

