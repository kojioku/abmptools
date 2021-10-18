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
import argparse
# Matrix operation


if __name__ == "__main__":
    ## -- user setting --
    parser = argparse.ArgumentParser(
                prog='addsolvfrag.py', # program name
                usage='''e.g.)
                python addsolvfrag.py -temp 6lu7-covneu-nowat-hinagata0514.ajf -solv HOH WAT NA -i *.pdb
                ''', # program usage
                description='description',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--incoord',
                        help='input frag info',
                        nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-temp', '--template',
                        help='template file',
                        required=True)

    parser.add_argument('-solv', '--solvmol',
                        help='output config file',
                        nargs='*',
                        default=['HOH', 'WAT', 'NA'])

    ###################################

    parser.add_argument('-pb', '--solvation',
                        help='parameter file',
                        action='store_true')

    parser.add_argument('-th', '--thicnv',
                        help='threshold',
                        type=float,
                        default=1.0E-5
                        )

    parser.add_argument('-arad', '--atmrad',
                        help='atmrad',
                        default='delphi'
                        )

    parser.add_argument('-ajfv', '--ajfversion',
                        help='ajf version',
                        default='rev23'
                        )

    parser.add_argument('-np', '--npro',
                        help='ajf version',
                        type=int,
                        default=4
                        )

    parser.add_argument('-nopieda', '--nopieda',
                        help='pieda',
                        action='store_false'
                        )

    parser.add_argument('-cmm', '--cmm',
                        help='cmm',
                        action='store_true'
                        )

    parser.add_argument('-nocpf', '--nocpf',
                        help='cpfflag',
                        action='store_false'
                        )

    parser.add_argument('-cpfv', '--cpfver',
                        help='cpf version',
                        default='10'
                        )

    parser.add_argument('-basis', '--basisset',
                        help='basis',
                        default='6-31G*',
                        )

    parser.add_argument('-m', '--method',
                        help='method',
                        default='MP2',
                        )

    parser.add_argument('-ml', '--mldat',
                        help='mldat flag',
                        action='store_true'
                        )

    parser.add_argument('-mll', '--mllimit',
                        help='mldat fraglimit',
                        type=int,
                        default=0
                        )

    parser.add_argument('-disp', '--disp',
                        help='flag disp',
                        action='store_true')

    # WriteMLdata='wstr-1E08_HIS_ES.new2.cmm5.mldat'
    # MLfraglimit=921

    parser.add_argument('-dg', '--dgemm',
                        help='dgemm',
                        action='store_true',
                        )

    parser.add_argument('-rp', '--resp',
                        help='resp',
                        action='store_true',
                        )

    parser.add_argument('-nonbo', '--nonbo',
                        help='nonbo',
                        action='store_false',
                        )

    parser.add_argument('-mem', '--memory',
                        help='memory',
                        default='3000',
                        )

    parser.add_argument('-lc', '--ligandcharge',
                        help='ligand charge',
                        nargs=2,
                        action='append',
                        # default='',
                        )

    parser.add_argument('-rs', '--rsolv',
                        help='rsolv',
                        nargs=2,
                        action='append',
                        # default='',
                        )

    parser.add_argument('-ma', '--manual',
                        help='manual table',
                        action='store_true',
                        )

    parser.add_argument('-bsse', '--bsse',
                        help='bsse',
                        action='store_true',
                        )

    # get args
    args = parser.parse_args()

    print('coord(pdb) =', args.incoord)
    print('solv = ', args.solvation)
    print('pbcnv = ', args.thicnv)
    print('atmrad =', args.atmrad)
    print('ajfversion =', args.ajfversion)
    print('np =', args.npro)
    print('pieda =', args.nopieda)
    print('cmm =', args.cmm)
    print('cpf =', args.nocpf)
    print('cpfver =', args.cpfver)
    print('basis =', args.basisset)
    print('method =', args.method)
    print('dgemm', args.dgemm)
    print('resp', args.resp)
    print('nbo', args.nonbo)
    print('memory', args.memory)
    print('ligand charge', args.ligandcharge)
    print('rsolv', args.rsolv)
    print('manual', args.manual)
    print('bsse', args.bsse)
    print('mldat', args.mldat)
    print('disp', args.disp)


    ####################################

    args = parser.parse_args()
    # print('output =', args.output)

    ajfname = args.template
    solvname = args.solvmol

    ## -- setting end --
    # read info
    mode = 'resnum' #rfile, resnum
    assignmolname = False
    refreshatmtype = False
    refreshresid = False

    aobj = ampt.setfmo()
    aobj.getmode = mode
    aobj.assignmolname = assignmolname
    aobj.refreshatmtype = refreshatmtype
    aobj.refreshresid = refreshresid

    # gen ajf file
    aobj.ajf_method = args.method
    aobj.ajf_basis_set = args.basisset
    aobj.abinit_ver = args.ajfversion
    aobj.autofrag = True
    aobj.piedaflag = args.nopieda
    aobj.cpfflag = args.nocpf
    aobj.cpfver = args.cpfver
    aobj.cmmflag = args.cmm
    aobj.npro = args.npro
    aobj.pbmolrad = args.atmrad
    # aobj.readgeom = args.incoord
    aobj.solv_flag = args.solvation
    aobj.pbcnv = args.thicnv
    aobj.dgemm = args.dgemm
    aobj.resp = args.resp
    aobj.nbo = args.nonbo
    aobj.memory = args.memory
    aobj.ligchg = args.ligandcharge
    aobj.rsolv = args.rsolv
    aobj.bsseflag = args.bsse
    aobj.disp = args.disp

    if args.mldat:
        aobj.mldatfrag = args.mldat
        aobj.mllimit = args.mllimit


    if args.manual:
        aobj.autofrag = False

    # main
    pdbnames = []
    for argv in sys.argv:
        if os.path.splitext(argv)[-1] == '.pdb':
            pdbnames.append(argv)

    for pdbname in pdbnames:

        print('infile:', pdbname)
        print('inajf', ajfname)

        # get pdbinfo
#         totalMol, atomnameMol, self.resnames, anummols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = aobj.getpdbinfowrap(pdbname)
        aobj = aobj.getpdbinfowrap(pdbname)

        # get tgt solvate mol info
        molname = solvname
        atomnumsets = []
        tgtmolsets = []
        for i in range(len(molname)):
            for j in range(aobj.totalRes):
            # for j in range(1):
                # print (molname[i], molnamelist_orig[j])
                if molname[i] == aobj.resnames[j]:
                    atomnumsets.append(len(aobj.posRes[j]))
                    tgtmolsets.append(aobj.resnames[j])
                    break
        print ('atomnumsets', atomnumsets)
        print('tgtmolsets', tgtmolsets)

        # get ajfinfo
        # fatomnumsets, fchgs, fbaas, fatminfos, connects = aobj.getajfinfo(ajfname)
        # get fraginfo
        aobj.getfragdict([ajfname], 'segment_data.dat')

        nameidMol = []
        for i in range(len(aobj.resnames)):
            for j in range(len(tgtmolsets)):
                if aobj.resnames[i] == tgtmolsets[j]:
                    nameidMol.append(j)
        # print(nameidMol)

        # add solute info
        tgtmolsets.append(os.path.splitext(ajfname)[0].split('/')[-1])
        nameidMol.insert(0, len(tgtmolsets)-1)
        atomnumsets.append(0)

        # fatomnums, fchgs, fbaas, fatminfos, connects = aobj.getfragtable(tgtmolsets, atomnumsets, nameidMol)
        aobj = aobj.getfragtable(tgtmolsets, atomnumsets, nameidMol)
        # print (frag_atoms, frag_charges)

        # gen ajf file
        # aobj.ajf_method = "MP2"
        # aobj.ajf_basis_set = "6-31G*"
        # aobj.abinit_ver = 'rev23'
        # aobj.piedaflag = True
        # aobj.npro = 1
        # aobj.para_job = 1
        fname = os.path.splitext(pdbname)[0].split('/')[-1]
        ohead = fname + '_forabmp'
        aobj.readgeom = ohead + '.pdb'
        aobj.writegeom = fname + '-' + aobj.ajf_method + '-' + aobj.ajf_basis_set.replace('*', 'd') + ".cpf'"
        if aobj.mllimit == 0:
            aobj.mldatname = "'" + fname + ".mldat'"
        else:
            aobj.mldatname = "'" + fname + "limit" + str(aobj.mllimit) + ".mldat'"

        opath = 'for_abmp'
        if os.path.exists(opath) is False:
            print(opath)
            subprocess.call(["mkdir", opath])

        ajf_oname = opath + '/' + ohead + '.ajf'
        aobj.saveajf(ajf_oname)

        # exportpdb
        index = [i for i in range(len(aobj.posRes))]
        aobj.exportardpdbfull(opath + '/' + ohead, index)


