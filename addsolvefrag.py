import numpy as np
import sys
import os
import math
import subprocess
import re
import time
import copy
import ampt.cutset_fmo as cutf
import collections
# Matrix operation


if __name__ == "__main__":
    ## -- user setting --
    # read info
    mode = 'resnum' #rfile, resnum
    assignmolname = False
    refreshatmtype = False
    ajfname = '6lu7orig_hip_nowat_hinagata.ajf'
    solvname = ['WAT', 'NA']

    ## -- setting end --

    obj = cutf.cutset_fmo()
    obj.getmode = mode
    obj.assignmolname = assignmolname
    obj.refreshatmtype = refreshatmtype

    # main
    pdbnames = []
    for argv in sys.argv:
        if os.path.splitext(argv)[-1] == '.pdb':
            pdbnames.append(argv)

    for pdbname in pdbnames:

        print('infile:', pdbname)
        print('inajf', ajfname)

        # get pdbinfo
        totalMol, atomnameMol, molnameMol, anummols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfowrap(pdbname)

        # get tgt solvate mol info
        molname = solvname
        atomnumsets = []
        tgtmolsets = []
        for i in range(len(molname)):
            for j in range(totalMol):
            # for j in range(1):
                # print (molname[i], molnamelist_orig[j])
                if molname[i] == molnameMol[j]:
                    atomnumsets.append(len(posMol[j]))
                    tgtmolsets.append(molnameMol[j])
                    break
        print ('atomnumsets', atomnumsets)
        print('tgtmolsets', tgtmolsets)


        # get ajfinfo
        # fatomnumsets, fchgs, fbaas, fatminfos, connects = obj.getajfinfo(ajfname)
        # get fraginfo
        obj.getfragdict([ajfname], 'segment_data.dat')

        nameidMol = []
        for i in range(len(molnameMol)):
            for j in range(len(tgtmolsets)):
                if molnameMol[i] == tgtmolsets[j]:
                    nameidMol.append(j)
        # print(nameidMol)

        # add solute info
        tgtmolsets.append(os.path.splitext(ajfname)[0])
        nameidMol.insert(0, len(tgtmolsets)-1)
        atomnumsets.append(0)

        fatomnums, fchgs, fbaas, fatminfos, connects = obj.getfragtable(tgtmolsets, atomnumsets, nameidMol)
        # print (frag_atoms, frag_charges)

        # gen ajf file
        obj.ajf_method = "MP2"
        obj.ajf_basis_set = "6-31G*"
        obj.abinit_ver = 'rev15'
        obj.abinitmp_path = 'abinitmp'
        obj.pbcnv = 1.0
        obj.piedaflag = False
        obj.cpfflag = False
        obj.cmmflag = True
        obj.npro = 1
        obj.para_job = 1


        oname = os.path.splitext(pdbname)[0].split('/')[-1] + '_forabmp'
        readgeom = oname + '.pdb'
        ajf_body = obj.gen_ajf_bodywrap(fatomnums, fchgs, fbaas, connects, fatminfos, readgeom, ajfname)

        # export
        # ajf
        opath = 'for_abmp'
        ajf_oname = opath + '/' + os.path.splitext(pdbname)[0].split('/')[-1] + '-forabmp.ajf'

        print('ajf_oname:', ajf_oname)
        ajf_file = open(ajf_oname, 'w')

        print(ajf_body, file=ajf_file)

        # pdb
        if os.path.exists(opath) is False:
            print(opath)
            subprocess.call(["mkdir", opath])

        index = [i for i in range(len(posMol))]
        obj.exportardpdbfull(opath + '/' + oname, index, posMol, atomnameMol, molnameMol, heads, labs, chains, resnums, codes, occs, temps, amarks, charges)


    # print('ajf_fragment', ajf_fragment)
    # print('num_fragment', num_fragment)
    # print('ajf_charge', ajf_charge)

    # obj.bindfraginfo
    # obj.make_abinput_rmap(tgtmolnames, molnames, oname, opath, atomnums)

    # get ajfinfo
    # fatomnums, fchgs, fbaas, fatminfos, connects = obj.getajfinfo(ajfname)
    # print(fatomnums, fchgs, fbaas, fatminfos, connects)

#     param = [fatomnums], [fchgs], [fbaas], [connects], [fatminfos]
#     ajf_fragment = obj.get_fragsection(param)[:-1]
#     # print(ajf_fragment)
#
#     ajf_charge = sum(fchgs)
#     ajf_parameter = [ajf_charge, "", "", ajf_fragment, len(fatomnums), 0]
#
#     # gen ajf file
#     ajf_oname = ajfname + 'fragmod.ajf'
#     print('oname:', ajf_oname)
#     ajf_file = open(ajf_oname, 'w')
#
#     obj.ajf_method = "MP2"
#     obj.ajf_basis_set = "6-31G*"
#     obj.abinit_ver = 'rev10'
#     obj.abinitmp_path = 'abinitmp'
#     obj.pbmolrad = 'vdw'
#     obj.pbcnv = 1.0
#     obj.piedaflag = False
#
#     basis_str = obj.ajf_basis_set
#     print(basis_str)
#     if basis_str == '6-31G*':
#         basis_str = '6-31Gd'
#     ajf_parameter[1] = "'" +  pdbname + "'"
#     ajf_parameter[2] = "'" +  ajfname + '-' + obj.ajf_method + '-' + basis_str + ".cpf'"
#     ajf_body = obj.gen_ajf_body(ajf_parameter)
#     print(ajf_body, file=ajf_file)


