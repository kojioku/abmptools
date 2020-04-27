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
    ## -- setting end --

    obj = cutf.cutset_fmo()
    obj.getmode = mode
    obj.assignmolname = assignmolname
    obj.refreshatmtype = refreshatmtype

    # main
    argvs = sys.argv
    # fname = str(argvs[1])

    pdbname = sys.argv[1]
    ajfname = sys.argv[2]

    print('infile:', pdbname)

    # get pdbinfo
    totalMol, atomnameMol, molnames, anummols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges = obj.getpdbinfowrap(pdbname)
    # get ajfinfo
    fatomnums, fchgs, fbaas, fatminfos, connects = obj.getajfinfo(obj, ajfname)
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

    fatomnums, fchgs, fbaas, connects, fatminfos = obj.modifyfragparam(totalMol, atomnameMol, molnames, anummols, posMol, heads, labs, chains
                                                                   ,resnums ,codes ,occs ,temps ,amarks ,charges, fatomnums, fchgs, fbaas, fatminfos, connects, bridgeds, doubles, nagatoms)

#     print('new_totalres', len(fatomnums))
#     print(len(fchgs))
#     print(len(fbaas))
#     print(len(fatminfos))
#     print(sum(fbaas))
#     print(len(connects))

    param = [fatomnums], [fchgs], [fbaas], [connects], [fatminfos]
    ajf_fragment = obj.get_fragsection(param)[:-1]
    # print(ajf_fragment)

    ajf_charge = sum(fchgs)
    ajf_parameter = [ajf_charge, "", "", ajf_fragment, len(fatomnums), 0]

    # gen ajf file
    ajf_oname = os.path.splitext(ajfname)[0] + 'fragmod.ajf'
    print('oname:', ajf_oname)
    ajf_file = open(ajf_oname, 'w')

    obj.ajf_method = "MP2"
    obj.ajf_basis_set = "6-31G*"
    obj.abinit_ver = 'rev10'
    obj.abinitmp_path = 'abinitmp'
    obj.pbmolrad = 'vdw'
    obj.pbcnv = 1.0
    obj.piedaflag = False

    basis_str = obj.ajf_basis_set
    print(basis_str)
    if basis_str == '6-31G*':
        basis_str = '6-31Gd'
    ajf_parameter[1] = "'" +  pdbname + "'"
    ajf_parameter[2] = "'" +  os.path.splitext(ajfname)[0] + '-' + obj.ajf_method + '-' + basis_str + ".cpf'"
    ajf_body = obj.gen_ajf_body(ajf_parameter)
    print(ajf_body, file=ajf_file)


