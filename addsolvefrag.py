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

    aobj = ampt.setfmo()
    aobj.getmode = mode
    aobj.assignmolname = assignmolname
    aobj.refreshatmtype = refreshatmtype

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
            for j in range(self.totalRes):
            # for j in range(1):
                # print (molname[i], molnamelist_orig[j])
                if molname[i] == self.resnames[j]:
                    atomnumsets.append(len(self.posRes[j]))
                    tgtmolsets.append(self.resnames[j])
                    break
        print ('atomnumsets', atomnumsets)
        print('tgtmolsets', tgtmolsets)


        # get ajfinfo
        # fatomnumsets, fchgs, fbaas, fatminfos, connects = aobj.getajfinfo(ajfname)
        # get fraginfo
        aobj = aobj.getfragdict([ajfname], 'segment_data.dat')

        nameidMol = []
        for i in range(len(self.resnames)):
            for j in range(len(tgtmolsets)):
                if self.resnames[i] == tgtmolsets[j]:
                    nameidMol.append(j)
        # print(nameidMol)

        # add solute info
        tgtmolsets.append(os.path.splitext(ajfname)[0])
        nameidMol.insert(0, len(tgtmolsets)-1)
        atomnumsets.append(0)

        # fatomnums, fchgs, fbaas, fatminfos, connects = aobj.getfragtable(tgtmolsets, atomnumsets, nameidMol)
        aobj = aobj.getfragtable(tgtmolsets, atomnumsets, nameidMol)
        # print (frag_atoms, frag_charges)

        # gen ajf file
        aobj.ajf_method = "MP2"
        aobj.ajf_basis_set = "6-31G*"
        aobj.abinit_ver = 'rev15'
        aobj.piedaflag = False
        aobj.cpfflag = False
        aobj.cmmflag = True
        aobj.npro = 1
        aobj.para_job = 1
        ohead = os.path.splitext(pdbname)[0].split('/')[-1] + '_forabmp'
        aobj.readgeom = ohead + '.pdb'

        basis_str = self.ajf_basis_set
        print(basis_str)
        if basis_str == '6-31G*':
            basis_str = '6-31Gd'

        aobj.writegeom = os.path.splitext(ajfname)[0] + '-' + self.ajf_method + '-' + basis_str + ".cpf'"

        ajf_body = aobj.gen_ajf_bodywrap(aobj)

        # export
        # ajf
        opath = 'for_abmp'
        ajf_oname = opath + '/' + ohead + '.ajf'
        if os.path.exists(opath) is False:
            print(opath)
            subprocess.call(["mkdir", opath])

        print('ajf_oname:', ajf_oname)
        ajf_file = open(ajf_oname, 'w')

        print(ajf_body, file=ajf_file)

        # pdb

        index = [i for i in range(len(self.posMol))]
        aobj.exportardpdbfull(opath + '/' + self.readgeom, index)


    # print('ajf_fragment', ajf_fragment)
    # print('num_fragment', num_fragment)
    # print('ajf_charge', ajf_charge)

    # aobj.bindfraginfo
    # aobj.make_abinput_rmap(tgtmolnames, molnames, oname, opath, atomnums)

    # get ajfinfo
    # fatomnums, fchgs, fbaas, fatminfos, connects = aobj.getajfinfo(ajfname)
    # print(fatomnums, fchgs, fbaas, fatminfos, connects)

#     param = [fatomnums], [fchgs], [fbaas], [connects], [fatminfos]
#     ajf_fragment = aobj.get_fragsection(param)[:-1]
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
#     aobj.ajf_method = "MP2"
#     aobj.ajf_basis_set = "6-31G*"
#     aobj.abinit_ver = 'rev10'
#     aobj.abinitmp_path = 'abinitmp'
#     aobj.pbmolrad = 'vdw'
#     aobj.pbcnv = 1.0
#     aobj.piedaflag = False
#
#     basis_str = aobj.ajf_basis_set
#     print(basis_str)
#     if basis_str == '6-31G*':
#         basis_str = '6-31Gd'
#     ajf_parameter[1] = "'" +  pdbname + "'"
#     ajf_parameter[2] = "'" +  ajfname + '-' + aobj.ajf_method + '-' + basis_str + ".cpf'"
#     ajf_body = aobj.gen_ajf_body(ajf_parameter)
#     print(ajf_body, file=ajf_file)


