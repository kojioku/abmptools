import sys
import os
scrdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scrdir)

from multiprocessing import Pool
import copy
import random
import numpy as np
import math
import re
import subprocess
from ctypes import *
# import setparam as sp
import mol_io as mi
try:
    import collections
except:
    pass
try:
    import pandas as pd
except:
    pass
try:
    import itertools
except:
    pass


class abinit_io(mi.mol_io):

    def __init__(self):
        # print('## load abinit io')
        super().__init__()
        # gen_rand: gen_coord

        # submit param
        self.npro = 4
        self.memory = "1800"
        self.queue = None
        self.solv_flag = False
        self.abinitmp_path = 'abinitmp'
        self.mpipath = 'mpirun'

        # abinitmp param
        self.ajf_method = "MP2"
        self.ajf_basis_set = "6-31G*"
        self.abinit_ver = 'rev22' #rev10, 11, 15, 17, mizuho
        self.pbmolrad = 'vdw'
        self.pbcnv = 1.0
        self.piedaflag = True
        self.cpfflag = False
        self.cpfver='23'
        self.cmmflag = False
        self.nofzc_flag = False
        self.bsseflag = False
        self.readgeom = ""
        self.writegeom = ""
        self.submit_system = None
        self.autofrag = False
        self.dgemm = False
        self.nbo = True
        self.resp = True
        self.ligchg = None
        # distlitname

        # frag table
        self.fatomnums = []
        self.fchgs = []
        self.fbaas = []
        self.fatminfos = []
        self.connects = []

        return

    def chkdepth(self, k):
        if not k:
            return 0
        else:
            if isinstance(k, list):
                return 1 + max(self.chkdepth(i) for i in k)
            else:
                return 0

    def get_fragsection(self):

        if self.chkdepth(self.fatomnums) == 1:
            frag_atom = [self.fatomnums]
            frag_charge = [self.fchgs]
            frag_connect_num = [self.fbaas]
            frag_connect = [self.connects]
            seg_info = [self.fatminfos]
        else:
            frag_atom = self.fatomnums
            frag_charge = self.fchgs
            frag_connect_num = self.fbaas
            frag_connect = self.connects
            seg_info = self.fatminfos

        ajf_fragment = ''
        # fragment atom
        for i in range(len(frag_atom)):
            icount = 0
            for j in range(len(frag_atom[i])):
                icount += 1
                ajf_fragment += '%8d' % (frag_atom[i][j])
                if icount % 10 == 0:
                    icount = 0
                    ajf_fragment += '\n'
        if icount % 10 != 0:
            ajf_fragment += '\n'

        # fragment charge
        for i in range(len(frag_charge)):
            icount = 0
            for j in range(len(frag_charge[i])):
                icount += 1
                ajf_fragment += '%8d' % (frag_charge[i][j])
                if icount % 10 == 0:
                    icount = 0
                    ajf_fragment += '\n'
        if icount % 10 != 0:
            ajf_fragment += '\n'

        # fragment connect num
        for i in range(len(frag_connect_num)):
            icount = 0
            for j in range(len(frag_connect_num[i])):
                icount += 1
                ajf_fragment += '%8d' % (frag_connect_num[i][j])
                if icount % 10 == 0:
                    icount = 0
                    ajf_fragment += '\n'
        if icount % 10 != 0:
            ajf_fragment += '\n'

        # fragment body
        atom_count = 0
        # i: mol j: frag, k:frag atom
        for i in range(len(seg_info)):
            for j in range(len(seg_info[i])):
                icount = 0
                for k in range(len(seg_info[i][j])):
                    icount += 1
                    ajf_fragment += '%8d' % (seg_info[i][j][k] + atom_count)
                    # print(seg_info[i][j][k] + atom_count, 'icount', icount)
                    if icount % 10 == 0:
                        icount = 0
                        ajf_fragment += '\n'

                if icount % 10 != 0:
                    ajf_fragment += '\n'
            atom_count += sum(frag_atom[0])

        atom_count = 0
        # connect info
        for i in range(len(frag_connect)):
            for j in range(len(frag_connect[i])):
                ajf_fragment += '%8d' % (frag_connect[i][j][0] + atom_count)
                ajf_fragment += '%8d' % (frag_connect[i][j][1] + atom_count)
                ajf_fragment += '\n'
            atom_count += sum(frag_atom[0])
        return ajf_fragment


    def gen_ajf_bodywrap(self, inname):
        # i, j ,k (mol, frag, component_perfrag) dimension
        if self.autofrag == True:
            ajf_fragment = ""
        else:
            ajf_fragment = self.get_fragsection()[:-1]

        ajf_charge = sum(self.fchgs)
        ajf_parameter = [ajf_charge, ajf_fragment, len(self.fatomnums), 0]

        if self.writegeom == "":
            self.writegeom = "'" + os.path.splitext(inname)[0] + '-' + self.ajf_method + '-' + self.ajf_basis_set.replace('*', 'd') + ".cpf'"

        ajf_body = self.gen_ajf_body(ajf_parameter)

        return ajf_body


    def gen_ajf_body(self, param, fzctemp="YES"):
        ajf_charge, ajf_fragment, num_fragment, MOME_flag = param

        if self.autofrag == True:
            auto = 'ON'
        else:
            auto = 'OFF'


        if MOME_flag is True:
            frag_section = """
NF= """ + str(num_fragment) + """
LMOTYP='MOME' """
        else:
            frag_section = """
NF= """ + str(num_fragment) + """
LMOTYP='ANO' """

        if num_fragment == 1:
            use_FMO = 'OFF'
        else:
            use_FMO = 'ON'

        if self.submit_system in ['K', 'OFP']:
            np = '1'
        else:
            np = str(int(self.npro/self.para_job))

        # CNTRL NAMELIST
        ajf_body = """&CNTRL
ElecState='S1'
Method='""" + str(self.ajf_method) + """'
Nprint=3
Memory=""" + str(self.memory) + """
Natom=0
Charge=""" + str(ajf_charge) + """
ReadGeom='""" + str(self.readgeom) + "'\n"

        if self.cpfflag == True:
            ajf_body += "WriteGeom=" + str(self.writegeom) + "\n"
        ajf_body += """Gradient='NO'
Vector='OFF'
CPFBIN='NO'
CPFVER=""" + str(self.cpfver) + """
THOVL=1.0E-12
E_THSWZ=1.0E-12
G_THSWZ=1.0E-12
/

&FMOCNTRL
FMO='""" + str(use_FMO) + """'
NBody=2
AutoFrag='"""+ auto + """' """ + str(frag_section) + """
Laoc=0.0
Lptc=2.0
esp_ptc_multipole='NO'
Ldimer=2.0
NP=""" + np + """
MaxSCCcyc=250
MaxSCCenergy=5.0E-7
"""
        if  self.ligchg != None:
            ligchgstr = "LigandCharge='"
            for i in range(len(self.ligchg)):
                ligchgstr += self.ligchg[i][0] + '=' + self.ligchg[i][1]
                if i == len(self.ligchg) - 1:
                    ligchgstr += "'"
                else:
                    ligchgstr += ','

            print(ligchgstr)
            ajf_body += ligchgstr + """
"""
        if self.cmmflag == True:
            ajf_body += """Dimer_es_multipole='YES'
Ldimer_cmm=5.0
/
"""
        else:
            ajf_body += """/
"""

        ajf_body += """
&SCF
MaxSCFenergy=1.0E-8
MaxSCFdensity=1.0E-6
MaxSCFcyc=150
DIISTYPE='C2_OLD'
THINTEG=1.0E-12
IFCD='NO'
/

&BASIS
BasisSet='""" + str(self.ajf_basis_set) + """'
DiffuseOn='NO'
/

&OPTCNTRL
OPT='OFF'
/

"""
        if self.abinit_ver in ['rev22', 'rev23']:
            ajf_body += """
&MFMO
/
"""

        ajf_body += """&SCZV
DimerResponseTerm='NO'
/
"""
        if self.abinit_ver in ['rev22', 'rev23']:
            ajf_body += """
&XUFF
/
"""

        ajf_body += """
&MP2"""

        if self.dgemm is True:
            ajf_body += """
MOD1ST='GEMM'
MOD2ND='GEMM'
MOD3RD='GEMM'
MOD4TH='GEMM'"""

        ajf_body += """
NP_MP2_IJ=1
NP_MP2_S=0
MemoryMP2=0
IFSCS='YES'
IFPRNM='YES'
OSSCAL=1.0
PSSCAL=1.0
NBODY=2
CHKFZC='""" + (fzctemp) + """'
LPRINT=2"""
        ajf_body += """
/

&MP2DNS
/

&MP2GRD
/

&MP3"""
        if self.ajf_method == 'MP3':
            ajf_body += """
  IFSCS='YES'
  IFMANU='YES'
  SCL2OS=1.0D0
  SCL2PS=1.0D0
  SCL3=0.5D0
"""

        ajf_body += """
/
&LMP2
/
"""
        # define new section
        if self.abinit_ver in ['rev15', 'rev17', 'rev22', 'rev23']:
            new_section = """
&LRD
/
"""
        else:
            new_section = ""

        if self.abinit_ver in ['rev10', 'rev11', 'rev15', 'rev17', 'rev22', 'rev23']:
            new_section += """
&DFT
/

&ANALYSIS """
            if self.piedaflag == True:
                new_section += """
PIEDA='YES'"""
            else:
                new_section += """
PIEDA='NO'"""

            new_section += """
/
"""

        elif self.abinit_ver == 'mizuho':
            new_section = """
&DFT
/

&PIEDA"""
            if self.piedaflag == True:
                new_section += """
EnergyDecomposition='YES'"""
            new_section += """
/
"""

        else:
            new_section = ""

        ajf_body += str(new_section) + """
&BSSE"""
        if self.bsseflag == True:
            ajf_body += """
CP='ON'"""
        ajf_body += """
/

&FRAGPAIR
/

"""
        # define solv section
        if self.solv_flag is True:
            solv_section = """
&SOLVATION
EFFECT='ON'
ITRMOD='Normal'
MAXITR=100
THICNV=""" + str(self.pbcnv) + """
IGSCC='Core'
INIEHF='OFF'
PRBRAD=1.4
EPSOUT=80.0
EPSIN=1.0
SRFTNS=0.0072
SRFOFF=0.0
NSPHER=1000"""
            if self.abinit_ver in ['rev10', 'rev11', 'rev15', 'mizuho']:
                solv_section += "\nSCREEN='ES+NP' \n"
        else:
            solv_section = """
&SOLVATION
EFFECT='OFF' """

        ajf_body += str(solv_section) + """
/

&PBEQ
"""
        if self.solv_flag is True:
            ajf_body += """
MAXITR=1000
JDGCNV='RMS'
THRCNV=1.0E-5
"""
            if self.abinit_ver in ['rev17', 'rev22', 'rev23']:
                ajf_body +="ATMRAD='" + self.pbmolrad + "'"
            else:
                ajf_body +="MOLRAD='" + self.pbmolrad + "'"
        ajf_body +="""
/

&POP
"""
        if self.resp == True:
            ajf_body += "ESPFIT='ON'\n"
            ajf_body += "ESPTYP='RESP'\n"
        if self.nbo == True:
            ajf_body +="NBOANL='ON'"
        ajf_body +="""
/

&GRIDCNTRL
GRID='NO'
/

&MCP
/

&CAFI
METLOC='PIPE'
IFLOC='OCC'
CHKFZC='NO'
LPRINT=2
/
"""
        # new section2
        if self.abinit_ver in ['rev11', 'rev15', 'rev17', 'rev22', 'rev23']:
            new_section2 = """
&CCPT
/
"""
        else:
            new_section2 = ""

        ajf_body += str(new_section2) + """
&XYZ
/

&FRAGMENT
""" + str(ajf_fragment) + """
/

&MDCNTRL
MD='OFF'
/

&VEL
/

&NHC
/

&TYPEFRAG
/"""
        return ajf_body


    def getfraginfo(self, seg1_conf, seg2_conf):
        frag1_num = len(seg1_conf['seg_info'])
        frag2_num = len(seg2_conf['seg_info'])

        fcount = 1
        frag = [[], []]
        for i in range(frag1_num):
            frag[0].append(fcount)
            fcount += 1
    #         print frag[0][i]

        for i in range(frag2_num):
            frag[1].append(fcount)
            fcount += 1
    #         print frag[1][i]

        return frag

    def getifie(self, target_dir, frag):
        ielist = []
        if self.pbflag is False:
            # print "**** 1.get ifie running*****"
            for i in range(1, self.total_num+1):
                target = target_dir + '/%05d' % i + ".out"
            #    if i % 1000 == 0:
            #        print target + " end"
                energy = self.read_ifie(target)
                ifiesum = self.getifiesum(energy, frag)
                ielist.append(ifiesum)
            # print energies
            out = target_dir + "/ielist_ifie" + self.mp2temp + self.solvtype + "_detail"
            f = open(out, 'w')
            # print ielist
            print("HF-IFIE", "MP2-IFIE", "PB-term", "PB-polar", "PB-nonpolar", file=f)
            for i in range(len(ielist)):
                print(ielist[i][0], ielist[i][1], 0.0, 0.0, 0.0, file=f)
            f.close()

            out = target_dir + "/ielist_ifie" + self.mp2temp + self.solvtype
            f = open(out, 'w')
            # print ielist
            for i in range(len(ielist)):
                print(ielist[i][0] + ielist[i][1], file=f)
            f.close()

        if self.pbflag is True:
            # print("*** 1.get ifie-pb running ***")
            for i in range(1, self.total_num+1):
                target = target_dir + '/%05d' % i + ".out"
                # if i % 1000 == 0:
                    # print(target)
                energy = self.read_ifiepb(target)
                ifiesum = self.getifiesum(energy[0], frag)
                pbsum = self.getifiesumpb(energy[1], frag)
                ielist.append([ifiesum, pbsum])
            # print ielist
            out = target_dir + "/ielist_ifie" + self.mp2temp + self.solvtype + "_detail"
            f = open(out, 'w')
            # print ielist
            print("HF-IFIE", "MP2-IFIE", "PB-term", "PB-polar", "PB-nonpolar", file=f)
            for i in range(len(ielist)):
                print(ielist[i][0][0], ielist[i][0][1], \
                        ielist[i][1][0] + ielist[i][1][1], \
                        ielist[i][1][0], ielist[i][1][1], file=f)
            f.close()

            out = target_dir + "/ielist_ifie" + self.mp2temp + self.solvtype
            f = open(out, 'w')
            # print ielist
            for i in range(len(ielist)):
                print(ielist[i][0][0] + ielist[i][0][1] + \
                        ielist[i][1][0] + ielist[i][1][1], file=f)
            f.close()

    def getifiesum(self, energy, frag):
        if len(energy) == 0:
            return [0, 0]
        ifiesum = []
        hf = 0.0
        mp2 = 0.0
        total_frag = len(frag[0]) + len(frag[1])
        pairnum = int((total_frag * (total_frag - 1)) / 2)
        for i in range(pairnum):
            for j in(len(frag[0]) + 1, total_frag + 1):  # seg2の要素を動かすループ
                for k in(1, len(frag[0]) + 1):  # seg1の要素を動かすループ
                    if int(energy[i][1]) == int(k) and int(energy[i][0]) == int(j):
                        hf += float(energy[i][4])
                        if self.PR_flag is True:
                            mp2 += float(energy[i][6])
                        else:
                            mp2 += float(energy[i][5])
        return [hf * 627.5095, mp2 * self.mp2fac * 627.5095]

    def getifiesumpb(self, energy, frag):
        # print energy
        if len(energy) == 0:
            return [0, 0]
        ifiesum = []
        pb_es = 0.0
        pb_np = 0.0
        total_frag = len(frag[0]) + len(frag[1])
        pairnum = int((total_frag * (total_frag - 1)) / 2)
        for i in range(pairnum):
            for j in(len(frag[0]) + 1, total_frag + 1):  # seg2の要素を動かすループ
                for k in(1, len(frag[0]) + 1):  # seg1の要素を動かすループ
                    if int(energy[i][2]) == int(k) and int(energy[i][1]) == int(j):
                        pb_es += float(energy[i][4])
                        pb_np += float(energy[i][5])
        return [pb_es * 627.5095, pb_np * 627.5095]

    def read_ifie(self, fname):
        ifie = []
        count = 0
        try:
            f = open(fname, "r")
            text = f.readlines()
            f.close()
        except:
            print("can't open", fname)
            return ifie
        flag = False
        # print text
        for i in range(len(text)):
            itemList = text[i][:-1].split()
            # print itemList
            if len(itemList) < 2:
                continue
            if itemList[1] == 'MP2-IFIE':
                flag = True
                # head.append(itemList)
                continue
            if itemList[1] == 'Mulliken' or itemList[1] =='PIEDA':
                # flag = False
                break
            if flag is True:
                count += 1
            if flag is True and count > 2:
                ifie.append(itemList)

        if flag is False:
            try:
                print("can't read ifie", fname.split("/")[1])
            except:
                pass

        for i in range(len(ifie)):
            if float(ifie[i][4]) < -2 or float(ifie[i][5]) < -2:
                ifie[i][4] = 0.0
                ifie[i][5] = 0.0
                ifie[i][6] = 0.0
        # print ifie
        return ifie

    def read_ifiepb(self, fname):
        ifie = []
        pbterm = []
        count = 0
        count2 = 0
        hit = 0
        try:
            f = open(fname, "r")
            text = f.readlines()
            f.close()
        except:
            print("can't open", fname)
            return [ifie, pbterm]
        flag = False
        flag2 = False
        # print text
        for i in range(len(text)):
            itemList = text[i][:-1].split()
            # print itemList
            if len(itemList) < 2:
                continue
            if itemList[1] == 'MP2-IFIE':
                hit += 1
                if hit == 2:
                    flag = True
                # head.append(itemList)
                continue
            if flag is True and itemList[1] == 'Mulliken':
                flag = False
                # break
            if flag is True:
                count += 1
            if flag is True and count > 2:
                ifie.append(itemList)

            if len(itemList) >= 4:
                if itemList[1] == 'SUMMARY' and itemList[3] == 'ELECTROSTATIC':
                    flag2 = True
                    continue
            if len(itemList) >= 4:
                if itemList[1] == 'ESTIMATION' and itemList[3] == 'NON-POLAR':
                    flag2 = False
                    break
            if flag2 is True:
                count2 += 1
            if flag2 is True and count2 >= 4:
                pbterm.append(itemList)
        if hit != 2 and flag2 is False:
            print("can't read pb", fname.split("/")[1])

        for i in range(len(ifie)):
            if float(ifie[i][4]) < -2 or float(ifie[i][5]) < -2:
                ifie[i][4] = 0.0
                ifie[i][5] = 0.0
                ifie[i][6] = 0.0
                pbterm[i][4] = 0.0
                pbterm[i][5] = 0.0

        return [ifie, pbterm]

    def getTE(self, target_dir, molname, mode,  fzcflag):
        elistname = target_dir + "/energylist_" + self.solvtype
        if mode == "batch":
            energies = []
            if self.pbflag is False:
                # print("**** 1.capt te running*****")
                for i in range(1, self.total_num+1):
                    target = target_dir + '/%05d' % i + ".out"
                    # if i % 1000 == 0:
                        # print(target + " end")
                    energy = self.captfmomp2e(target)
                    energies.append(energy)
                # print energies
                f = open(elistname, "w")
                for i in energies:
                    try:
                        print(i[0], i[1], file=f)
                    except TypeError:
                        print(0, 0, file=f)
                f.close()
                # print("create", elistname)

            if self.pbflag is True:
                # print("*** 1.capt bepb running ***")
                for i in range(1, self.total_num+1):
                    target = target_dir + '/%05d' % i + ".out"
                    # if i % 1000 == 0:
                        # print(target)
                    # print(target)
                    energy = self.getfmopbenergy(target)
                    energies.append(energy)
                # print energies
                f = open(elistname, "w")
                for i in energies:
                    try:
                        print(i[0], i[1], i[2], i[3], i[4], file=f)
                    except TypeError:
                        print(0, 0, 0, 0, 0, file=f)

                f.close()
                # print("*****create", elistname)

        if mode == "single":
            target = target_dir + "/" + molname + ".out"
            if fzcflag == True:
                # print (molname + ": nofzc")
                target = target_dir + "/nofzc/" + molname + ".out"
            if self.pbflag is False:
                self.captmom_single(target)
            elif self.pbflag is True:
                self.captpb_single(target)

        return

    def captmom_single(self, target):
        out = os.path.splitext(target)[0]
        fmoflag = self.getmo_or_fmo(target)  # FMO計算かMO計算かの判定
        # print(target)
        # print(fmoflag)
        if fmoflag is True:
            hf, mp2 = self.captfmomp2e(target)
        else:
            hf, mp2 = self.getmomp2ene(target)
        f = open(out + "_te" + self.solvtype + ".dat", "w")
        # print("hf:", hf, "mp2:", mp2)
        print(hf, mp2, file=f)
        f.close()

    def captpb_single(self, target):
        out = os.path.splitext(target)[0]
        fmoflag = self.getmo_or_fmo(target)  # FMO計算かMO計算かの判定
        if fmoflag is True:
            eg, cor, dg, dges, dgnp = self.getfmopbenergy(target)
        else:
            eg, cor, dg, dges, dgnp = self.getmopbenergy(target)
        f = open(out + "_te" + self.solvtype + ".dat", "w")
        # print("eg:" + eg, "cor:" + cor, "dg:" + dg, "dges:" + dges, "dgnp:" + dgnp)
        print(eg, cor, dg, dges, dgnp, file=f)
        f.close()

    def getmo_or_fmo(self, target):
        f = open(target, "r")
        text = f.readlines()
        f.close()
        fmoflag = False
        for i in range(len(text)):
            itemList = text[i][:-1].split()
            # print itemList
            if "FMO" in itemList and "ON" in itemList:
                fmoflag = True  # FMO mom

                # print(target + " is FMO")
                break
        return fmoflag

    def captfmomp2e(self, target):
        mp2 = 0
        hf = 0
        try:
            f = open(target, "r")

            text = f.readlines()
            f.close()
            for i in range(len(text)):
                itemList = text[i][:-1].split()
                if itemList == ['##', 'FMO', 'TOTAL', 'ENERGY']:
                    hf = text[i + 6].split()
                    mp2 = text[i + 8].split()

                    return [eval(hf[3]), eval(mp2[2])]
        except:
            print("Warning: can't get result:", target)
            return [0, 0]
            # sys.exit()

    def getmomp2ene(self, target):
        mp2 = 0
        hf = 0
        try:
            f = open(target, "r")
            text = f.readlines()
            f.close()
        except:
            print("Error! can't open monomer file:", target)
            return [0,0]
            # sys.exit()

        hfflag = False
        for i in range(len(text)):
            itemList = text[i][:-1].split()
            # print itemList
            if itemList == ['SCF', 'COMPLETED']:
                # print itemList
                if hfflag is False:
                    hf = text[i + 4].split()
                    hfflag = True
            if itemList == ['##', 'MP2', 'ENERGY']:
                # print itemList
                mp2 = text[i + 2].split()
                break

        try:
            return [eval(hf[3]), eval(mp2[2])]
        except:
            print("Error! can't get monomer result:", target)
            return [0,0]
            # sys.exit()



    def getfmopbenergy(self, target):
        try:
            f = open(target, "r")
        except:
            print ("can't open " + target)
            return 0,0,0,0,0
        if self.abinit_ver in ['rev11', 'rev15', 'mizuho']:
            text = f.readlines()
            f.close()
            index = 0
            pbdone =False
            for i in range(len(text)):
                itemList = text[i][:-1].split()
                if len(itemList) == 0:
                    continue
                if itemList[0:4] == ['Energy', 'in', 'gas', 'phase']:
                    eg = itemList  # [5]
                if itemList[0:4] == ['ElectroStatic', '[', 'DGes=Es-Eg-DE', ']']:
                    dges = itemList
                if itemList[0:4] == ['Non-polar', '[', 'DGnp', ']']:
                    dgnp = itemList
                if itemList[0:4] == ['Total', '[', 'DG=DGes+DGnp', ']']:
                    dg = itemList  # [5]
                    pbdone =True
                if itemList == ['##', 'FMO', 'TOTAL', 'ENERGY']:
                    index = i
            if  pbdone is False:
                return 0, 0, 0, 0, 0
            cor = text[index+8].split()
            # print eg, cor, dg
            try:
                return eg[8], cor[2], dg[5], dges[5], dgnp[5]
            except:
                print("Warning: can't get result:", target)
                return 0, 0, 0, 0, 0

        if self.abinit_ver in ['rev17', 'rev22', 'rev23']:
            text = f.readlines()
            f.close()
            index = 0
            pbdone =False
            getflag = False
            for i in range(len(text)):
                itemList = text[i][:-1].split()
                if len(itemList) < 3:
                    continue
                if itemList[0:4] == ['##', 'SOLUTE', 'TOTAL', 'ENERGY']:
                    getflag = True
                if  getflag == True:
                    if itemList[0:3] == ['FMO2', 'in', 'vacuo']:
                        # print('vacuo')
                        eg = itemList  # [5]
                    if itemList[0:2] == ['in', 'solvent']:
                        # print('in solv')
                        esol = itemList
                    if itemList[0] == 'difference':
                        # print('diff')
                        dif = itemList
                    if itemList[0] == 'Electrostatic':
                        # print('ES')
                        dges = itemList
                    if itemList[0] == 'Nonpolar':
                        # print('NP')
                        dgnp = itemList
                    if itemList[0] == 'Total':
                        # print('total')
                        dg = itemList
                        pbdone =True
                if pbdone:
                    break
            # egas, cor(insolv), dgtotal, dges, dgnp
            # print(eg[3], esol[3], dg[2], dges[2], dgnp[2])
            try:
                return eg[3], esol[3], dg[2], dges[2], dgnp[2]
            except:
                print("Warning: can't get result:", target)
                return 0, 0, 0, 0, 0

            '''
            ## SOLUTE TOTAL ENERGY

                                       SCF / Hartree   MP2 corr. / Hartree       Total / Hartree
               FMO2 in vacuo           -236.52210938           -0.80641280         -237.32852218
                    in solvent         -236.52238631           -0.80642693         -237.32881324
                    difference           -0.00027693           -0.00001413           -0.00029106
                  ( ES correction        -0.00028778 )


             ## SOLVATION FREE ENERGY

                                       FMO2 / Hartree       FMO2 / kcal/mol
                Electrostatic             -0.04769392          -29.92838623
                Nonpolar                   0.00548501            3.44189687
                Total                     -0.04220891          -26.48648936
            '''


    def getmopbenergy(self, target):
        count = 0
        try:
            f = open(target, "r")
            text = f.readlines()
            f.close()
        except:
            print ("can't open " + target)
            return 0, 0, 0, 0, 0

        if self.abinit_ver in ['rev11', 'rev15', 'mizuho']:

            for i in range(len(text)):
                itemList = text[i][:-1].split()
        #         print itemList
                if len(itemList) == 0:
                    continue
                if itemList[0:4] == ['Energy', 'in', 'gas', 'phase']:
                    eg = itemList  # [5]
                if itemList[0:4] == ['ElectroStatic', '[', 'DGes=Es-Eg-DE', ']']:
                    dges = itemList
                if itemList[0:4] == ['Non-polar', '[', 'DGnp', ']']:
                    dgnp = itemList
                if itemList[0:4] == ['Total', '[', 'DG=DGes+DGnp', ']']:
                    dg = itemList  # [5]
                if itemList == ['##', 'MP2', 'ENERGY']:
                    count += 1
                    if count == 2:
                        cor = text[i+2].split()
            # print eg, cor, dg
            try:
                return eg[8], cor[2], dg[5], dges[5], dgnp[5]
            except:
                print("Error! can't get monomer result:", target)
                # sys.exit()
                return 0, 0, 0, 0, 0

        if self.abinit_ver in ['rev17', 'rev22', 'rev23']:
            a, b, c, d, e = self.getfmopbenergy(target)
            return a, b, c, d, e



    def readajf(self, fname):
        lines = open(fname, 'r').readlines()

        flag = False
        baseline = 1

        datas = []
        fatomnums = []
        fchgs = []
        fbaas = []
        nfcount = 0
        fatminfo = []
        fatminfos = []
        connects = []
        for i in range(len(lines)):
            itemlist = lines[i].split()
            if len(itemlist) == 0:
                continue
            if itemlist[0][:2].upper() == 'NF':
                self.nf = int(itemlist[0].split('=')[-1])
                nf = self.nf
                print('NF =', self.nf)
                if self.nf > 10:
                    baseline = math.ceil(self.nf/10)
                    print('n_line', baseline)
            if itemlist[0] == '&FRAGMENT':
                flag = True
                typcount = 0
                lcount = 0
                continue
            # fatom section
            if flag == True and typcount == 0:
                fatomnums.append(itemlist)
                lcount += 1
                if lcount == baseline:
                    typcount += 1
                    lcount = 0
                    fatomnums = self.flatten(fatomnums)
                    continue
            # chg section
            if flag == True and typcount == 1:
                fchgs.append(itemlist)
                lcount += 1
                if lcount == baseline:
                    typcount += 1
                    lcount = 0
                    fchgs = self.flatten(fchgs)
                    continue

            # bda section
            if flag == True and typcount == 2:
                fbaas.append(itemlist)
                lcount += 1
                if lcount == baseline:
                    typcount += 1
                    lcount = 0
                    fbaas = self.flatten(fbaas)
                    continue

            # seginfo section
            if flag == True and typcount == 3:
                # print('typ', typcount)
                # print(nfcount)
                fatomnum = fatomnums[nfcount]
                base2 = math.ceil(int(fatomnum)/10)
                fatminfo.append(itemlist)
                lcount += 1
                if lcount == base2:
                    nfcount += 1
                    lcount = 0
                    fatminfo = self.flatten(fatminfo)
                    fatminfos.append(fatminfo)
                    fatminfo = []
                    if  nfcount > nf - 1:
                        typcount += 1
                        lcount = 0
                        continue

            # bda-baa atom section
            if flag == True and itemlist[0] == '/':
                break

            if flag == True and typcount == 4:
                connects.append(itemlist)
                connects = self.functor(int, connects)
                # datas.append(itemlist)

        self.fatomnums = fatomnums
        self.fchgs = fchgs
        self.fbaas = fbaas
        self.fatminfos = fatminfos
        self.connects = connects
        return self
        # return fatomnums, fchgs, fbaas, fatminfos, connects


    def modifyfragparam(self, totalMol, atomnameMol, molnames, anummols, posMol, heads, labs, chains ,resnums ,codes ,occs ,temps ,amarks ,charges, fatomnums, fchgs, fbaas, fatminfos, connects, bridgeds, doubles, nagatoms, nagmolids, nagbdas):
        # modify start
        for i in range(len(nagatoms)):
            for j in range(len(fatminfos)):
                # i: nag mol index
                # j: ajf fragment index

                # print(fatminfos[j])
                # print(nagatoms[i][0])
                dendflag = False
                dmidflag = False
                if nagatoms[i][0] in fatminfos[j]: # search (asn + nag)
                    print('NAG', i, 'frag', j, 'start atom', nagatoms[i][0])

                    # check baa asn atom id
                    tgtid = j
                    for bridged in bridgeds:
                        if tgtid >= bridged:
                            tgtid += 1
                    print(molnames[tgtid], resnums[tgtid][0])
                    print('atomnameMol', tgtid, atomnameMol[tgtid])
                    baaidx = atomnameMol[tgtid].index(' ND2')
                    print('asn_baaidx', anummols[tgtid][baaidx])


                    tgtbaaid = sum(fbaas[:j+1])
                    print('tgtbaaid', tgtbaaid)

                    #fatomnums
                    # j: asn frag id
                    fatomnums[j] -= len(nagatoms[i])
                    fatomnums.append(len(nagatoms[i]))

                    print('doubles', doubles)
                    for double in doubles:
                        if double[1] == nagmolids[i]:
                            dendflag = True
                        if double[0] == nagmolids[i]:
                            dmidflag = True

                    if dendflag == True:
                        baaidx = atomnameMol[nagmolids[i-1]].index(' O4 ')
                        print('baaidx double', baaidx, nagmolids[i-1])
                        baaatomid = anummols[nagmolids[i-1]][baaidx]
                        print('baaaaaaaaaaaaaa', baaatomid)

                    # for old(asn) frag
                    if dendflag == False:
                        fchgs[j] -= 1
                        fbaas[j] += 1

                    # for new frag
                    if dmidflag == True:
                        fchgs.append(0)
                        fbaas.append(1)
                    else:
                        fchgs.append(1)
                        fbaas.append(0)

                    #fatminfos
                    delid = []

                    # get delete atom index in mol j
                    for k in range(len(nagatoms[i])):
                        # k: nagatom(in mol) index
                        for l in range(len(fatminfos[j])):
                            # l: atom(in fragment) index
                            if nagatoms[i][k] == fatminfos[j][l]:
                                delid.append(l)
                                break

                    print('delid', delid)
                    dellist = lambda items, indexes: [item for index, item in enumerate(items) if index not in indexes]
                    fatminfos[j] = dellist(fatminfos[j], delid)
                    # atomnameMol[tgtid] = dellist(atomnameMol[j], delid)
                    print('atom-len-check frag', j, ':',  len(fatminfos[j]), '-- mol', j, ':',  len(atomnameMol[tgtid]))
                    fatminfos.append(nagatoms[i])

                    #connects
                    print(len(connects))
                    if dendflag == True:
                        connects.append([int(nagbdas[i]), int(baaatomid)])
                        print('connects', nagbdas[i], baaatomid)

                    else:
                        connects.insert(tgtbaaid, [int(nagbdas[i]), int(anummols[tgtid][baaidx])])
                        print('connects', nagbdas[i], anummols[tgtid][baaidx])
                    print(' ')
        return fatomnums, fchgs, fbaas, connects, fatminfos

    def flatten(self, nested_list):
        """2重のリストをフラットにする関数"""
        return [int(e) for inner_list in nested_list for e in inner_list]

    # def listtoint(nested_list):
    #     """リストの中をintegerに"""
    #     return [int(e) for inner_list in nested_list for e in inner_list]

    def functor(self, f, l):
        if isinstance(l,list):
            return [self.functor(f,i) for i in l]
        else:
            return f(l)


    def getfragdict(self, fnames, ofile):
        f = open(ofile, 'w')

        print("seg_data = [", file=f)
        for fname in fnames:
            head, ext = os.path.splitext(fname)
            # fatomnums, fchgs, fbaas, fatminfos, connects, = self.readajf(fname)
            self = self.readajf(fname)
            print("    {", file=f)
            print("    'name': '" + head.split('/')[-1] + "',", file=f)
            print("    'atom': ", self.fatomnums, ',', file=f)
            print("    'charge': ", self.fchgs, ',', file=f)
            print("    'connect_num': ", self.fbaas, ',', file=f)
            print("    'seg_info': ", self.fatminfos, ',', file=f)
            print("    'connect': ", self.connects, ',', file=f)
            end =  """    'nummol_seg': [1],
    'repeat': [1],
    'pair_file': [],
    'multi_xyz': 'none'
    },"""
            print(end, file=f)

        print(']', file=f)

        '''
                seg_conf = {'name': seg_name,
                            'atom': [seg_atom],
                            'charge': [0],
                            'connect_num': [0],
                            'connect': [],
                            'seg_info': [[i + 1 for i in range(seg_atom)]],
                            'nummol_seg': [1],
                            'repeat': [1],
                            'pair_file': [],
                            'multi_xyz': 'none'}
        '''

    def config_read(self, seg_name, seg_atom=0):
        fdata = {}
        exec(open(self.mainpath + "/segment_data.dat", "r").read(), fdata)
        seg_data = fdata['seg_data']
        # print (seg_data)
        # default data
        seg_conf = {'name': seg_name,
                    'atom': [seg_atom],
                    'charge': [0],
                    'connect_num': [0],
                    'connect': [],
                    'seg_info': [[i + 1 for i in range(seg_atom)]],
                    'nummol_seg': [1],
                    'repeat': [1],
                    'pair_file': [],
                    'multi_xyz': 'none'}

        for i in range(len(seg_data)):
            if seg_name == seg_data[i]['name']:
                for key, data in seg_data[i].items():
                    seg_conf[key] = data

            # print "seg_conf =", seg_conf
        return seg_conf

    def getmb_frag_seclists(self, param, nameid):
        frag_atom = param[0]
        frag_charge = param[1]
        frag_connect_num = param[2]
        frag_connect = param[3]
        seg_info = param[4]

        # fragment atom
        count = 0

        frag_atoms = []
        for i in range(len(nameid)):
            for j in range(len(frag_atom[nameid[i]])):
                    frag_atoms.append(frag_atom[nameid[i]][j])
                    count += 1

        # fragment charge
        count = 0
        frag_charges = []
        for i in range(len(nameid)):
            for j in range(len(frag_charge[nameid[i]])):
                frag_charges.append(frag_charge[nameid[i]][j])
                count += 1

        # fragment connect num
        count = 0
        frag_baanums = []
        for i in range(len(nameid)):
            for j in range(len(frag_connect_num[nameid[i]])):
                frag_baanums.append(frag_connect_num[nameid[i]][j])
                count += 1

        # fragment body
        atom_count = 0
        # i mol  j frag k atom
        frag_atmlabss = []
        for i in range(len(nameid)):
            for j in range(len(seg_info[nameid[i]])):
                icount = 0
                frag_atmlabs = []
                for k in range(len(seg_info[nameid[i]][j])):
                    icount += 1
                    frag_atmlabs.append(seg_info[nameid[i]][j][k] + atom_count)
                frag_atmlabss.append(frag_atmlabs)
            atom_count += sum(frag_atom[nameid[i]])

        atom_count = 0
        # connect info
        frag_connects = []
        frag_connectss = []
        for i in range(len(nameid)):
            for j in range(len(frag_connect[nameid[i]])):
                frag_connects = []
                frag_connects.append(frag_connect[nameid[i]][j][0] + atom_count)
                frag_connects.append(frag_connect[nameid[i]][j][1] + atom_count)
                frag_connectss.append(frag_connects)
            atom_count += sum(frag_atom[nameid[i]])

        self.fatomnums = frag_atoms
        self.fchgs = frag_charges
        self.fbaas = frag_baanums
        self.fatminfos = frag_atmlabss
        self.connects = frag_connectss
        return  self

        # return frag_atoms, frag_charges, frag_baanums, frag_atmlabss, frag_connectss


    def writemb_frag_section(self, param, nameid):
        frag_atom = param[0]
        frag_charge = param[1]
        frag_connect_num = param[2]
        frag_connect = param[3]
        seg_info = param[4]

        ajf_fragment = ''
        # fragment atom
        count = 0
        for i in range(len(nameid)):
            for j in range(len(frag_atom[nameid[i]])):
                    ajf_fragment += '%8d' % (frag_atom[nameid[i]][j])
                    count += 1
                    if count % 10 == 0:
                        ajf_fragment += '\n'
        if count % 10 != 0:
            ajf_fragment += '\n'

        # fragment charge
        count = 0
        for i in range(len(nameid)):
            for j in range(len(frag_charge[nameid[i]])):
                ajf_fragment += '%8d' % (frag_charge[nameid[i]][j])
                count += 1
                if count % 10 == 0:
                    ajf_fragment += '\n'
        if count % 10 != 0:
            ajf_fragment += '\n'

        # fragment connect num
        count = 0
        for i in range(len(nameid)):
            for j in range(len(frag_connect_num[nameid[i]])):
                ajf_fragment += '%8d' % (frag_connect_num[nameid[i]][j])
                count += 1
                if count % 10 == 0:
                    ajf_fragment += '\n'
        if count % 10 != 0:
            ajf_fragment += '\n'

        # fragment body
        atom_count = 0
        for i in range(len(nameid)):
            for j in range(len(seg_info[nameid[i]])):
                icount = 0
                for k in range(len(seg_info[nameid[i]][j])):
                    icount += 1
                    ajf_fragment += '%8d' % (
                            seg_info[nameid[i]][j][k] + atom_count)
                    if icount % 10 == 0:
                        icount = 0
                        ajf_fragment += '\n'
                if icount % 10 != 0:
                    ajf_fragment += '\n'
            atom_count += sum(frag_atom[nameid[i]])

        print ('atom_count', atom_count)
        atom_count = 0
        # connect info
        print("frag_atom", frag_atom)
        for i in range(len(nameid)):
            for j in range(len(frag_connect[nameid[i]])):
                ajf_fragment += '%8d' % (
                        frag_connect[nameid[i]][j][0] + atom_count)
                ajf_fragment += '%8d' % (
                        frag_connect[nameid[i]][j][1] + atom_count)
                ajf_fragment += '\n'
            atom_count += sum(frag_atom[nameid[i]])

        return ajf_fragment


    def saveajf(self, oname=None):
        if oname == None:
            oname = os.path.splitext(self.readgeom)[0] + '-' + self.ajf_method + '-' + self.ajf_basis_set.replace('*', 'd') + '.ajf'
            ohead = os.path.splitext(self.readgeom)[0] + '-' + self.ajf_method
        else:
            ohead = os.path.splitext(oname)[0]
        ajf_body = self.gen_ajf_bodywrap(ohead)

        print('ajf_oname:', oname)
        ajf_file = open(oname, 'w')

        print(ajf_body, file=ajf_file)

