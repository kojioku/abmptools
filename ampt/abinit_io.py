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


class abinit_io(mi.mol_io):

    def __init__(self):
        super().__init__()
        # gen_rand: gen_coord
        self.solv_flag = False

        # submit param
        self.submit_system = 'none'
        self.npro = 4
        self.nnode = 0
        self.para_job = 2
        self.memory = "1800"
        self.queue = None

        # abinitmp param
        self.ajf_method = "MP2"
        self.ajf_basis_set = "6-31G*"
        self.abinit_ver = 'rev10'
        self.abinitmp_path = 'abinitmp'
        self.mpipath = 'mpirun'
        self.pbmolrad = 'vdw'
        self.pbcnv = 1.0
        self.piedaflag = False
        # distlitname

        # others
        self.python_path = 'python'
        self.pynp = 4
        self.babelflag = False
        self.targets = []
        self.mainpath = '.'
        self.nofzc_flag = False
        self.fzcnum = 0
        self.ksubcmd = 'pjsub'

        return


    def get_fragsection(self, param):
        frag_atom, frag_charge, frag_connect_num, frag_connect, seg_info = param

        ajf_fragment = ''
        # fragment atom
        for i in range(len(frag_atom)):
            icount = 0
            for j in range(len(frag_atom[i])):
                icount += 1
                ajf_fragment += '%8d' % (frag_atom[i][j])
                if icount % 10 is 0:
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
                if icount % 10 is 0:
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
                if icount % 10 is 0:
                    icount = 0
                    ajf_fragment += '\n'
            if icount % 10 != 0:
                ajf_fragment += '\n'

        # fragment body
        atom_count = 0
        for i in range(len(seg_info)):
            for j in range(len(seg_info[i])):
                icount = 0
                for k in range(len(seg_info[i][j])):
                    icount += 1
                    ajf_fragment += '%8d' % (seg_info[i][j][k] + atom_count)
                    if icount % 10 is 0:
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

    def gen_ajf_body(self, param, fzctemp="YES"):
        ajf_charge, read_geom, cpf, ajf_fragment, num_fragment, MOME_flag = param

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
            if self.abinit_ver != 'rev17':
                solv_section += """
SCREEN='ES+NP' """
        else:
            solv_section = """
&SOLVATION
EFFECT='OFF' """

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

        if self.abinit_ver == 'rev15' or self.abinit_ver == 'rev17':
            new_section = """
&LRD
/
"""
        else:
            new_section = ""
        if self.abinit_ver == 'rev10' or self.abinit_ver == 'rev11' or self.abinit_ver == 'rev15' or self.abinit_ver =='rev17':
            new_section += """
&DFT
/

&ANALYSIS """
            if self.piedaflag == True:
                new_section += """
PIEDA='YES'"""
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

        if self.abinit_ver == 'rev11' or self.abinit_ver == 'rev15' or self.abinit_ver == 'rev17':
            new_section2 = """
&CCPT
/
"""
        else:
            new_section2 = ""

        if self.submit_system == 'K':
            np = '1'
        else:
            np = str(int(self.npro/self.para_job)) 

        ajf_body = """&CNTRL
ElecState='S1'
Method='""" + str(self.ajf_method) + """'
Nprint=3
Memory=""" + str(self.memory) + """
Natom=0
Charge=""" + str(ajf_charge) + """
ReadGeom=""" + str(read_geom) + "\n"

        if self.cpfflag == True:
            ajf_body += "WriteGeom=" + str(cpf) + "\n"
        ajf_body += """Gradient='NO'
Vector='OFF'
CPFBIN='NO'
THOVL=1.0E-12
E_THSWZ=1.0E-12
G_THSWZ=1.0E-12
/

&FMOCNTRL
FMO='""" + str(use_FMO) + """'
NBody=2
AutoFrag='OFF' """ + str(frag_section) + """
Laoc=0.0
Lptc=2.0
esp_ptc_multipole='NO'
Ldimer=2.0
Dimer_es_multipole='NO'
NP=""" + np + """
MaxSCCcyc=250
MaxSCCenergy=5.0E-7
/

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

&SCZV
DimerResponseTerm='NO'
/

&MP2
NP_MP2_IJ=1
NP_MP2_S=0
MemoryMP2=0
IFSCS='YES'
IFPRNM='YES'
OSSCAL=1.0
PSSCAL=1.0
NBODY=2
CHKFZC='""" + (fzctemp) + """'
LPRINT=2
/

&MP2DNS
/

&MP2GRD
/

&MP3
/

&LMP2
/
""" + str(new_section) + """
&BSSE
/

&FRAGPAIR
/
""" + str(solv_section) + """
/

&PBEQ
MAXITR=1000
JDGCNV='RMS'
THRCNV=1.0E-5"""
        if self.abinit_ver == 'rev17':
            ajf_body +="""
ATMRAD='""" + self.pbmolrad + """'"""
        else:
            ajf_body +="""
MOLRAD='""" + self.pbmolrad + """'"""

        ajf_body +="""
/

&POP
ESPFIT='ON'
ESPTYP='RESP'
NBOANL='ON'
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
""" + str(new_section2) + """
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



    def ajf_make_monomer(self, se1_conf, seg1_name, MOME_flag, nofzc=False):
        frag_atom = [se1_conf['atom']]
        frag_charge = [se1_conf['charge']]
        frag_connect_num = [se1_conf['connect_num']]
        frag_connect = [se1_conf['connect']]
        seg_info = [se1_conf['seg_info']]

        # fragment body
        ajf_fragment = self.get_fragsection([frag_atom, frag_charge, frag_connect_num, frag_connect, seg_info])
        ajf_fragment = ajf_fragment[:-1]

        # fragment section
        ajf_charge = sum(frag_charge[0])
        num_fragment = len(frag_atom[0])
        ajf_parameter = [ajf_charge, "", "", ajf_fragment, num_fragment, MOME_flag]
        # gen ajf file
        if nofzc == False:
            ajf_file_name = 'monomer/' + seg1_name + ".ajf"
            fzctemp = "YES"
        else:
            ajf_file_name = 'monomer/nofzc/' + seg1_name + ".ajf"
            fzctemp = "NO"

        ajf_file = open(ajf_file_name, 'w', newline = "\n")

        ajf_parameter[1] = "'pdb/" + seg1_name + ".pdb'"
        ajf_parameter[2] = "'" + seg1_name + '-' + \
            self.ajf_method + '-' + self.ajf_basis_set + ".cpf'"
        ajf_body = self.gen_ajf_body(ajf_parameter, fzctemp)
        print(ajf_body, file=ajf_file)
        return

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
            if itemList[1] == 'Mulliken':
                # flag = False
                break
            if flag is True:
                count += 1
            if flag is True and count > 2:
                ifie.append(itemList)

        if flag is False:
            print("can't read ifie", fname.split("/")[1])

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


    def getmopbenergy(self, target):
        count = 0
        try:
            f = open(target, "r")
            text = f.readlines()
            f.close()
        except:
            print ("can't open " + target)
            return 0, 0, 0, 0, 0

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

