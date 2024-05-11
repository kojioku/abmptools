import sys
import os
scrdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scrdir)

from multiprocessing import Pool
import copy
import random
import math
import re
import subprocess
import csv
import pdb_io as pdio
import time
from ctypes import *

try:
    import numpy as np
except:
    pass
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

class anlfmo(pdio.pdb_io):
    def __init__(self):
        super().__init__()
        self.cpfflag = True
        self.solv_flag = False  # True -> in water , False -> in vacuum
        self.abinit_ver = True
        self.memory = 3000
        self.npro = 8
        self.para_job = 1
        self.cutmode = 'sphere'
        self.abinit_ver = 'rev11'
        self.piedaflag = True
        self.molname = []
        self.criteria = []
        self.tgtpos = []
        self.icolumn = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE',
                        'PR-TYPE1', 'GRIMME', 'JUNG', 'HILL']
        self.bicolumn = ['I', 'J', 'DIST', 'DIMER-ES-BSSE', 'HF-BSSE',
                         'MP2-BSSE', 'PR-T1-BSSE', 'GRIMME-BSSE',
                         'JUNG-BSSE', 'HILL-BSSE']
        self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']
        self.pcolumnpb = ['I', 'J', 'ES', 'EX', 'CT-mix', 'Solv(ES)',
                          'DI(MP2)', 'q(I=>J)']
        self.ifdfsumcolumn = ['HF-IFIE', 'MP2-IFIE', 'PR-TYPE1',
                              'GRIMME', 'JUNG', 'HILL', 'ES',
                              'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']
        self.logMethod = 'MP2'

        self.anlmode = 'frag'  # frag, 'mol', 'fraginmol', 'multi'
        self.fragmode = 'auto'  # 'hybrid', 'auto', 'manual'
        self.dist = 1000.0
        self.tgt1frag = None

        self.rpdbflag = False
        self.pdbname = None
        self.is_disp = False

        # -- for mol mode or multi mode--
        self.tgt2type = 'frag'  # frag: mol-frag, mol: mol-mol

        # mol - mol mode
        self.selecttype = 'molid'
        self.tgtmolid = None

        # -- fraginmol mode --
        self.tgt1_lofrag = None
        self.tgt2_lofrag = None
        self.tgt2molname = None
        # ------ multi mode ------
        # if tgt2type == 'frag':
        self.ifdfs = []
        self.pidfs = []
        self.tgt1frag = None
        self.tgt2frag = None

        # if tgt2type == 'molname':
        self.tgt2dist = None

        # multi file setting
        self.ilog_head = None
        self.ilog_tail = None
        self.pdb_head = None
        self.pdb_tail = None

        self.start = None
        self.end = None
        self.interval = None

        # --hybrid mode
        self.hyfrag = None  # 320
        self.hynum = None

        self.pynp = 3
        self.addresinfo = True
        self.writeresnamecsv = True

        mydir = os.path.dirname(os.path.abspath(__file__))
        self.f90sofile = mydir + '/f90/bin/readifiepiedalib.so'
        self.f90soflag = True

        # for svd
        self.matrixtype = 'normal'
        self.exceptfrag = []

        self.pbflag = False
        self.nf = 0
        pass

    def read_fraginfo(self, fname):
        frags = []
        count = 0
        with open(fname, 'r') as file:
            flag = False
            while True:
                # for i in range(len(text)):
                # itemList = text[i][:-1].split()
                itemList = file.readline().strip().split()
                print(itemList)
                if len(itemList) < 2:
                    continue
                if itemList[1] == 'AUTOMATIC' or itemList[1] == 'HYBRID':
                    flag = True
                    continue
                if itemList[1] == 'MANUAL':
                    manflag = True
                if itemList[1] == 'system' or itemList[0] == 'Ions':
                    print('read end')
                    break
                if flag is True:
                    count += 1
                if flag is True and count > 2:
                    if self.fragmode == 'hybrid':
                        frags.append(itemList[3] + itemList[1])
                    elif self.fragmode == 'auto':
                        frags.append(itemList[2] + itemList[0])

        return frags

    def read_pieda(self, fname):
        ifie = []
        count = 0
        text = open(fname, "r").readlines()
        flag = False
        # print text
        for i in range(len(text)):
            itemList = text[i][:-1].split()
            # print itemList
            if len(itemList) < 2:
                continue
            if itemList[1] == 'PIEDA':
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

        return ifie

    def getconnect(self, idxs, molfrags, df, tgtid):
        '''getconnect

        get connect info

        Args:
            idxs (list): index list
            molfrags (list): molfrag list
            df (DataFrame): ifie DataFrame
            tgtid (int): target id
        Returns:
            newfrags (list): new frag list
            molfrags (list): molfrag list
        '''

        neighbors_list = []
        # print(idxs)
        for idx in idxs:
            tgtdf = df[df['I'] == idx]
            # tgtdf = tgtdf.append(df[df['J'] == idx])
            tgtdf = pd.concat([tgtdf, df[df['J'] == idx]])
            tgtdf_zero = tgtdf[tgtdf['DIST'] == 0.]
            # print(tgtdf_zero)
            neighbor_i = [index for index, row in tgtdf_zero.groupby("I")]
            neighbor_j = [index for index, row in tgtdf_zero.groupby("J")]
            neighbors = set(neighbor_i + neighbor_j)
            # print('connect_idx', neighbors)
            neighbors_list.append(neighbors)
        neighbors_flat = list(itertools.chain.from_iterable(neighbors_list))
        # print(idx)

        newfrags = []
        for idx in neighbors_flat:
            if idx == tgtid:
                if self.anlmode != 'fraginmol':
                    continue
            if idx not in molfrags:
                molfrags.append(idx)
                newfrags.append(idx)

        return newfrags, molfrags

    def getmolfrags(self, tgtid, df):
        '''getmolfrags

        get molfrag list

        Args:
            tgtid (int): target id
            df (DataFrame): ifie DataFrame
        Returns:
            molfrags (list): molfrag list
        '''

        molfrags  = [tgtid]
        newfrags = [tgtid]
        while True:
            newfrags, molfrags = self.getconnect(newfrags, molfrags, df, tgtid)
            # print (newfrags, molfrags)
            if len(newfrags) == 0:
                break
        molfrags.sort()
        # print('aaa', molfrags)

        return molfrags

    def getallmolfrags(self, logname, df, nf):
        # fragment connect is judged by checking frag-frag distance
        alfrags = []
        molfragss = []
        for i in range(1, nf+1):
            if i in alfrags:
                # print(i, 'already')
                continue
            molfrags = self.getmolfrags(i, df)
            molfragss.append(molfrags)
            for j in molfrags:
                alfrags.append(j)
        # print(molfragss)
        return molfragss

    def getlognf(self, logname, fragmode):
        # print('fragmode', fragmode)
        if fragmode == 'manual':
            text = open(logname, "r").readlines()
            for i in range(len(text)):
                itemList = text[i][:-1].split()
                # print(itemList)
                if len(itemList) < 2:
                    continue
                if itemList[:2]== ['NF', '=']:
                   nf = int(itemList[2])
                   break
        if fragmode == 'auto':
            f = open(logname, 'r')
            readflag = False
            autoreadflag = False
            fragdata = []
            fragdatas = []
            elecs = []
            seqnos = []
            fragnos = []
            residuestr = []
            logreadGeom = []
            pdbabs = ""

            for line in f:
                items = line[1:].split()
                chains = line[0]
                if len(items) == 0:
                    continue

                if items[0] == 'ReadGeom':
                     logreadGeom = items[2]

                if items[0] == 'AutoFrag':
                    if items[2] == 'ON':
                        self.fragmode = 'auto'
                    else:
                        self.fragmode = 'manual'

                # read frag table
                # if len(items) == 3:
                    # print (items)
                if items[0:3] == ['Frag.', 'Elec.', 'ATOM']:
                    readflag = True
                    # print('# readflag ON #')
                    continue
                if items[0:2] == ["ALL", "ELECTRON"]:
                    fragdatas.append(fragdata)
                    readflag = False
                if items[0:2] == ["ALL", "ATOM"]:
                    natom = int(items[3])
                if readflag is True:
                    if line[0:21] == "                     ":
                        # print(line)
                        fragdata = fragdata + items
                    else:
                        if len(fragdata) != 0:
                            fragdatas.append(fragdata)
                            elecs.append(int(elec))
                        fragdata = []
                        elec = items[1]
                        fragdata = fragdata + items[2:]

                if items [0:2] == ['START', 'FRAGMENT']:
                    break

                ## AUTOMATIC FRAGMENTATION
                if items[0:3] == ['Seq.', 'Frag.', 'Residue']:
                    autoreadflag = True
                    continue

                if autoreadflag == True and items[0:2] == ['The', 'system']:
                    autoreadflag = False
                    continue

                if autoreadflag == True and items[0] == 'Ions':
                    autoreadflag = False
                    continue

                if autoreadflag == True:
                   # print(items)
                   seqnos.append(items[0])
                   fragnos.append(items[1])
                   residuestr.append(items[2])

            nf = len(fragdatas)
            # print('-------nf--------', nf)
            # print(fragdatas)

        return nf

    def readlog(self, logname, fragmode):
        # print('fragmode', fragmode)
        if fragmode == 'manual':
            text = open(logname, "r").readlines()
            for i in range(len(text)):
                itemList = text[i][:-1].split()
                # print(itemList)
                if len(itemList) < 2:
                    continue
                if itemList[:2]== ['NF', '=']:
                   nf = int(itemList[2])
                   break
        if fragmode == 'auto':
            f = open(logname, 'r')
            readflag = False
            autoreadflag = False
            fragdata = []
            fragdatas = []
            elecs = []
            seqnos = []
            fragnos = []
            residuestr = []
            logreadGeom = []
            pdbabs = ""

            for line in f:
                items = line[1:].split()
                chains = line[0]
                if len(items) == 0:
                    continue

                if items[0] == 'ReadGeom':
                     logreadGeom = items[2]

                if items[0] == 'AutoFrag':
                    if items[2] == 'ON':
                        self.fragmode = 'auto'
                    else:
                        self.fragmode = 'manual'

                # read frag table
                # if len(items) == 3:
                    # print (items)
                if items[0:3] == ['Frag.', 'Elec.', 'ATOM']:
                    readflag = True
                    # print('# readflag ON #')
                    continue
                if items[0:2] == ["ALL", "ELECTRON"]:
                    fragdatas.append(fragdata)
                    readflag = False
                if items[0:2] == ["ALL", "ATOM"]:
                    natom = int(items[3])
                if readflag == True:
                    if line[0:21] == "                     ":
                        # print(line)
                        fragdata = fragdata + items
                    else:
                        if len(fragdata) != 0:
                            fragdatas.append(fragdata)
                            elecs.append(int(elec))
                        fragdata = []
                        elec = items[1]
                        fragdata = fragdata + items[2:]

                if items [0:2] == ['START', 'FRAGMENT']:
                    break

                ## AUTOMATIC FRAGMENTATION
                if items[0:3] == ['Seq.', 'Frag.', 'Residue']:
                    autoreadflag = True
                    continue

                if autoreadflag == True and items[0:2] == ['The', 'system']:
                    autoreadflag = False
                    continue

                if autoreadflag == True and items[0] == 'Ions':
                    autoreadflag = False
                    continue

                if autoreadflag == True:
                   # print(items)
                   seqnos.append(items[0])
                   fragnos.append(items[1])
                   residuestr.append(items[2])

            nf = len(fragdatas)
            print('-------nf--------', nf)

        return natom

    def getlognatom(self, fname):
        acount = 0
        flag = False
        f =open(fname, "r")
        text = f.readlines()
        for i in range(len(text)):
            itemList = text[i].split()
            if itemList[0:4] == ['##', 'READ', 'MOLECULAR', 'STRUCTURE']:
                flag = True
                continue
            if flag:
                if itemList[0:3] == ['##', 'Molecular', 'formula']:
                    break
                elif len(itemList) <= 1:
                    continue
                else:
                    acount += 1
        return acount


    def getlogchg(self, fname, natom):
        chgs = []
        f =open(fname, "r")
        text = f.readlines()
        for i in range(len(text)):
            itemList = text[i].split()
            # MO
            if itemList == ['TWO-STAGE', 'RESP', 'FITTING:', 'SECOND', 'STAGE']:
                for j in range(int(natom)):
                    chgval = text[i+20+j].split()
                    chgs.append(float(chgval[2]))

            # FMO
            if itemList == ['##', 'ESP-FITTING', 'TYPE:', 'RESP']:
                for j in range(int(natom)):
                    chgval = text[i+6+j].split()
                    chgs.append(float(chgval[4]))
        return chgs


    def getlogchgall(self, fname, natom, chgtype):
        alabs = []
        elems = []
        ress = []
        frags = []
        chgs = []
        pops = []
        f =open(fname, "r")
        text = f.readlines()
        if chgtype == 'nbo':
            for i in range(len(text)):
                itemList = text[i].split()
                # FMO
                if itemList == ['##', 'NATURAL', 'ATOMIC', 'POPULATIONS']:
                    for j in range(int(natom)):
                        chgval = text[i+6+j].split()
                        alabs.append(int(chgval[0]))
                        elems.append(str(chgval[1]))
                        ress.append(int(chgval[2]))
                        frags.append(int(chgval[3]))
                        chgs.append(float(chgval[4]))
                        pops.append(float(chgval[5]))
                    chgdf=pd.DataFrame({'AtomLabel':alabs,
                                       'Elem': elems,
                                       'Res': ress,
                                       'Frag': frags,
                                       'Charge': chgs,
                                       'Pop': pops})
        else:
            print('Options except A are not supported yet.')
            sys.exit()
        return chgdf

        '''
         ========================================================
           ## NATURAL POPULATION ANALYSIS -- Ver.2.73(20131003)
         ========================================================


           ## NATURAL ATOMIC POPULATIONS

         -----------------------------------------
             Atom  Res Frag     Charge      Pop
                               FMO2-HF    FMO2-HF
         -----------------------------------------
             1 N     1    1  -0.801541   7.801541
             2 C     1    1  -0.108078   6.108078
             3 C     1    2   0.847514   5.152486
             4 O     1    2  -0.776786   8.776786
             5 C     1    1  -0.669410   6.669410
             6 H     1    1   0.461111   0.538889
        '''

    def getlogorpdbfrag(self, ifile):

        f = open(ifile, 'r')
        readflag = False
        autoreadflag = False
        fragdata = []
        fragdatas = []
        elecs = []
        seqnos = []
        fragnos = []
        residuestr = []
        logreadGeom = []
        pdbabs = ""

        for line in f:
            items = line[1:].split()
            chains = line[0]
            if len(items) == 0:
                continue

            if items[0] == 'ReadGeom':
                 logreadGeom = items[2]

            if items[0] == 'AutoFrag':
                if items[2] == 'ON':
                    self.fragmode = 'auto'
                else:
                    self.fragmode = 'manual'

            # read frag table
            if items[0:3] == ['Frag.', 'Elec.', 'ATOM']:
                readflag = True
                continue
            if items[0:2] == ["ALL", "ELECTRON"]:
                fragdatas.append(fragdata)
                readflag = False
            if items[0:2] == ["ALL", "ATOM"]:
                natom = int(items[3])
            if readflag == True:
                if line[0:21] == "                     ":
                    # print(line)
                    fragdata = fragdata + items
                else:
                    if len(fragdata) != 0:
                        fragdatas.append(fragdata)
                        elecs.append(int(elec))
                    fragdata = []
                    elec = items[1]
                    fragdata = fragdata + items[2:]

            if items [0:2] == ['START', 'FRAGMENT'] or items[0] == 'Ions':
                break

            ## AUTOMATIC FRAGMENTATION
            if items[0:3] == ['Seq.', 'Frag.', 'Residue']:
                autoreadflag = True
                continue

            if autoreadflag == True and items[0:2] == ['The', 'system']:
                autoreadflag = False
                continue

            if autoreadflag == True:
               # print(items)
               seqnos.append(items[0])
               fragnos.append(items[1])
               residuestr.append(items[2])


        nf = len(fragdatas)
        # print('natom', natom)
        # print('nf', len(fragdatas))

        na_nfrag = []
        for i in range(len(fragdatas)):
            na_nfrag.append(len(fragdatas[i]))
        # print('na_nfrag\n', na_nfrag)

        # print('Elec.\n', elecs)
        # print('logreadGeom:', logreadGeom)

        # print(pdbabs)
        # print('Frag Atom number\n', fragdatas)

        resname_perfrag = []
        resnamenonum_perfrag = []
        if self.fragmode == 'manual':
            logabsitems = os.path.abspath(ifile).split('/')
            logabsitems[-1] = logreadGeom
            # print(logabsitems)
            pdbabs = ""
            for logabsitem in logabsitems:
                pdbabs = pdbabs + logabsitem + '/'

            pdbabs = pdbabs[:-1]

        ## manual

            self.getpdbinfowrap(pdbabs)
            # print(self.resnameRes)
            tgts = []
            for fragdata in fragdatas:
                try:
                    tgts.append(fragdata[2])
                except:
                    tgts.append(fragdata[0])

            for tgt in tgts:
                for i in range(len(self.gatmlabRes)):
                    gatmlabs = self.gatmlabRes[i]
                    # print(gatmlabs)
                    tgtstr = str(tgt).rjust(5)
                    if tgtstr in gatmlabs:
                        headid = [i, gatmlabs.index(tgtstr)]
                        resname_perfrag.append(self.resnameRes[headid[0]][headid[1]] + self.resnumRes[headid[0]][headid[1]].strip())
                        resnamenonum_perfrag.append(self.resnameRes[headid[0]][headid[1]])


            # print(resname_perfrag)

        ## auto
        # AUTOMATIC FRAGMENTATION
        if self.fragmode == 'auto':
            alreadys = []
            for i in range(len(fragnos)):
                if fragnos[i] in alreadys:
                    continue
                else:
                    resname_perfrag.append(residuestr[i] + seqnos[i])
                    resnamenonum_perfrag.append(residuestr[i])
                    alreadys.append(fragnos[i])
#            print('resname_perfrag', resname_perfrag)

        self.resnamenonum_perfrag = resnamenonum_perfrag
        return resname_perfrag, pdbabs


    def getifiedf(self, ifie, solv=[]):
        '''get ifie data frame from ifie file

        get ifie data frame from ifie file

        Args:
            ifie (str): ifie file name
            solv (list): solvent name list

        Returns:
            df (pandas.DataFrame): ifie data frame
        '''

        df = pd.DataFrame(ifie, columns=self.icolumn)
        df['I'] = df['I'].astype(int)
        df['J'] = df['J'].astype(int)
        df['DIST'] = df['DIST'].astype(float)
        df['HF-IFIE'] = df['HF-IFIE'].astype(float) * 627.5095

        if self.logMethod == 'MP2':
            df['MP2-IFIE'] = df['MP2-IFIE'].astype(float) * 627.5095
            df['PR-TYPE1'] = df['PR-TYPE1'].astype(float) * 627.5095
            df['GRIMME'] = df['GRIMME'].astype(float) * 627.5095
            df['JUNG'] = df['JUNG'].astype(float) * 627.5095
            df['HILL'] = df['HILL'].astype(float) * 627.5095

        # print('solv', solv)
        if len(solv) != 0:
            solvdf = pd.DataFrame(solv, columns=['I', 'J', 'Solv(ES)'])
            solvdf['I'] = solvdf['I'].astype(int)
            solvdf['J'] = solvdf['J'].astype(int)
            solvdf['Solv(ES)'] = solvdf['Solv(ES)'].astype(float)

            # print(df.head())
            # print(solvdf.head())
            df = pd.merge(df, solvdf, on=['I', 'J'], how='left')
        print(df.head())
        return df


    def getbssedf(self, ifie, solv=[]):
        '''get bsse data frame from ifie file

        get bsse data frame from ifie file

        Args:
            ifie (str): ifie file name
            solv (list): solvent name list

        Returns:
            df (pandas.DataFrame): bsse data frame
        '''

        df = pd.DataFrame(ifie, columns=self.bicolumn)
        df['I'] = df['I'].astype(int)
        df['J'] = df['J'].astype(int)
        df['DIST'] = df['DIST'].astype(float)
        df['HF-BSSE'] = df['HF-BSSE'].astype(float) * 627.5095

        if self.logMethod == 'MP2':
            df['MP2-BSSE'] = df['MP2-BSSE'].astype(float) * 627.5095
            df['PR-T1-BSSE'] = df['PR-T1-BSSE'].astype(float) * 627.5095
            df['GRIMME-BSSE'] = df['GRIMME-BSSE'].astype(float) * 627.5095
            df['JUNG-BSSE'] = df['JUNG-BSSE'].astype(float) * 627.5095
            df['HILL-BSSE'] = df['HILL-BSSE'].astype(float) * 627.5095

        df = df.drop(columns='DIST')
        # print(df.head())
        return df

    def getpiedadf(self, pieda):
        # print('l669', pieda[1])
        pidf = pd.DataFrame(pieda, columns=self.pcolumn)
        pidf['I'] = pidf['I'].astype(int)
        pidf['J'] = pidf['J'].astype(int)
        pidf['ES'] = pidf['ES'].astype(float)
        pidf['EX'] = pidf['EX'].astype(float)
        pidf['CT-mix'] = pidf['CT-mix'].astype(float)
        if self.abinit_ver == 'rev16' or self.abinit_ver == 'rev17':
            pidf['Solv(ES)'] = pidf['Solv(ES)'].astype(float)
        pidf['DI(MP2)'] = pidf['DI(MP2)'].astype(float)
        pidf['q(I=>J)'] = pidf['q(I=>J)'].astype(float)
        return pidf

    def getmomenedf(self, momene):
        df = pd.DataFrame(momene, columns=['Frag.', 'HF', 'MP2'])
        df['Frag.'] = df['Frag.'].astype(int)
        df['HF'] = df['HF'].astype(float)
        df['MP2'] = df['MP2'].astype(float)
        return df

    def getdimenedf(self, dimene):
        df = pd.DataFrame(dimene, columns=['I', 'J', 'DIMER-HF', 'DIMER-MP2'])
        df['I'] = df['I'].astype(int)
        df['J'] = df['J'].astype(int)
        df['DIMER-HF'] = df['DIMER-HF'].astype(float)
        df['DIMER-MP2'] = df['DIMER-MP2'].astype(float)
        return df

    def getpbpiedadf(self, pieda):
        pidf = pd.DataFrame(pieda, columns=self.pcolumnpb)
        pidf['I'] = pidf['I'].astype(int)
        pidf['J'] = pidf['J'].astype(int)
        pidf['ES'] = pidf['ES'].astype(float)
        pidf['EX'] = pidf['EX'].astype(float)
        pidf['CT-mix'] = pidf['CT-mix'].astype(float)
        pidf['Solv(ES)'] = pidf['Solv(ES)'].astype(float)
        pidf['DI(MP2)'] = pidf['DI(MP2)'].astype(float)
        pidf['q(I=>J)'] = pidf['q(I=>J)'].astype(float)
        pidf.drop(columns=['Solv(ES)'], inplace=True)

        return pidf


    def gettgtpidf_n2ffmatrix(self, mydf=None, is_pb=False):
        print('\n--- generate pieda', str(self.tgt1frag), str(self.tgt2frag), 'ffmatrix ---\n')
        esdf = pd.DataFrame(index=self.tgt2frag)
        exdf = pd.DataFrame(index=self.tgt2frag)
        ctdf = pd.DataFrame(index=self.tgt2frag)
        count = 0

        if mydf == None:
            df = self.pidf
        else:
            df = mydf
        for f1 in self.tgt1frag:
            wodimesapr_id = []
            tgtdf1 = df[(df['I'] == f1)].rename(columns={'I':'J', 'J':'I'})
            tgtdf2 = df[(df['J'] == f1)]
            # tgtdf = tgtdf1.append(tgtdf2)
            tgtdf = pd.concat([tgtdf1, tgtdf2])
            # print(tgtdf)

            tgtfrags = copy.deepcopy(self.tgt2frag)
            try:
                tgtfrags.remove(f1)
            except:
                pass
            tgtdf_filter = tgtdf[(tgtdf['I'].isin(tgtfrags)) | (tgtdf['J'].isin(tgtfrags))]

            fragis = tgtdf_filter['I'].values.tolist()
            fragjs = tgtdf_filter['J'].values.tolist()

            # print(tgtdf_filter)
            # pickup fragids from I and J(without dimer es approximation ID)
            for i in range(len(fragis)):
                if fragis[i] != f1:
                    wodimesapr_id.append(fragis[i])
                else:
                    wodimesapr_id.append(fragjs[i])
            # print('wodimesapr_id', wodimesapr_id)


            hfifie = 0
            mp2corr = 0
            prmp2corr = 0

            esbuf = tgtdf_filter['ES'].values.tolist()
            exbuf = tgtdf_filter['EX'].values.tolist()
            ctbuf = tgtdf_filter['CT-mix'].values.tolist()

            # complement values from ifdf
            es = []
            ex = []
            ct = []

            print('check frag', f1, 'pieda info')
            for i in range(len(self.tgt2frag)):
                tgtid = self.tgt2frag[i]
                # if tgtid  == f1:
                    # print('error!! target frag1 and target 2 is duplicate!!')
                    # sys.exit()
                if tgtid in wodimesapr_id:
                    es.append(esbuf[wodimesapr_id.index(tgtid)])
                    ex.append(exbuf[wodimesapr_id.index(tgtid)])
                    ct.append(ctbuf[wodimesapr_id.index(tgtid)])
                else:
                    if is_pb:
                        es.append(self.hfdf.loc[tgtid, str(f1)] - self.solvesdf.loc[tgtid, str(f1)])
                    else:
                        # print(tgtid, 'is no data')
                        es.append(self.hfdf.loc[tgtid, str(f1)])
                    ex.append(0.0)
                    ct.append(0.0)

            esdf[str(f1)] = es
            exdf[str(f1)] = ex
            ctdf[str(f1)] = ct

            count += 1

        print ('ES\n', esdf.head())
        print ('EX\n', exdf.head())
        print ('CT\n', ctdf.head())

        self.esdf = esdf
        self.exdf = exdf
        self.ctdf = ctdf

        return


    def gettgtdf_n2ffmatrix(self, mydf=None):
        # generate frags-frags matrix
        print('\n--- generate ifie', str(self.tgt1frag), str(self.tgt2frag), 'ffmatrix---\n')
        hfdf = pd.DataFrame(index=self.tgt2frag)
        distdf = pd.DataFrame(index=self.tgt2frag)
        solvesdf = pd.DataFrame(index=self.tgt2frag)

        if mydf == None:
            df = self.ifdf
        else:
            df = mydf
        count = 0

        if self.logMethod in ['MP2', 'HF+D']:
            mp2corrdf = pd.DataFrame(index=self.tgt2frag)
            prmp2corrdf = pd.DataFrame(index=self.tgt2frag)
            mp2tdf =  pd.DataFrame(index=self.tgt2frag)
            prmp2tdf = pd.DataFrame(index=self.tgt2frag)

            for f1 in self.tgt1frag:
                fragids = []
                tgtdf1 = df[(df['I'] == f1)].rename(columns={'I':'J', 'J':'I'})
                tgtdf2 = df[(df['J'] == f1)]
                # tgtdf = tgtdf1._append(tgtdf2)
                tgtdf = pd.concat([tgtdf1, tgtdf2])
                # print(tgtdf)

                tgtfrags = copy.deepcopy(self.tgt2frag)
                try:
                    tgtfrags.remove(f1)
                except:
                    pass
                tgtdf_filter = tgtdf[(tgtdf['I'].isin(tgtfrags)) | (tgtdf['J'].isin(tgtfrags))]

                # print(tgtdf_filter.columns.tolist())
                if f1 in self.tgt2frag:
                    # print ([i for i in range(len(tgtdf_filter.columns))])
                    adddf = pd.DataFrame(
                        [f1, f1] +
                        [0 for i in range(len(tgtdf_filter.columns)-2)],
                        index=tgtdf_filter.columns).T
                    # tgtdf_filter = tgtdf_filter.append(adddf).sort_values('I')
                    tgtdf_filter = pd.concat([tgtdf_filter, adddf]).sort_values('I')
                print(tgtdf_filter.head())

                hfifie = 0
                mp2corr = 0
                prmp2corr = 0
                hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
                mp2corr = tgtdf_filter['MP2-IFIE'].values.tolist()
                prmp2corr = tgtdf_filter['PR-TYPE1'].values.tolist()
                dist = tgtdf_filter['DIST'].values.tolist()

                mp2total = []
                prmp2total  = []
                for i in range(len(hfifie)):
                    mp2total.append(hfifie[i] + mp2corr[i])
                    prmp2total.append(hfifie[i] + prmp2corr[i])

                hfdf[str(f1)] = hfifie
                mp2corrdf[str(f1)] = mp2corr
                prmp2corrdf[str(f1)] = prmp2corr
                mp2tdf[str(f1)] = mp2total
                prmp2tdf[str(f1)] = prmp2total
                distdf[str(f1)] = dist

                if 'Solv(ES)' in tgtdf_filter.columns:
                    solves = tgtdf_filter['Solv(ES)'].values.tolist()
                    solvesdf[str(f1)] = solves

                count += 1

            print ('HF\n', hfdf.head())
            print ('MP2corr\n', mp2corrdf.head())
            print ('PRMP2corr\n', prmp2corrdf.head())
            print ('MP2total\n', mp2tdf.head())
            print ('PRMP2total\n', prmp2tdf.head())

            self.hfdf = hfdf
            self.mp2corrdf = mp2corrdf
            self.prmp2corrdf = prmp2corrdf
            self.mp2tdf = mp2tdf
            self.prmp2tdf = prmp2tdf
            self.distdf = distdf

            if 'Solv(ES)' in tgtdf_filter.columns:
                self.solvesdf = solvesdf

        elif self.logMethod == 'HF':
            for f1 in self.tgt1frag:
                fragids = []
                tgtdf = df[(df['I'] == f1) | (df['J'] == f1)]
                tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

                hfifie = 0
                hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
                dist = tgtdf_filter['DIST'].values.tolist()

                hfdf[str(f1)] = hfifie
                distdf[str(f1)] = dist

                if 'Solv(ES)' in tgtdf_filter.columns:
                    solves = tgtdf_filter['Solv(ES)'].values.tolist()
                    solvesdf[str(f1)] = solves

                count += 1

            print ('HF\n', hfdf.head())

            self.hfdf = hfdf
            self.distdf = distdf

            if 'Solv(ES)' in tgtdf_filter.columns:
                self.solvesdf = solvesdf

        elif self.logMethod == 'MP3':
            mp2corrdf = pd.DataFrame(index=self.tgt2frag)
            mp3corrdf = pd.DataFrame(index=self.tgt2frag)
            mp25corrdf = pd.DataFrame(index=self.tgt2frag)
            usermp3corrdf = pd.DataFrame(index=self.tgt2frag)
            mp2tdf =  pd.DataFrame(index=self.tgt2frag)
            mp3tdf =  pd.DataFrame(index=self.tgt2frag)
            mp25tdf =  pd.DataFrame(index=self.tgt2frag)
            usermp3tdf = pd.DataFrame(index=self.tgt2frag)

            for f1 in self.tgt1frag:
                fragids = []
                tgtdf = df[(df['I'] == f1) | (df['J'] == f1)]
                tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

                print(tgtdf_filter)
                hfifie = 0
                mp2corr = 0
                mp3corr = 0
                usermp3corr = 0
                hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
                mp2corr = tgtdf_filter['MP2-IFIE'].values.tolist()
                mp3corr = tgtdf_filter['MP3-IFIE'].values.tolist()
                usermp3corr = tgtdf_filter['USER-MP3'].values.tolist()
                dist = tgtdf_filter['DIST'].values.tolist()

                mp2total = []
                mp3total = []
                usermp3total  = []
                mp25corr = []
                mp25total = []
                for i in range(len(hfifie)):
                    mp2total.append(hfifie[i] + mp2corr[i])
                    mp3total.append(hfifie[i] + mp3corr[i])
                    usermp3total.append(hfifie[i] + usermp3corr[i])
                    mp25total.append(hfifie[i] + (mp2corr[i] + mp3corr[i])*0.5)
                    mp25corr.append((mp2corr[i] + mp3corr[i])*0.5)

                # print('hfifie', hfifie)
                # print('tgtfrag', self.tgt2frag)

                hfdf[str(f1)] = hfifie
                mp2corrdf[str(f1)] = mp2corr
                mp3corrdf[str(f1)] = mp3corr
                mp25corrdf[str(f1)] = mp25corr
                usermp3corrdf[str(f1)] = usermp3corr
                mp2tdf[str(f1)] = mp2total
                mp3tdf[str(f1)] = mp3total
                mp25tdf[str(f1)] = mp25total
                usermp3tdf[str(f1)] = usermp3total
                distdf[str(f1)] = dist


                count += 1

            print ('HF\n', hfdf.head())
            print ('MP3corr\n', mp3corrdf.head())
            print ('USER-MP3corr\n', usermp3corrdf.head())
            print ('MP3total\n', mp3tdf.head())
            print ('USER-MP3total\n', usermp3tdf.head())

            self.hfdf = hfdf
            self.mp2corrdf = mp2corrdf
            self.mp3corrdf = mp3corrdf
            self.mp25corrdf = mp25corrdf
            self.usermp3corrdf = usermp3corrdf
            self.mp2tdf = mp2tdf
            self.mp3tdf = mp3tdf
            self.mp25tdf = mp25tdf
            self.usermp3tdf = usermp3tdf
            self.distdf = distdf

        elif self.logMethod == 'CCPT':
            mp2corrdf = pd.DataFrame(index=self.tgt2frag)
            mp3corrdf = pd.DataFrame(index=self.tgt2frag)
            mp4corrdf = pd.DataFrame(index=self.tgt2frag)
            mp25corrdf = pd.DataFrame(index=self.tgt2frag)
            mp35corrdf = pd.DataFrame(index=self.tgt2frag)
            mp2tdf =  pd.DataFrame(index=self.tgt2frag)
            mp3tdf =  pd.DataFrame(index=self.tgt2frag)
            mp4tdf =  pd.DataFrame(index=self.tgt2frag)
            mp25tdf =  pd.DataFrame(index=self.tgt2frag)
            mp35tdf =  pd.DataFrame(index=self.tgt2frag)

            for f1 in self.tgt1frag:
                fragids = []
                tgtdf = df[(df['I'] == f1) | (df['J'] == f1)]
                tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

                hfifie = 0
                mp3corr = 0
                usermp3corr = 0
                hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
                mp2corr = tgtdf_filter['MP2-IFIE'].values.tolist()
                mp3corr = tgtdf_filter['MP3-IFIE'].values.tolist()
                mp4corr = tgtdf_filter['MP4-IFIE'].values.tolist()
                dist = tgtdf_filter['DIST'].values.tolist()


                mp2total = []
                mp3total = []
                mp4total  = []
                mp25corr = []
                mp35corr = []
                mp25total = []
                mp35total = []
                for i in range(len(hfifie)):
                    mp2total.append(hfifie[i] + mp2corr[i])
                    mp3total.append(hfifie[i] + mp3corr[i])
                    mp4total.append(hfifie[i] + mp4corr[i])
                    mp25corr.append((mp2corr[i] + mp3corr[i])*0.5)
                    mp35corr.append((mp2corr[i] + mp4corr[i])*0.5)
                    mp25total.append(hfifie[i] + (mp2corr[i] + mp3corr[i])*0.5)
                    mp35total.append(hfifie[i] + (mp2corr[i] + mp4corr[i])*0.5)

                # print('hfifie', hfifie)
                # print('tgtfrag', self.tgt2frag)

                hfdf[str(f1)] = hfifie
                mp2corrdf[str(f1)] = mp2corr
                mp3corrdf[str(f1)] = mp3corr
                mp4corrdf[str(f1)] = mp4corr
                mp25corrdf[str(f1)] = mp25corr
                mp35corrdf[str(f1)] = mp35corr
                mp2tdf[str(f1)] = mp2total
                mp3tdf[str(f1)] = mp3total
                mp4tdf[str(f1)] = mp4total
                mp25tdf[str(f1)] = mp25total
                mp35tdf[str(f1)] = mp35total
                distdf[str(f1)] = dist

                count += 1

#             print ('HF\n', hfdf.head())
#             print ('MP3corr\n', mp3corrdf.head())
#             print ('USER-MP3corr\n', usermp3corrdf.head())
#             print ('MP3total\n', mp3tdf.head())
#             print ('USER-MP3total\n', usermp3tdf.head())

            self.hfdf = hfdf
            self.mp2corrdf = mp2corrdf
            self.mp3corrdf = mp3corrdf
            self.mp4corrdf = mp4corrdf
            self.mp25corrdf = mp25corrdf
            self.mp35corrdf = mp35corrdf
            self.mp2tdf = mp2tdf
            self.mp3tdf = mp3tdf
            self.mp4tdf = mp4tdf
            self.mp25tdf = mp25tdf
            self.mp35tdf = mp35tdf
            self.distdf = distdf

        return

    def depth(self, k):
        if not k:
            return 0
        else:
            if isinstance(k, list):
                return 1 + max(self.depth(i) for i in k)
            else:
                return 0

    def gettgtdf_n2tfmatrix(self, i, df, f1):
        # gettgtdf_normal to times-frags
        if self.depth(self.tgt2frag) >= 2:
            tgt2frag = self.tgt2frag[i]
        else:
            tgt2frag = self.tgt2frag

        hfdf = pd.DataFrame(index=tgt2frag)
        mp2corrdf = pd.DataFrame(index=tgt2frag)
        prmp2corrdf = pd.DataFrame(index=tgt2frag)
        mp2tdf =  pd.DataFrame(index=tgt2frag)
        prmp2tdf = pd.DataFrame(index=tgt2frag)
        distdf = pd.DataFrame(index=tgt2frag)

        fragids = []
        tgtdf1 = df[(df['I'] == f1)].rename(columns={'I':'J', 'J':'I'})
        tgtdf2 = df[(df['J'] == f1)]
        tgtdf = pd.concat([tgtdf1, tgtdf2])

        tgtfrags = copy.deepcopy(tgt2frag)
        try:
            tgtfrags.remove(f1)
        except:
            pass
        tgtdf_filter = tgtdf[(tgtdf['I'].isin(tgtfrags)) | (tgtdf['J'].isin(tgtfrags))]


        if f1 in tgt2frag:
            # print ([i for i in range(len(tgtdf_filter.columns))])
            adddf = pd.DataFrame([f1, f1]+ [0 for j in range(len(tgtdf_filter.columns)-2)], index=tgtdf_filter.columns).T
            # tgtdf_filter = tgtdf_filter.append(adddf).sort_values('I')
            tgtdf_filter = pd.concat([tgtdf_filter, adddf]).sort_values('I')
        # print(tgtdf_filter.head())

        hfifie = 0
        mp2corr = 0
        prmp2corr = 0
        hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
        mp2corr = tgtdf_filter['MP2-IFIE'].values.tolist()
        prmp2corr = tgtdf_filter['PR-TYPE1'].values.tolist()
        dist = tgtdf_filter['DIST'].values.tolist()

        mp2total = []
        prmp2total  = []
        for j in range(len(hfifie)):
            mp2total.append(hfifie[j] + mp2corr[j])
            prmp2total.append(hfifie[j] + prmp2corr[j])

        # print('hfifie', len(hfifie))
        # print('hfifie', hfifie)
        hfdf[self.tgttimes[i]] = hfifie
        mp2corrdf[self.tgttimes[i]] = mp2corr
        prmp2corrdf[self.tgttimes[i]] = prmp2corr
        mp2tdf[self.tgttimes[i]] = mp2total
        prmp2tdf[self.tgttimes[i]] = prmp2total
        distdf[self.tgttimes[i]] = dist


        # print (hfdf.head())
        # print (mp2corrdf.head())
        # print (prmp2corrdf.head())
        # print (mp2tdf.head())
        # print (prmp2tdf.head())
        #print(distdf.head())

        return hfdf, mp2corrdf, prmp2corrdf, mp2tdf, prmp2tdf, distdf


    def gettgtpidf_n2tfmatrix(self, i, df, hfdf, f1):
        # define tgtfrag for molname - frags
        if self.depth(self.tgt2frag) >= 2:
            tgt2frag = self.tgt2frag[i]
        else:
            tgt2frag = self.tgt2frag

        esdf = pd.DataFrame(index=tgt2frag)
        exdf = pd.DataFrame(index=tgt2frag)
        ctdf = pd.DataFrame(index=tgt2frag)
        didf = pd.DataFrame(index=tgt2frag)

        tgtdf1 = df[(df['I'] == f1)].rename(columns={'I':'J', 'J':'I'})
        tgtdf2 = df[(df['J'] == f1)]
        # tgtdf = tgtdf1.append(tgtdf2)
        tgtdf = pd.concat([tgtdf1, tgtdf2])
            # print(tgtdf)

        tgtfrags = copy.deepcopy(tgt2frag)
        try:
            tgtfrags.remove(f1)
        except:
            pass
        tgtdf_filter = tgtdf[(tgtdf['I'].isin(tgtfrags)) | (tgtdf['J'].isin(tgtfrags))]

        fragis = tgtdf_filter['I'].values.tolist()
        fragjs = tgtdf_filter['J'].values.tolist()

        # pickup wodimesapr_id from I and J
        wodimesapr_id = []
        for j in range(len(fragis)):
            if fragis[j] != f1:
                wodimesapr_id.append(fragis[j])
            else:
                wodimesapr_id.append(fragjs[j])
        # print(wodimesapr_id)

        hfifie = 0
        mp2corr = 0
        prmp2corr = 0

        esbuf = tgtdf_filter['ES'].values.tolist()
        exbuf = tgtdf_filter['EX'].values.tolist()
        ctbuf = tgtdf_filter['CT-mix'].values.tolist()
        dibuf = tgtdf_filter['DI(MP2)'].values.tolist()


        # complement values from ifdf
        es = []
        ex = []
        ct = []
        di = []

        for j in range(len(tgt2frag)):
            tgtid = tgt2frag[j]
#             if tgtid == f1:
#                 continue
            if tgtid in wodimesapr_id:
                es.append(esbuf[wodimesapr_id.index(tgtid)])
                ex.append(exbuf[wodimesapr_id.index(tgtid)])
                ct.append(ctbuf[wodimesapr_id.index(tgtid)])
                di.append(dibuf[wodimesapr_id.index(tgtid)])
            else:
                es.append(hfdf.loc[tgtid, str(self.tgttimes[i])])
                ex.append(0.0)
                ct.append(0.0)
                di.append(0.0)

        esdf[self.tgttimes[i]] = es
        exdf[self.tgttimes[i]] = ex
        ctdf[self.tgttimes[i]] = ct
        didf[self.tgttimes[i]] = di

#         print (esdf.head())
#         print (exdf.head())
#         print (ctdf.head())

        # print('esdf\n', esdf)
        return esdf, exdf, didf, ctdf

        # return hfdf, mp2corrdf, prmp2corrdf, mp2tdf, prmp2tdf

    def gettgtdf_fd(self, df):

        # filter tgt1frag IFIE
        tgtdf = df[(df['I'] == self.tgt1frag[0]) | (df['J'] == self.tgt1frag[0])]
        if self.tgt2type == 'dist':
            print('--- ifie around tgt ', self.tgt1frag, self.dist, 'angstrom ---')
            tgtdf_filter = tgtdf[tgtdf['DIST'] < self.dist]
        elif self.tgt2type == 'dimer-es':
            print('--- ifie around tgt ', self.tgt1frag, 'without Dimer-es approximation ---')
            if self.f90soflag == True:
                tgtdf_filter = tgtdf[tgtdf['DIMER-ES'] == 0]
            else:
                tgtdf_filter = tgtdf[tgtdf['DIMER-ES'] == 'F']

        return tgtdf, tgtdf_filter

    def gettgtdf_ff(self, df, frag1, frag2):
        # print('--- ifie frag ', frag1, frag2, '----')
        frag1 = int(frag1)
        frag2 = int(frag2)
        tgtdf_filter = df[((df['I'] == frag1) & (df['J'] == frag2)) | ((df['I'] == frag2) & (df['J'] == frag1))]
        # tgtdf_filter = df[((df['I'] == frag1) & (df['J'] == frag2)) | ((df['I'] == frag2) & (df['J'] == frag1))]
        return tgtdf_filter


    def gettgtdf_ffs(self, df, frag1, frag2):
        # print('--- ifie frag ', frag1, frag2, '----')
        print(frag1, frag2)
        tgtdf_filter = df[((df['I'] == frag1) & (df['J'].isin(frag2))) | ((df['I'].isin(frag2)) & (df['J'] == frag1))]
        # tgtdf_filter = df[((df['I'] == frag1) & (df['J'] == frag2)) | ((df['I'] == frag2) & (df['J'] == frag1))]
        return tgtdf_filter


    def getifiesummol(self, df, molfrags, molid):
        tgtdf_filters = pd.DataFrame(columns=self.icolumn)
        print(self.dist)
        tgtdf = df[df['I'].isin(molfrags) | df['J'].isin(molfrags)]
        tgtdf_filters = tgtdf[tgtdf['DIST'] < self.dist]

        print(tgtdf_filters.head())

        # get I column value
        neighbor_i = [index for index, row in tgtdf_filters.groupby("I")]
        print(neighbor_i)
        # get J coulum value
        neighbor_j = [index for index, row in tgtdf_filters.groupby("J")]
        neighbors= list(set(neighbor_i + neighbor_j))
        neighbors.sort()
        alreadys = copy.deepcopy(molfrags)
        contactmolfrags = []
        for i in neighbors:
            if i in alreadys:
                continue
            molfrag_new = self.getmolfrags(i, df)
            contactmolfrags.append(molfrag_new)
            alreadys = alreadys + molfrag_new
            # print(alreadys)
        print('contactmolfrags\n', contactmolfrags)

        # print('-- ifie frag_mol --')
        '''
        ifie_frag_mol: df frag-frag in each mol
        ifiesums: mol-mol ifie for each mol
        '''

        ifdf_frag_mols = []
        for contactmolfrag in contactmolfrags:
            ifie_frag_mol = df[(df['I'].isin(contactmolfrag) & df['J'].isin(molfrags)) | (df['J'].isin(contactmolfrag) & df['I'].isin(molfrags))]
            # ifdf_frag_mols.append(ifie_frag_mol)
            ifdf_frag_mols = pd.concat([ifdf_frag_mols, ifie_frag_mol])

        #pieda
        for i in range(len(ifdf_frag_mols)):
            ifdf_frag_mols[i] = pd.merge(ifdf_frag_mols[i], self.pidf, on=['I', 'J'], how='left')
        print(ifdf_frag_mols[0])


        count = 0
        if self.abinit_ver == 'rev17' or self.abinit_ver == 'rev16':
            self.ifdfsumcolumn = [['HF-IFIE', 'MP2-IFIE', 'ES', 'EX', 'Solv(ES)', 'CT-mix', 'DI(MP2)', 'q(I=>J)']]

        HF_IFIE_sums = []
        MP2_IFIE_sums = []
        PR_TYPE1_sums = []
        GRIMME_sums = []
        JUNG_sums = []
        HILL_sums = []
        ES_sums = []
        EX_sums = []
        CT_sums = []
        DI_sums = []
        q_sums = []

        for datadf in ifdf_frag_mols:
            HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum = self.getsumdf(datadf)
            HF_IFIE_sums.append(HF_IFIE_sum)
            MP2_IFIE_sums.append(MP2_IFIE_sum)
            PR_TYPE1_sums.append(PR_TYPE1_sum)
            GRIMME_sums.append(GRIMME_sum)
            JUNG_sums.append(JUNG_sum)
            HILL_sums.append(HILL_sum)
            ES_sums.append(ES_sum)
            EX_sums.append(EX_sum)
            CT_sums.append(CT_sum)
            DI_sums.append(DI_sum)
            q_sums.append(q_sum)

        ifdf_mol_mol = pd.DataFrame(columns=self.ifdfsumcolumn).astype(float)
        # self.ifdfsumcolumn = ['HF-IFIE', 'MP2-IFIE', 'PR-TYPE1', 'GRIMME', 'JUNG', 'HILL', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']

        ifdf_mol_mol['I'] = contactmolfrags
        ifdf_mol_mol['J'] = [molfrags for i in range(len(HF_IFIE_sums))]
        ifdf_mol_mol['HF-IFIE'] = HF_IFIE_sums
        ifdf_mol_mol['MP2-IFIE'] = MP2_IFIE_sums
        ifdf_mol_mol['PR_TYPE1'] = PR_TYPE1_sums
        ifdf_mol_mol['GRIMME'] = GRIMME_sums
        ifdf_mol_mol['JUNG'] = JUNG_sums
        ifdf_mol_mol['HILL'] = HILL_sums
        ifdf_mol_mol['ES'] = ES_sums
        ifdf_mol_mol['EX'] = EX_sums
        ifdf_mol_mol['CT-mix'] = CT_sums
        ifdf_mol_mol['DI(MP2)'] = DI_sums
        ifdf_mol_mol['q(I=>J)'] = q_sums

        HF_IFIE_molsum = sum(HF_IFIE_sums)
        MP2_IFIE_molsum = sum(MP2_IFIE_sums)
        PR_TYPE1_molsum = sum(PR_TYPE1_sums)
        GRIMME_molsum = sum(GRIMME_sums)
        JUNG_molsum = sum(JUNG_sums)
        HILL_molsum = sum(HILL_sums)
        ES_molsum = sum(ES_sums)
        EX_molsum = sum(EX_sums)
        CT_molsum = sum(CT_sums)
        DI_molsum = sum(DI_sums)
        q_molsum = sum(q_sums)

        ifdf_molsum = pd.Series([HF_IFIE_molsum, MP2_IFIE_molsum, PR_TYPE1_molsum, GRIMME_molsum, JUNG_molsum, HILL_molsum, ES_molsum, EX_molsum, CT_molsum, DI_molsum, q_molsum], index=self.ifdfsumcolumn, name='mol'+str(molid))

        return contactmolfrags, ifdf_frag_mols, ifdf_mol_mol, ifdf_molsum

    def getsumdf(self, df):
        HF_IFIE_sum = df['HF-IFIE'].sum()
        MP2_IFIE_sum = df['MP2-IFIE'].sum()
        PR_TYPE1_sum = df['PR-TYPE1'].sum()
        GRIMME_sum = df['GRIMME'].sum()
        JUNG_sum = df['JUNG'].sum()
        HILL_sum = df['HILL'].sum()

        ES_sum = df['ES'].sum()
        EX_sum = df['EX'].sum()
        CT_sum = df['CT-mix'].sum()
        DI_sum = df['DI(MP2)'].sum()
        q_sum = df['q(I=>J)'].sum()

        return HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum


    def getpitgtdf(self, pidf, ifdf_filter):
        print('--- pieda near tgt ', self.dist, 'angstrom ----')
        # print('--- pieda for tgt frag ----')
        tgt1frag = self.tgt1frag[0]
        tgtdf_filter = pidf[(pidf['I'] == tgt1frag) |(pidf['J'] == tgt1frag)]

#         fragids = []
#         tgtdf = df[(df['I'] == self.tgt1frag) | (df['J'] == self.tgt1frag)]
#         tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

        fragids = []
        fragis = tgtdf_filter['I'].values.tolist()
        fragjs = tgtdf_filter['J'].values.tolist()
        for i in range(len(fragis)):
            if fragis[i] != tgt1frag:
                fragids.append(fragis[i])
            else:
                fragids.append(fragjs[i])
        print('inside Dimter-ES region', fragids)

        allids = []
        allis = ifdf_filter['I'].values.tolist()
        alljs = ifdf_filter['J'].values.tolist()

        for i in range(len(allis)):
            if allis[i] != tgt1frag:
                allids.append(allis[i])
            else:
                allids.append(alljs[i])
        # print('all pair', allids)

        # complement values from ifdf

        df = ifdf_filter
        pitgtdf_filters = pd.DataFrame(columns=self.pcolumn)
        count = 0
        for i in range(len(allids)):
            pitgtdf_filter = self.gettgtdf_ff(tgtdf_filter, tgt1frag, allids[i])
            if len(pitgtdf_filter) == 0:
                esdata = df[(df['I'] == allids[i]) |(df['J'] == allids[i])]['HF-IFIE'].values.tolist()[0]
                pitgtdf_filters.loc[str(count)] = [int(tgt1frag), int(allids[i]), esdata, 0.0, 0.0, 0.0, 0.0]
            else:
                pitgtdf_filters = pitgtdf_filters.append(pitgtdf_filter)
            count +=1

        # print(pitgtdf_filters.head())

        return pitgtdf_filters


    def read_ifiepieda(self, fname):
        ''' read ifie and pieda data from file

        read ifie and pieda data from file

        Args:
            fname (str): file name

        Returns:
            ifie (list): ifie data
            pieda (list): pieda data
        '''
        # print('start getifiepieda')
        ifie = []
        count = 0
        pieda = []
        pcount = 0
        momcount = 0
        dimcount = 0
        momene = []
        dimene = []
        pflag = False
        bsseflag = False
        bssecount = 0
        bsse = []

        if not os.path.exists(fname):
            print("can't open", fname)
            return ifie, pieda, momene, dimene

        file = open(fname, "rt")

        flag = False
        momflag = False
        dimflag = False
        for line in file:
            # itemList = text[i][:-1].split()
            # itemList = file.readline().strip().split()
            itemList = line.strip().split()

            # print itemList
            if len(itemList) < 2:
                continue
            if itemList[1:3] == ['MONOMER', 'ENERGY']:
                momflag = True
                continue
            if itemList[1:3] == ['DIMER', 'ENERGY']:
                dimflag = True
                continue
            if itemList[1:3] == ['DIMER', '<S^2>']:
                dimflag = False
            if dimflag is True:
                dimcount += 1
            if dimflag is True and dimcount > 2:
                # print('DIMER Energy', itemList)
                dimene.append(itemList)

            if momflag:
                # print(itemList)
                momcount += 1
                if momcount == 1 + int(self.tgt1frag[0]):
                    momene.append(itemList)
                    momflag = False

            if itemList[1] == 'MP2-IFIE' or itemList[1] == 'HF-IFIE':
                flag = True
                # head.append(itemList)
                continue
            if itemList[1] == 'PIEDA':
                flag = False
                pflag = True
                # print('pieda start!!')
                continue
            if flag is True:
                count += 1
            if flag is True and count > 2:
                if self.logMethod == 'HF':
                    ifie.append(itemList[:-2])
                else:
                    ifie.append(itemList)
            if len(itemList) < 2:
                continue

            # after pieda or BSSE (break)
            if itemList[1] == 'Mulliken':
                # flag = False
                break
            # for BSSE
            if pflag is True and itemList[:5] == ['##','BSSE', 'for','non-bonding','MP2-IFIE']:
                pflag = False
                # print('pieda end! next is BSSE')
                continue
            if itemList[:4] == ['##','BSSE', 'for', 'MP2-IFIE']:
                bsseflag = True
                # print('BSSE start!')
                continue
            if bsseflag is True:
                bssecount += 1
            if bsseflag is True and bssecount > 2:
                bsse.append(itemList)

            # for pieda
            if pflag is True:
                pcount += 1
            if pflag is True and pcount > 2:
                pieda.append(itemList)

        if not flag and not pflag and not bsseflag:
            print("can't read ifie", fname.split("/")[-1])
            return [], [], [], []

        for i in range(len(ifie)):
            if float(ifie[i][4]) < -2:
                ifie[i][4] = 0.0
                if self.logMethod != 'HF':
                    ifie[i][5] = 0.0
                    ifie[i][6] = 0.0

        file.close()
        return ifie, pieda, momene, dimene, bsse


    def read_pbifiepieda(self, fname):
        ifie = []
        count = 0
        pieda = []
        pcount = 0
        pflag = False
        sflag = False

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
            if itemList[1] == 'MP2-IFIE' or itemList[1] == 'HF-IFIE':
                flag = True
                ifie = []
                pieda = []
                # head.append(itemList)
                continue
            if itemList[1] == 'PIEDA':
                flag = False
                pflag = True
                continue
            if flag is True:
                count += 1
            if flag is True and count > 2:
                if self.logMethod == 'HF':
                    ifie.append(itemList[:-2])
                else:
                    ifie.append(itemList)
            # print itemList
            if len(itemList) < 2:
                continue

            if itemList[1] == 'PIEDA':
                pflag = True
                # head.append(itemList)
                continue
            if itemList[1] == 'Mulliken':
                flag = False
                pflag = False
                count = 0
                pcount = 0
                continue
            if pflag is True:
                pcount += 1
            if pflag is True and pcount > 2:
                pieda.append(itemList)

            if itemList[1] == 'SOLVENT-SCREENING':
                sflag = True
                solv = []
                scount = 0
            if sflag == True:
                scount += 1

            if itemList[1:3] == ['NONPOLAR', 'CONTRIBUTION']:
                sflag = False
                break

            if sflag is True and scount > 5:
                solv.append([itemList[1], itemList[2], itemList[12]])

        if not flag and not pflag:
            try:
                print("can't read ifie", fname.split("/")[1])
            except:
                pass

        for i in range(len(ifie)):
            # print(ifie[i][4])
            if float(ifie[i][4]) < -2:
                ifie[i][4] = 0.0
                if self.logMethod != 'HF':
                    ifie[i][5] = 0.0
                    ifie[i][6] = 0.0

        return ifie, pieda, solv


    def read_ifiepiedas(self, tgtlog):
        '''read ifie and pieda from log file
        Args:
            tgtlog (str): log file name
        Returns:
            ifdfs (list): ifie dataframe
            pidfs (list): pieda dataframe
            momenedf (list): momene dataframe
            dimenedf (list): dimene dataframe
            bssedfs (list): bsse dataframe
        '''
        # tgtlog, tgttime = args
        print('read', tgtlog)

        # getifie
        # HF or MP2 is supported on py function
        ifie, pieda, momene, dimene, bsse = self.read_ifiepieda(tgtlog)
        # print('l1613', 'ifie', ifie[0], 'pieda', pieda[0], 'mom', momene)

        # get dataframe list
        ifdfs = self.getifiedf(ifie)
        pidfs = self.getpiedadf(pieda)

        self.is_bsse = True
        if self.is_bsse:
            bssedfs = self.getbssedf(bsse)
        else:
            bssedfs = []
        if self.is_momdimene:
            momenedf = self.getmomenedf(momene)
            dimenedf = self.getdimenedf(dimene)
        else:
            momenedf = []
            dimenedf = []
        # print('ifdfs', ifdfs)
        # print('pidfs', pidfs)

        return [ifdfs, pidfs, momenedf, dimenedf, bssedfs]


    def getfiltifpiff(self, i, ifdf, pidf):
        '''get filtered ifie and pieda dataframe

        Args:
            i (int): index of rec
            ifdf (dataframe): ifie dataframe
            pidf (dataframe): pieda dataframe
        Returns:
            ifdf_filter (dataframe): filtered ifie dataframe
            pidf_filter (dataframe): filtered pieda dataframe
        '''

        # in (class var: molnames(rec i), tgt2frag, tgt1frag, tgttimes, matrixtype
        #     local var:i,  ifdf, pidf
        # out tgtdf_filter, pitgtdf_filter

        # ifie
        ifdf_filter = pd.DataFrame()
        ifdf_filter = self.gettgtdf_ff(ifdf, self.tgt1frag[0], self.tgt2frag[0])
        # print(ifdf_filter)
        ifdf_filter = pd.merge(ifdf_filter, pidf, on=['I', 'J'], how='left')
        ifdf_filter = ifdf_filter.rename(index={ifdf_filter.index[0]:self.tgttimes[i]})

        # pieda
#         pitgtdf_filter = self.gettgtdf_ff(pidf, self.tgt1frag[0], self.tgt2frag[0])
#         if len(pitgtdf_filter) == 0:
#             esdata = ifdf_filter['HF-IFIE'][0]
#             tmpps = pd.Series([self.tgt1frag[0], self.tgt2frag[0], esdata, 0.0, 0.0, 0.0, 0.0], index=self.pcolumn, name =self.tgttimes[i])
#             pitgtdf_filter = pitgtdf_filter.append(tmpps)
#         else:
#             pitgtdf_filter = pitgtdf_filter.rename(index={pitgtdf_filter.index[0]:self.tgttimes[i]})
#
#             # print(pitgtdf_filter)

        return ifdf_filter


    def getfiltifpifm(self, i, ifdf, pidf):
        '''get filtered ifie and pieda dataframe

        Args:
            i (int): index of rec
            ifdf (dataframe): ifie dataframe
            pidf (dataframe): pieda dataframe

        Returns:
            ifdf_filter (dataframe): filtered ifie dataframe
            pidf_filter (dataframe): filtered pieda dataframe
        '''

        # in (class var: log(rec i), nf(rec i), molnames(rec i), tgt2molname, tgt1frag,
        #     local var: ifdf, pidf
        # out tgtdf_filter, tgtifdfsum
        #     pitgtdf_filter, pitgtdfsum

        print('### read frag info ###')

        # molfragss
        molfragss = self.getallmolfrags(self.tgtlogs[i], ifdf, self.nfs[i])
        # molfragss fragment ids per mol (fragment connect is judged by checking frag-frag distance)
        print('molfragss', molfragss)
        print('len molfragss', len(molfragss))

        # get tgt frag id
        # IFIE
        molnames_inrec = self.molnames_perrec[i]
        rname_perfraginrec = self.resnamenonums_perfrag[i]

        tgtmolfrags = []
        # print('molnames', molnames_inrec)
        # print('len molnames', len(molnames_inrec))

        print('molnames_perfrag', rname_perfraginrec)
        print('len molnames_perfrag', len(rname_perfraginrec))

        rname_permolinrec = copy.deepcopy(molfragss)
        for j in range(len(rname_permolinrec)):
            # if type(molfragss[j]) == list:
            for k in range(len(molfragss[j])):
                num = copy.deepcopy(molfragss[j][k]) - 1
                rname_permolinrec[j][k] = rname_perfraginrec[num]

        for j in range(len(rname_permolinrec)):
            try:
                if self.tgt2molname in rname_permolinrec[j]:
                    tgtmolfrags += molfragss[j]
            except:
                continue
        print(tgtmolfrags)

        frag1 = self.tgt1frag
        if type(frag1) == int:
            frag1s = [frag1]
            print('tgtfrag1s', frag1s)
        else:
            frag1s = copy.deepcopy(frag1)

        self.frag1s = frag1s
        ifdf_filter = pd.DataFrame()
        ifdf_filters = []
        ifdfsums = []
        for frag1p in frag1s:
            ifdf_filter = self.gettgtdf_ffs(ifdf, frag1p, tgtmolfrags)
            print('ifdf_filter\n', ifdf_filter.head())

            # merge ifie and pieda
            ifdf_filter = pd.merge(ifdf_filter, pidf, on=['I', 'J'], how='left')

            ## screening dist
            if self.dist != 1000.0:
                ifdf_filter = ifdf_filter[ifdf_filter['DIST'] < self.dist]

            HF_IFIE_sum = ifdf_filter['HF-IFIE'].sum()
            MP2_IFIE_sum = ifdf_filter['MP2-IFIE'].sum()
            PR_TYPE1_sum = ifdf_filter['PR-TYPE1'].sum()
            GRIMME_sum = ifdf_filter['GRIMME'].sum()
            JUNG_sum = ifdf_filter['JUNG'].sum()
            HILL_sum = ifdf_filter['HILL'].sum()

            ES_sum = ifdf_filter['ES'].sum()
            EX_sum = ifdf_filter['EX'].sum()
            CT_sum = ifdf_filter['CT-mix'].sum()
            DI_sum = ifdf_filter['DI(MP2)'].sum()
            q_sum = ifdf_filter['q(I=>J)'].sum()

            ifdf_filter['TIMES'] = self.tgttimes[i]

            ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum], index=self.ifdfsumcolumn, name=self.tgttimes[i])
            print('ifdfsum\n', ifdfsum)
            ifdf_filters.append(ifdf_filter)
            ifdfsums.append(ifdfsum)

        return ifdf_filters, ifdfsums


    def getfiltifpifd(self, i, ifdf, pidf, momenedf=None, dimenedf=None, bssedf=None):
        # in class var: tgttime
        #    local var: i, ifdf,pidf
        # out local var: tgtifdfsum, tgtdf_filter
        #                pitgtdfsum, pitgtedf

        # get tgt frag id
        ifdf, ifdf_filter = self.gettgtdf_fd(ifdf)

        # merge ifiedf and pieda df
        ifdf_filter = pd.merge(ifdf_filter, pidf, on=['I', 'J'], how='left')

        # merge filtered-ifie data and dimer energy data
        if self.is_momdimene:
            ifdf_filter = pd.merge(ifdf_filter, dimenedf, on=['I', 'J'], how='left')

        # merge ifiedf and pieda df
        if self.is_bsse:
            ifdf_filter = pd.merge(ifdf_filter, bssedf, on=['I', 'J'], how='left')

        print(ifdf_filter.head())
        ifdf_filter['TIMES'] = self.tgttimes[i]
        HF_IFIE_sum = ifdf_filter['HF-IFIE'].sum()
        MP2_IFIE_sum = ifdf_filter['MP2-IFIE'].sum()
        PR_TYPE1_sum = ifdf_filter['PR-TYPE1'].sum()
        GRIMME_sum = ifdf_filter['GRIMME'].sum()
        JUNG_sum = ifdf_filter['JUNG'].sum()
        HILL_sum = ifdf_filter['HILL'].sum()
        ES_sum = ifdf_filter['ES'].sum()
        EX_sum = ifdf_filter['EX'].sum()
        CT_sum = ifdf_filter['CT-mix'].sum()
        DI_sum = ifdf_filter['DI(MP2)'].sum()
        q_sum = ifdf_filter['q(I=>J)'].sum()

        try:
            HF_BSSE_sum = ifdf_filter['HF-BSSE'].sum()
            MP2_BSSE_sum = ifdf_filter['MP2-BSSE'].sum()
        except:
            HF_BSSE_sum = 0
            MP2_BSSE_sum = 0

        # data for bsse
        try:
            momenetgtdf = momenedf[momenedf['Frag.'] == self.momfrag]
            momene_tgt = momenetgtdf['HF'][0] + momenetgtdf['MP2'][0]
            momlabel = 'MonomerEnergy(' + str(self.momfrag) + ')'
        except:
            momene_tgt = None
            momlabel = 'MonomerEnergy(' + str(self.momfrag) + ')'

        try:
            # dimene_tgt = momenedf['HF'][0] + momenedf['MP2'][0]
            dimenetgtdf = dimenedf[((dimenedf['I'] == self.dimfrag1) & (dimenedf['J'] == self.dimfrag2)) | ((dimenedf['I'] == self.dimfrag2) & (dimenedf['J'] == self.dimfrag1))]
            dimene_tgt = dimenetgtdf['DIMER-HF'][0] + dimenetgtdf['DIMER-MP2'][0]
            dimlabel = 'DimerEnergy(' + str(self.dimfrag1) + '-' + str(self.dimfrag2) + ')'

        except:
            dimene_tgt = None
            dimlabel = 'DimerEnergy(' + str(self.dimfrag1) + '-' + str(self.dimfrag2) + ')'


        ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum, momene_tgt, dimene_tgt, HF_BSSE_sum, MP2_BSSE_sum], index=self.ifdfsumcolumn + [momlabel, dimlabel] + ['HF-BSSE', 'MP2-BSSE'], name=self.tgttimes[i])
        # print(ifdfsum)

        # pieda
#         frags = self.frags
#         if self.fragmode != 'manual':
#             # print('len_frags', len(frags))
#             #assign resname(e.g. Gly6)
#             for j in range(1, len(frags) + 1):
#                 ifdf_filter.I = ifdf_filter.I.replace(j, frags[j-1])
#                 ifdf_filter.J = ifdf_filter.J.replace(j, frags[j-1])

        # print(ifdf_filter)
        return ifdf_filter, ifdfsum

    def read_ifpif90(self, tgtlog):
        '''read ifpif90.so and call readifiepieda_
        Args:
            tgtlog (str): target log file name
        Returns:
            ifdf_filters (list): list of ifdf_filter
        '''

        print('read', tgtlog)
        if not os.path.exists(tgtlog):
            print('Warning:', tgtlog, 'is not exist: skip data')
            return []
        f = np.ctypeslib.load_library(self.f90sofile, ".")

        f.readifiepieda_.argtypes = [
            c_char_p,
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.float64),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            ]

        instr = tgtlog
        enc_str = instr.encode('utf-8')
        instr = create_string_buffer(enc_str)

        ifi = np.empty(100000000, dtype=np.int32)
        ifj = np.empty(100000000, dtype=np.int32)
        pii = np.empty(100000000, dtype=np.int32)
        pij = np.empty(100000000, dtype=np.int32)
        dist = np.empty(100000000, dtype=np.float64)
        hfifie = np.empty(100000000, dtype=np.float64)
        mp2ifie = np.empty(100000000, dtype=np.float64)
        prtype1 = np.empty(100000000, dtype=np.float64)
        grimme = np.empty(100000000, dtype=np.float64)
        jung = np.empty(100000000, dtype=np.float64)
        hill = np.empty(100000000, dtype=np.float64)
        es = np.empty(100000000, dtype=np.float64)
        ex = np.empty(100000000, dtype=np.float64)
        ct = np.empty(100000000, dtype=np.float64)
        di = np.empty(100000000, dtype=np.float64)
        erest = np.empty(100000000, dtype=np.float64)
        qval = np.empty(100000000, dtype=np.float64)
        fdimesint = np.empty(100000000, dtype=np.int32)
        ifpair = np.empty(1, dtype=np.int32)
        pipair = np.empty(1, dtype=np.int32)

        f.readifiepieda_(instr, ifi, ifj, pii, pij, dist, hfifie, mp2ifie, prtype1, grimme, jung, hill, es, ex, ct, di, erest, qval, fdimesint, ifpair, pipair)

        if ifpair == 0:
            print('Warning:', tgtlog, 'is not converged: skip data')
            return []
        #check
        # print ('i', ifi[0], ifj[0])
        # print ('j', ifi[1], ifj[1])
        # print ('dist', dist[0:2])
        # print ('dist', hfifie[0:2])
        # print ('fdimesint', fdimesint[0:2])
        print ('ifpair', ifpair, 'pipair', pipair)

        # fdimesstr = []
        # for i in range(ifpair[0]):
        #     if fdimesint[i] == 1:
        #         fdimesstr.append('F')
        #     elif fdimesint[i] == 2:
        #         fdimesstr.append('T')

        # MP3 case
        if self.logMethod == 'MP3':
            self.icolumn = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE', 'USER-MP2', 'MP3-IFIE','USER-MP3', 'PADE[2/1]' ]

            ifdf = pd.DataFrame(columns=self.icolumn)
            ifdf['I'] = copy.deepcopy(ifi[:ifpair[0]])
            ifdf['J'] = copy.deepcopy(ifj[:ifpair[0]])
            ifdf['DIST'] = copy.deepcopy(dist[:ifpair[0]])
            ifdf['DIMER-ES'] = copy.deepcopy(fdimesint[:ifpair[0]])
            ifdf['HF-IFIE'] = copy.deepcopy(hfifie[:ifpair[0]])
            ifdf['MP2-IFIE'] = copy.deepcopy(mp2ifie[:ifpair[0]])
            ifdf['USER-MP2'] = copy.deepcopy(prtype1[:ifpair[0]])
            ifdf['MP3-IFIE'] = copy.deepcopy(grimme[:ifpair[0]])
            ifdf['USER-MP3'] = copy.deepcopy(jung[:ifpair[0]])
            ifdf['PADE[2/1]'] = copy.deepcopy(hill[:ifpair[0]])

        # CCPT case
        elif self.logMethod == 'CCPT':
            self.icolumn = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE', 'GRIMME-MP2', 'MP3-IFIE','GRIMME-MP3', 'MP4-IFIE' ]

            ifdf = pd.DataFrame(columns=self.icolumn)
            ifdf['I'] = copy.deepcopy(ifi[:ifpair[0]])
            ifdf['J'] = copy.deepcopy(ifj[:ifpair[0]])
            ifdf['DIST'] = copy.deepcopy(dist[:ifpair[0]])
            ifdf['DIMER-ES'] = copy.deepcopy(fdimesint[:ifpair[0]])
            ifdf['HF-IFIE'] = copy.deepcopy(hfifie[:ifpair[0]])
            ifdf['MP2-IFIE'] = copy.deepcopy(mp2ifie[:ifpair[0]])
            ifdf['GRIMME-MP2'] = copy.deepcopy(prtype1[:ifpair[0]])
            ifdf['MP3-IFIE'] = copy.deepcopy(grimme[:ifpair[0]])
            ifdf['GRIMME-MP3'] = copy.deepcopy(jung[:ifpair[0]])
            ifdf['MP4-IFIE'] = copy.deepcopy(hill[:ifpair[0]])

        # MP2 case
        else:
            ifdf = pd.DataFrame(columns=self.icolumn)
            ifdf['I'] = copy.deepcopy(ifi[:ifpair[0]])
            ifdf['J'] = copy.deepcopy(ifj[:ifpair[0]])
            ifdf['DIST'] = copy.deepcopy(dist[:ifpair[0]])
            ifdf['DIMER-ES'] = copy.deepcopy(fdimesint[:ifpair[0]])
            ifdf['HF-IFIE'] = copy.deepcopy(hfifie[:ifpair[0]])
            ifdf['MP2-IFIE'] = copy.deepcopy(mp2ifie[:ifpair[0]])
            ifdf['PR-TYPE1'] = copy.deepcopy(prtype1[:ifpair[0]])
            ifdf['GRIMME'] = copy.deepcopy(grimme[:ifpair[0]])
            ifdf['JUNG'] = copy.deepcopy(jung[:ifpair[0]])
            ifdf['HILL'] = copy.deepcopy(hill[:ifpair[0]])

        # PIEDA
        pidf = pd.DataFrame(columns=self.pcolumn)
        pidf['I'] = copy.deepcopy(pii[:pipair[0]])
        pidf['J'] = copy.deepcopy(pij[:pipair[0]])
        pidf['ES'] = copy.deepcopy(es[:pipair[0]])
        pidf['EX'] = copy.deepcopy(ex[:pipair[0]])
        pidf['CT-mix'] = copy.deepcopy(ct[:pipair[0]])
        pidf['DI(MP2)'] = copy.deepcopy(di[:pipair[0]])
        # LRD
        if self.is_disp:
            pidf['Erest'] = copy.deepcopy(erest[:pipair[0]])
        pidf['q(I=>J)'] = copy.deepcopy(qval[:pipair[0]])

        return [ifdf, pidf]


    def read_ifpimulti(self, i):
        ''' read ifpi and pieda from multi log file
        Args:
            i: log number
        Returns:
            ifdf: ifpi dataframe
            pidf: pieda dataframe
        '''

        # i: log number
        # read ifpi f90
        momenedf = None
        if self.f90soflag:
            # get ifie using f90 module
            print("use fortran library")
            dfs = self.read_ifpif90(self.tgtlogs[i])
            if len(dfs) == 0:
                return None
            ifdf, pidf = dfs

        # note: nof90(py) mode only can get momnomer energy
        else:
            dfs = self.read_ifiepiedas(self.tgtlogs[i])
            if len(dfs) == 0:
                return None
            ifdf, pidf, momenedf, dimenedf, bssedf = dfs

        # IFIE pieda filter
        tgt2type = self.tgt2type
        print('tgttimes', self.tgttimes[i])
        print('tgt2type', tgt2type)

        # frag mode
        if tgt2type == 'frag':
            if self.matrixtype == 'times-frags':
                hfdfs = []
                mp2corrdfs =[]
                prmp2corrdfs =[]
                mp2tdfs =[]
                prmp2tdfs =[]
                esdfs = []
                exdfs = []
                ctdfs = []
                didfs = []
                distdfs = []
                erestdfs = []
                for tmptgt1 in self.tgt1frag:
                    hfdf, mp2corrdf, prmp2corrdf, mp2tdf, prmp2tdf, distdf = self.gettgtdf_n2tfmatrix(i, ifdf, tmptgt1)
                    esdf, exdf, didf, ctdf = self.gettgtpidf_n2tfmatrix(i, pidf, hfdf, tmptgt1)
                    hfdfs.append(hfdf)
                    mp2corrdfs.append(mp2corrdf)
                    prmp2corrdfs.append(prmp2corrdf)
                    mp2tdfs.append(mp2tdf)
                    prmp2tdfs.append(prmp2tdf)
                    distdfs.append(distdf)
                    esdfs.append(esdf)
                    exdfs.append(exdf)
                    ctdfs.append(ctdf)
                    didfs.append(didf)
                    erestdfs.append(prmp2corrdf - didf)
                return hfdfs, mp2corrdfs, prmp2corrdfs, mp2tdfs, prmp2tdfs, esdfs, exdfs, ctdfs, distdfs, didfs, erestdfs

            else:
                ifdf_filter = self.getfiltifpiff(i, ifdf, pidf)
                iftgtdfsum = None
                return ifdf_filter

        # molname mode
        if tgt2type == 'molname':
            ifdf_filters, ifdfsums = self.getfiltifpifm(i, ifdf, pidf)
            return ifdf_filters, ifdfsums

        # dist or dimer-es mode
        if tgt2type in ['dist', 'dimer-es']:
            if self.is_momdimene:
                if self.is_bsse:
                    ifdf_filter, ifdfsum  = self.getfiltifpifd(i, ifdf, pidf, momenedf, dimenedf, bssedf)
                    return ifdf_filter, ifdfsum
                else:
                    ifdf_filter, ifdfsum  = self.getfiltifpifd(i, ifdf, pidf, momenedf, dimenedf)
                    return ifdf_filter, ifdfsum
            else:
                if self.is_bsse:
                    ifdf_filter, ifdfsum  = self.getfiltifpifd(i, ifdf, pidf, None, None, bssedf)
                    return ifdf_filter, ifdfsum
                else:
                    ifdf_filter, ifdfsum  = self.getfiltifpifd(i, ifdf, pidf)
                    return ifdf_filter, ifdfsum

        else:
            print('tgt2type error!!')
            return

    def readmultiifie(self):
        # check disp(LRD)
        self.is_disp = self.getisdisp(self.tgtlogs[0])
        print('## read multi mode')
        st = time.time()
        p = Pool(self.pynp)

        # main-read module
        ifpidfs = p.map(self.read_ifpimulti, [i for i in range(len(self.tgtlogs))])
        # i: time j:type, k:tgt1frag
        # [[hfdf[time 1][frag 1, frag2...], , mp2df[time 1], ...], [hfdf[time 2], mp2df[time 2], ...], ...]
        pd.set_option('display.width', 500)
        print('valid data num:', len(ifpidfs))
        # print(ifpidfs)

        # delete label of unfinished log
        delids = []
        for i in range(len(ifpidfs)):
            try:
                if ifpidfs[i] is None:
                    delids.append(i)
            except:
                pass
        dellist = lambda items, indexes: [item for index, item in enumerate(items) if index not in indexes]
        ifpidfs = dellist(ifpidfs, delids)

        # filter section
        if self.tgt2type == 'frag':
            if self.matrixtype == 'times-frags':
                self.hfdfs = []
                self.mp2corrdfs = []
                self.prmp2corrdfs = []
                self.mp2tdfs = []
                self.prmp2tdfs = []
                self.esdfs = []
                self.exdfs = []
                self.ctdfs = []
                self.distdfs = []
                self.didfs = []
                self.erestdfs = []
                print('readable file num:', len(ifpidfs))
                for i in range(len(self.tgt1frag)):
                    hfdf = pd.concat([dfs[0][i] for dfs in ifpidfs], axis=1)
                    mp2corrdf = pd.concat([dfs[1][i] for dfs in ifpidfs], axis=1)
                    prmp2corrdf = pd.concat([dfs[2][i] for dfs in ifpidfs], axis=1)
                    mp2tdf = pd.concat([dfs[3][i] for dfs in ifpidfs], axis=1)
                    prmp2tdf = pd.concat([dfs[4][i] for dfs in ifpidfs], axis=1)
                    esdf = pd.concat([dfs[5][i] for dfs in ifpidfs], axis=1)
                    exdf = pd.concat([dfs[6][i] for dfs in ifpidfs], axis=1)
                    ctdf = pd.concat([dfs[7][i] for dfs in ifpidfs], axis=1)
                    distdf = pd.concat([dfs[8][i] for dfs in ifpidfs], axis=1)
                    didf = pd.concat([dfs[9][i] for dfs in ifpidfs], axis=1)
                    erestdf = pd.concat([dfs[10][i] for dfs in ifpidfs], axis=1)

                    self.hfdfs.append(hfdf)
                    self.mp2corrdfs.append(mp2corrdf)
                    self.prmp2corrdfs.append(prmp2corrdf)
                    self.mp2tdfs.append(mp2tdf)
                    self.prmp2tdfs.append(prmp2tdf)
                    self.esdfs.append(esdf)
                    self.exdfs.append(exdf)
                    self.ctdfs.append(ctdf)
                    self.distdfs.append(distdf)
                    self.didfs.append(didf)
                    self.erestdfs.append(erestdf)
            else:
                self.ifdf_filters = pd.DataFrame()
                for dfs in ifpidfs:
                    # self.ifdf_filters = self.ifdf_filters.append(dfs)
                    self.ifdf_filters = pd.concat([self.ifdf_filters, dfs])

        if self.tgt2type in ['dist', 'dimer-es']:
            self.ifdf_filters = pd.DataFrame()
            self.ifdfsum = pd.DataFrame(columns=self.ifdfsumcolumn).astype(float)

            for dfs in ifpidfs:
                self.ifdf_filters = pd.concat([self.ifdf_filters,
                                               dfs[0]])
                sumdf = dfs[1].to_frame().T
                sumdf = sumdf.dropna(axis=1, how='all')  # drop nacolum
                self.ifdfsum = pd.concat([self.ifdfsum, sumdf])

            print(self.ifdfsum)

        if self.tgt2type in ['molname']:
            self.ifdf_filters = []
            self.ifdfsum = []

            nfrag = len(ifpidfs[0][0])
            print('nfrag', nfrag)
            # i:time, j:frag
            for j in range(nfrag):
                ifdf_filter = pd.DataFrame()
                ifdfsum = pd.DataFrame(columns=self.ifdfsumcolumn).astype(float)
                for i in range(len(ifpidfs)):
                    # ifdf_filter = ifdf_filter.append(ifpidfs[i][0][j])
                    # ifdfsum = ifdfsum.append(ifpidfs[i][1][j])

                    ifdf_filter = pd.concat(
                        [ifdf_filter, ifpidfs[i][0][j]])

                    sumdf = ifpidfs[1][1][j].to_frame().T
                    sumdf = sumdf.dropna(axis=1, how='all')
                    ifdfsum = pd.concat([ifdfsum, sumdf])

                self.ifdf_filters.append(ifdf_filter)
                self.ifdfsum.append(ifdfsum)

#                self.ifdf_filters = pd.concat([self.ifdf_filters,
#                                               ifdf_filter], ignore_index=True)
#                self.ifdfsum = pd.concat([self.ifdfsum, ifdfsum],
#                                         ignore_index=True)

            # print(self.ifdf_filters)
            # print(self.ifdfsum)
            # print('len ifdf_filters', len(self.ifdf_filters))
            # print('len ifdfsum', len(self.ifdfsum))

        p.close()
        print('read elapsed', time.time() - st)

        return

    def readsingleifie(self):
        '''read single mode

        read ifie single mode

        Args:
            None
        Returns:
            None(self.ifdf, self.pidf, self.pbifdf, self.pbpidf)
        '''

        print('## read single mode')
        self.logMethod = self.getlogmethod(self.tgtlogs)
        if self.logMethod == 'HF':
            self.icolumn = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE']
        elif self.logMethod == 'MP3':
            self.icolumn = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE',
                            'USER-MP2', 'MP3-IFIE', 'USER-MP3', 'PADE[2/1]' ]
        elif self.logMethod == 'CCPT':
            self.icolumn = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE',
                            'GRIMME-MP2', 'MP3-IFIE', 'GRIMME-MP3', 'MP4-IFIE']
        self.getpbflag(self.tgtlogs)

        # self.logMethod = 'MP2'
        if self.matrixtype != 'frags-frags' \
                and (self.logMethod == 'MP3' or self.logMethod == 'CCPT'):
            print('Error: ' + self.logMethod + ' mode for this mode is unsupported yet.')
            sys.exit()

        # pb python-based capture only in this version.
        if self.pbflag:
            ifie, pieda, momene, dimene, bsse = \
                self.read_ifiepieda(self.tgtlogs)
            df = self.getifiedf(ifie)
            self.ifdf = df

            pidf = self.getpiedadf(pieda)
            self.pidf = pidf

            # pbpieda
            pbifie, pbpieda, solvterm = self.read_pbifiepieda(self.tgtlogs)
            pbifdf = self.getifiedf(pbifie, solvterm)
            self.pbifdf = pbifdf

            pbpidf = self.getpbpiedadf(pbpieda)
            self.pbpidf = pbpidf

            # debug write
            # print(self.pbifdf)
            # sys.exit()
            return self

        if self.f90soflag is True:
            print("use fortran library")
            ifpidfs = self.read_ifpif90(self.tgtlogs)
            self.ifdf = ifpidfs[0]
            self.pidf = ifpidfs[1]
            # print(self.ifdf)
            # print(self.pidf)

        else:
            ifie, pieda, momene, dimene, bsse = self.read_ifiepieda(self.tgtlogs)
            df = self.getifiedf(ifie)
            self.ifdf = df

            pidf = self.getpiedadf(pieda)
            self.pidf = pidf

        return


    def readifiewrap(self, item1=None, item2=None, item3=None):
        '''read ifie and pieda
        Args:
            item1: tgtlogs
            item2: tgt1frag
            item3: tgt2type
        Returns:
        '''

        # param setup
        self.setupreadparm(item1, item2, item3)

        # multi mode (read and filter)
        if self.anlmode == 'multi':
            self.readmultiifie()
        # single mode
        else:
            self.readsingleifie()

        return self

    def getisdisp(self, tgtlog):
        f = open(tgtlog, 'r')
        for line in f:
            Items = line.split()
            if len(Items) <= 2:
                continue
            if Items[0:3] == ['Disp', '=', 'ON']:
                print('MP2-LRD single shot mode')
                return True
        return False

    def getlogmethod(self, tgtlog):
        f = open(tgtlog, 'r')
        for line in f:
            Items = line.split()
            if len(Items) <= 2:
                continue
            if Items[0:2] == ['Method', '=']:
                print('logMethod =', Items[2])
                return Items[2]

    def getpbflag(self, tgtlog):
        self.pbflag = False
        f = open(tgtlog, 'r')
        for line in f:
            Items = line.split()
            if len(Items) <= 2:
                continue
            if Items[0:2] == ['EFFECT', '=']:
                if Items[2] == 'ON':
                    self.pbflag = True
                    print('PB effect =', Items[2])
            if Items[0:3] == ['##', 'CHECK', 'AVAILABLE']:
                break

        return

    def setupreadparm(self, item1=None, item2=None, item3=None):
        '''setup read parameter
        Args:
            item1: tgt1
            item2: tgt2
            item3: tgt3
        Returns:
            tgtlogs: list of log files
            tgtpdbs: list of pdb files
            tgttimes: list of time steps
        '''

        tgtlogs = []
        tgtpdbs = []
        tgttimes = []

        if self.anlmode == 'multi' and self.tgt2type == 'molname':
            self.rpdbflag = True

        # multi mode
        if self.anlmode == 'multi':
            # item1: tgt1
            # item2: tgt2
            print('tgt2type:', self.tgt2type)

            # setup tgttimes, logs, and pdbs
            for i in range(self.start, self.end+1, self.interval):
                tgttimes.append(str(i).zfill(self.zp))
                tgtlogs.append(self.ilog_head + str(i).zfill(self.zp) + self.ilog_tail)
            print('tgtlogs', tgtlogs)

            # setup tgt1frag
            if item1 is not None:
                print('type', type(item1))
                if type(item1) == str:
                    if '-' in item1:
                        tgt = item1.split('-')
                        print('tgt', tgt)
                        self.tgt1frag = [i for i in range(int(tgt[0]), int(tgt[1]) + 1)]
                        if self.tgt1frag in self.tgt1frag:
                            del self.tgt1frag[self.tgt1frag.index(self.tgt1frag)]
                    else:
                        self.tgt1frag = list(map(int, item1.split(',')))
                        print(self.tgt1frag)
                else:
                        self.tgt1frag = item1

            self.tgtlogs = tgtlogs
            self.tgttimes = tgttimes

            # setup tgt2frag
            if self.tgt2type == 'frag':
                if item2 != None:
                    print('type', type(item2))
                    if item2[-1] == '-':
                        print('check tgt2 frags')

                        self.tgt2frag = []
                        for i in range(len(tgtlogs)):
                            self.resname_perfrag, tgtpdb = self.getlogorpdbfrag(self.tgtlogs[i])
                            nf = self.getlognf(tgtlogs[i], self.fragmode)
                            tgt = item2.split('-')[0]
                            print('tgt', tgt)
                            tgt2frag = [i for i in range(int(tgt), nf + 1)]
                            self.tgt2frag.append(tgt2frag)

                    elif type(item2) == str:
                        if '-' in item2:
                            tgt = item2.split('-')
                            print('tgt', tgt)
                            self.tgt2frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                            if self.tgt1frag in self.tgt2frag:
                                del self.tgt2frag[self.tgt2frag.index(self.tgt1frag)]
                        else:
                            self.tgt2frag = list(map(int, item2.split(',')))
                            print(self.tgt2frag)
                    else:
                            self.tgt2frag = [item2]
                    if type(self.tgt2frag) == list:
                        for dfrag in self.exceptfrag:
                            try:
                                del self.tgt2frag[self.tgt2frag.index(dfrag)]
                                print('- Info: del frag', dfrag, 'from tgt2')
                            except:
                                pass
                print('tgt1frag, tgt2frag', self.tgt1frag, self.tgt2frag)


            if self.tgt2type == 'molname':
                if item2 != None:
                    self.tgt2molname = item2
                print('tgt1frag, tgt2mol', self.tgt1frag, self.tgt2molname)


            if self.rpdbflag == True:
                nfs = []
                molnames_perrec = []
                resnamenonums_perfrag = []
                self.assignmolname = False
                for i in range(len(tgtlogs)):
                    self.resname_perfrag, tgtpdb = self.getlogorpdbfrag(self.tgtlogs[i])
                    tgtpdbs.append(tgtpdb)
#                     if self.fragmode == 'auto':
#                         print('Error: auto fragment in mol mode is not suppoted yet.')
#                         sys.exit()
                    nf = self.getlognf(tgtlogs[i], self.fragmode)
                    nfs.append(nf)
                    molnames_perrec.append(self.resnames)
                    resnamenonums_perfrag.append(self.resnamenonum_perfrag)
                self.nfs = nfs
                self.molnames_perrec  = molnames_perrec
                self.resnamenonums_perfrag = resnamenonums_perfrag
                self.tgtpdbs = tgtpdbs

        # single mode
        # item1 log
        # item2 tgt1
        # item3 tgt2
        else:
            if item1 is not None:
                self.tgtlogs = item1
            if item2 is not None:
                if self.anlmode == 'mol' and self.selecttype == 'molid':
                    self.tgtmolid = int(item2)
                else:
                    print('item2', item2)
                    if '-' in item2:
                        tgt = item2.split('-')
                        print('tgt1', tgt)
                        self.tgt1frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]

#                     if type(eval(item2)) != list:
#                         self.tgt1frag = [item2]
#                     else:
#                         self.tgt1frag = eval(item2)

                    else:
                        self.tgt1frag = list(map(int, item2.split(',')))
                        print(self.tgt1frag)

            if self.anlmode == 'fraginmol' or self.anlmode == 'mol':
                if type(self.tgtmolid) == str:
                    if '-' in self.tgtmolid:
                        tgt = self.tgtmolid.split('-')
                        print('tgtmol', tgt)
                        self.tgtmolid = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                    else:
                        self.tgtmolid = [eval(self.tgtmolid)]

                elif type(self.tgtmolid) == int:
                    self.tgtmolid = [self.tgtmolid]

            if self.tgt2type == 'frag' or self.tgt2type == 'dist':
                if item2 is not None:
                    # print('type tgt2', type(item2))
                    if type(item2) == str:
                        if '-' in item2:
                            tgt = item2.split('-')
                            print('tgt1', tgt)
                            self.tgt1frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                        else:
                            print('check2')
                            self.tgt1frag = list(map(int, item2.split(',')))
                            print(self.tgt1frag)

#                             if type(eval(item2)) != list:
#                                 self.tgt1frag = [eval(item2)]
#                             else:
#                                 self.tgt1frag = eval(item2)
                    else:
                        self.tgt1frag = [item2]

            if self.tgt2type == 'frag':
                if item3 is not None:
                    # print('type tgt2', type(item3))
                    if type(item3) == str:
                        if '-' in item3:
                            tgt = item3.split('-')
                            print('tgt2', tgt)
                            self.tgt2frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                        else:
                            self.tgt2frag = list(map(int, item3.split(',')))
                            print(self.tgt2frag)
                            # self.tgt2frag = [eval(item3)]
                    else:
                        self.tgt2frag = [item3]

        if self.matrixtype == 'frags-frags':
            print(self.tgt1frag, self.tgt2frag)
#             dp = set(self.tgt1frag) & set(self.tgt2frag)
#             if len(dp) != 0:
#                 print('Error! tgt1 and tgt2 is duplicate')
#                 sys.exit()

        if type(self.tgt2frag) == list:
            for dfrag in self.exceptfrag:
                try:
                    del self.tgt2frag[self.tgt2frag.index(dfrag)]
                    print('- Info: del frag', dfrag, 'from tgt2')
                except:
                    pass

        # PIEDA
        if self.abinit_ver == 'rev16' or self.abinit_ver == 'rev17':
            self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'Solv(ES)', 'DI(MP2)', 'q(I=>J)']

        # get resname reference
        if self.addresinfo == True or self.rpdbflag == True or self.anlmode == 'fraginmol':
            print('\n## read reference resname')
            if type(self.tgtlogs) == list:
                for i in range(len(self.tgtlogs)):
                    self.resname_perfrag, self.tgtpdb = self.getlogorpdbfrag(self.tgtlogs[i])
                    if len(self.resname_perfrag) != 0:
                        break
                    else:
                        print('cannot read frag data:', self.tgtlogs[i])
            else:
                self.resname_perfrag, self.tgtpdb = self.getlogorpdbfrag(self.tgtlogs)

        # print section
        print ('\n## input summary')
        try:
            print('- tgtlogs:', self.tgtlogs)
            print('- anlmode:', self.anlmode)
            print('- tgtfrag:', self.tgt1frag, self.tgt2frag)
            print('- anlmode:' ,self.anlmode)
            print('- fragmode:', self.fragmode)
            print('- NP:', self.pynp)
            print('- addresflag', self.addresinfo)
        except:
            pass
        print('## input summary end\n')

        ### read fraginfo section
#         frags = []
#         if self.fragmode != 'manual':
#             frags = self.read_fraginfo(self.tgtlogs)
#             # print('frags', frags)
#
#         if self.fragmode == 'hybrid':
#             getf = frags.pop(hyfrag-1)
#             for i in range(self.hynum):
#                 frags.append(getf)
#             # print('frags', frags)
#
#         self.frags = frags

    ## filter section
    def filterifiewrap(self, dist=None, myifdf=None, mypidf=None, is_pb=False):

        tgt2type = self.tgt2type
        if dist is not None:
            self.dist = dist
        # frag mode
        if self.anlmode == 'frag':
            if self.matrixtype == 'frags-frags':
                self.gettgtdf_n2ffmatrix(myifdf)
                self.gettgtpidf_n2ffmatrix(mypidf, is_pb)

            else:
                if tgt2type in ['dist', 'dimer-es']:
                    print(self.ifdf)
                    tgtdf, ifdf_filter = self.gettgtdf_fd(self.ifdf)
                    self.ifdf_filter = pd.merge(ifdf_filter,
                                                self.pidf, on=['I', 'J'], how='left')
                    print(self.ifdf_filter)

                    if self.pbflag is True:
                        pbtgtdf, pbifdf_filter = self.gettgtdf_fd(self.pbifdf)
                        self.pbifdf_filter = pd.merge(
                            pbifdf_filter, self.pbpidf, on=['I', 'J'], how='left')
                        print(self.pbifdf_filter)

                elif tgt2type == 'frag':
                    self.ifdf_filters = []
                    self.pbifdf_filters = []

                    if type(self.tgt1frag) == int:
                        self.tgt1frag == [self.tgt1frag]
                    if type(self.tgt2frag) == int:
                        self.tgt1frag == [self.tgt2frag]
                    for tgt1 in self.tgt1frag:
                        ifdf_filter = self.gettgtdf_ffs(self. ifdf, tgt1, self.tgt2frag)
                        self.ifdf_filters.append(
                            pd.merge(ifdf_filter, self.pidf, on=['I','J'], how='left'))
                        if self.pbflag is True:
                            pbifdf_filter = self.gettgtdf_ffs(
                                self.pbifdf, tgt1, self.tgt2frag)
                            self.pbifdf_filters.append(
                                pd.merge(pbifdf_filter, self.pbpidf, on=['I', 'J'], how='left'))
                    print(self.ifdf_filters[0].head())
                    if self.pbflag is True:
                        print(self.pbifdf_filters[0].head())

        # mol-mol mode
        if self.anlmode == 'mol':
            # ifie
            df = self.ifdf
#             if self.selecttype == 'fragid':
#                 tgtmolfrags = self.getmolfrags(self.tgt1frag[0], df)
#                 print('target-frags:', molfrags)
            if self.selecttype == 'molid':
                nf = self.getlognf(self.tgtlogs, self.fragmode)
                molfragss = self.getallmolfrags(self.tgtlogs, df, nf)
                print('frags_permol\n', molfragss)
                tgtmolfrags = []
                for tgtmolid in self.tgtmolid:
                    tgtmolfrags.append(molfragss[tgtmolid-1])
            elif self.selecttype == 'molname':
                sys.exit()

            self.tgtmolfrags = tgtmolfrags
            print('self.tgtmolfrags:', self.tgtmolfrags)

            # IFIE and pieda
            ifdf_frag_mols = pd.DataFrame()
            ifdfmol_mols = pd.DataFrame(columns=['I', 'J'] + self.ifdfsumcolumn)
            ifdfmolsums = pd.DataFrame(columns=self.ifdfsumcolumn).astype(float)

            for i in range(len(self.tgtmolfrags)):
                contactmolfrags, ifdf_frag_mol, ifdfmol_mol, ifdfmolsum = \
                    self.getifiesummol(df, tgtmolfrags[i], self.tgtmolid[i])

                ifdf_frag_mol = pd.concat(ifdf_frag_mol)

                # self.contactmolfrags = contactmolfrags
                print('ifdf_frag_mol', ifdf_frag_mol)  # list
                print('ifdfmol_mol', ifdfmol_mol)
                print('ifdfmolsum', ifdfmolsum)

                print('ifdf_frag_mol', type(ifdf_frag_mol))
                print('ifdfmol_mol', type(ifdfmol_mol))
                print('ifdfmolsum', type(ifdfmolsum))

                # ifdf_frag_mols = ifdf_frag_mols.append(ifdf_frag_mol)
                # ifdfmol_mols = ifdfmol_mols.append(ifdfmol_mol)
                # ifdfmolsums = ifdfmolsums.append(ifdfmolsum)

                 # ifdf_frag_molsにifdf_frag_molを結合
                ifdf_frag_mols = pd.concat([ifdf_frag_mols, ifdf_frag_mol])
                print('ifdf_frag_mols', type(ifdf_frag_mols))

                # ifdfmol_molsにifdfmol_molを結合
                ifdfmol_mols = pd.concat([ifdfmol_mols, ifdfmol_mol])
                # ifdfmolsumsにifdfmolsumを結合
                ifdfmolsum = ifdfmolsum.to_frame().T
                ifdfmolsum = ifdfmolsum.dropna(axis=1, how='all')
                ifdfmolsums = pd.concat([ifdfmolsums, ifdfmolsum])

            self.ifdf_frag_mols = ifdf_frag_mols
            print('ifdf_frag_mols', ifdf_frag_mols)
            self.ifdfmol_mols = ifdfmol_mols
            print('ifdfmol_mols', ifdfmol_mols)
            self.ifdfmolsums = ifdfmolsums
            print('ifdfmolsums', ifdfmolsums)

        # fraginmol mode
        if self.anlmode == 'fraginmol':
            ifdf_filters = pd.DataFrame()
            ifdfsums = pd.DataFrame(columns=self.ifdfsumcolumn).astype(float)

            df = self.ifdf
            tgt1_lofrag = self.tgt1_lofrag
            tgt2_lofrag = self.tgt2_lofrag
            tgt2molname = self.tgt2molname
            nf = self.getlognf(self.tgtlogs, self.fragmode)
            print('nf', nf)
            molfragss = self.getallmolfrags(self.tgtlogs, df, nf)
            print('molfragss', molfragss)
            print('len_molfragss', len(molfragss))
            print('resnames', self.resnames)
            print('self.resnames', len(self.resnames))

            print('resnames_perfrag', self.resnamenonum_perfrag)
            print('len self.resnames', len(self.resnamenonum_perfrag))

            rname_perfrag = self.resnamenonum_perfrag

            tgtmolfrags = []

            rname_permol = copy.deepcopy(molfragss)
            for j in range(len(rname_permol)):
                # if type(molfragss[j]) == list:
                for k in range(len(molfragss[j])):
                    num = copy.deepcopy(molfragss[j][k]) - 1
                    rname_permol[j][k] = rname_perfrag[num]

            tgt2_glofrags = []
            # print('resnames', self.resnames)
            for i in range(len(rname_permol)):
                if tgt2molname in rname_permol[i]:
                    # print(self.resnames[i], tgt2molname)
                    tgt2frag = molfragss[i][tgt2_lofrag - 1]
                    tgt2_glofrags.append(tgt2frag)
            print('tgt2_glofrags', tgt2_glofrags)

            # tgtmol loop
            for tgtmol in self.tgtmolid:
                tgtmol = tgtmol - 1

                tgt1_glofrag = molfragss[tgtmol][tgt1_lofrag - 1]
                print('tgt1glofrag', tgt1_glofrag)
                print('centermolfrag:', tgt1_glofrag)
                print('tgt2molname', tgt2molname)
                tgtdf = df[df['I'] == tgt1_glofrag]
                tgtdf = tgtdf.append(df[df['J'] == tgt1_glofrag])
                tgtdf = tgtdf[tgtdf['DIST'] < self.dist]
                tgtdf = tgtdf[tgtdf['DIST'] != 0.0]

                # ifdf_filters = pd.DataFrame()
                tgtdf_filter = tgtdf[(tgtdf['I'].isin(tgt2_glofrags)) | (tgtdf['J'].isin(tgt2_glofrags))]
                # print('ifdf_filters\n', tgtdf_filter)

                # PIEDA
                ifdf_filter = pd.merge(tgtdf_filter, self.pidf, on=['I', 'J'], how='left')
                ifdf_filters = ifdf_filters.append(ifdf_filter)

                # print(ifdf_filter)
                HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, \
                    JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, \
                    q_sum = self.getsumdf(ifdf_filter)

                ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum,
                                     GRIMME_sum, JUNG_sum, HILL_sum, ES_sum,
                                     EX_sum, CT_sum, DI_sum, q_sum],
                                    index=self.ifdfsumcolumn, name='mol' + str(tgtmol+1))
                ifdfsums = ifdfsums.append(ifdfsum)

            # self.tgt1_glofrag = tgt1_glofrags
            # self.tgt2_glofrags = tgt2_glofrags
            self.ifdf_filters = ifdf_filters
            self.ifdfsums = ifdfsums

        return self

    def writecsvwrap(self, head=None, word='', pbwrite=False):
        '''writecsvwrap

        writecsv section

        Args:
            head (str): head of csv file name
            word (str): word of csv file name
            pbwrite (bool): write csv file or not

        Returns:
            self (object)
        '''

        print('## Write Section')
        path = 'csv'
        tgt2type = self.tgt2type
        if os.path.exists('csv') is False:
            os.mkdir('csv')

        if self.anlmode == 'frag':
            if head is None:
                head = os.path.splitext(self.tgtlogs)[0].split('/')[-1]

            if self.matrixtype == 'frags-frags':
                if self.logMethod == 'MP2':
                    datadfs = [
                                self.esdf,
                                self.exdf,
                                self.ctdf,
                                self.hfdf,
                                self.mp2corrdf,
                                self.prmp2corrdf,
                                self.mp2tdf,
                                self.prmp2tdf,
                                self.distdf
                            ]
                    names = [
                                'ES',
                                'EX',
                                'CT',
                                'HF',
                                'MP2corr',
                                'PRMP2corr',
                                'MP2total',
                                'PRMP2total',
                                'Distance'
                            ]
                elif self.logMethod == 'MP3':
                    datadfs = [
                                self.esdf,
                                self.exdf,
                                self.ctdf,
                                self.hfdf,
                                self.mp2corrdf,
                                self.mp3corrdf,
                                self.mp25corrdf,
                                self.usermp3corrdf,
                                self.mp2tdf,
                                self.mp3tdf,
                                self.mp25tdf,
                                self.usermp3tdf,
                                self.distdf
                            ]
                    names = [
                                'ES',
                                'EX',
                                'CT',
                                'HF',
                                'MP2corr',
                                'MP3corr',
                                'MP25corr',
                                'USER-MP3corr',
                                'MP2total',
                                'MP3total',
                                'MP25total',
                                'USER-MP3total',
                                'Distance'
                            ]
                elif self.logMethod == 'CCPT':
                    datadfs = [
                                self.esdf,
                                self.exdf,
                                self.ctdf,
                                self.hfdf,
                                self.mp2corrdf,
                                self.mp3corrdf,
                                self.mp4corrdf,
                                self.mp25corrdf,
                                self.mp35corrdf,
                                self.mp2tdf,
                                self.mp3tdf,
                                self.mp4tdf,
                                self.mp25tdf,
                                self.mp35tdf,
                                self.distdf
                           ]
                    names = [
                                'ES',
                                'EX',
                                'CT',
                                'HF',
                                'MP2corr',
                                'MP3corr',
                                'MP4corr',
                                'MP25corr',
                                'MP35corr',
                                'MP2total',
                                'MP3total',
                                'MP4total',
                                'MP25total',
                                'MP35total',
                                'Distance'
                            ]
                if self.logMethod == 'HF':
                    datadfs = [
                                self.esdf,
                                self.exdf,
                                self.ctdf,
                                self.hfdf,
                                self.distdf
                            ]
                    names = [
                                'ES',
                                'EX',
                                'CT',
                                'HF',
                                'Distance'
                            ]
                if self.logMethod == 'HF+D':
                    datadfs = [
                                self.esdf,
                                self.exdf,
                                self.ctdf,
                                self.hfdf,
                                self.mp2corrdf,
                                self.mp2tdf,
                                self.distdf
                            ]
                    names = [
                                'ES',
                                'EX',
                                'CT',
                                'HF',
                                'DILRDcorr',
                                'LRDtotal',
                                'Distance'
                            ]

                tgt1str = str(self.tgt1frag[0]) + '-'  + str(self.tgt1frag[-1])
                tgt2str = str(self.tgt2frag[0]) + '-'  + str(self.tgt2frag[-1])
                for i in range(len(datadfs)):
                    # rename index-columns
                    if self.addresinfo:
                        datadfs[i].rename(index=lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')',
                                          columns=lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')', inplace=True)

                    ocsv = head + '_frag' + str(tgt1str) + '-frag' + str(tgt2str) + '-' + word + names[i] + '-ffmatrix.csv'
                    datadfs[i].T.to_csv(path + '/' + ocsv, float_format='%.6f')
                    print(path + '/' + ocsv + ' was created.')

                if pbwrite:
                    if self.addresinfo:
                        self.solvesdf.rename(index=lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')',
                                             columns=lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')', inplace=True)

                    ocsv = head + '_frag' + str(tgt1str) + '-frag' + str(tgt2str) + '-' + word + 'SolvES-ffmatrix.csv'
                    self.solvesdf.T.to_csv(path + '/' + ocsv, float_format='%.6f')
                    print(path + '/' + ocsv + ' was created.')

            else:
                if self.tgt2type in ['dist', 'dimer-es']:
                    tgtid = self.tgt1frag[0]
                    try:
                        ohead = head + '-' + str(tgtid) + '-' + frags[tgtid - 1]
                    except:
                        ohead = head + '-' + str(tgtid)

                    if self.tgt2type == 'dist':
                        oifie = ohead + '-ifie_' + 'dist' + str(self.dist) + '.csv'
                    else:
                        oifie = ohead + '-ifie_dimer-es-false.csv'


                    if self.addresinfo:
                        for i in range(1, len(self.resname_perfrag)+1):
                            val1 = i
                            val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                            self.ifdf_filter.I = self.ifdf_filter.I.replace(val1, val2)
                            self.ifdf_filter.J = self.ifdf_filter.J.replace(val1, val2)

                        print(self.ifdf_filter.I)

                    self.ifdf_filter.to_csv(path + '/' + oifie, float_format='%.6f')
                    print(path + '/' + oifie, 'was generated.')

                    if self.pbflag:
                        try:
                            ohead = head + '-' + str(tgtid) + '-' + frags[tgtid - 1]
                        except:
                            ohead = head + '-' + str(tgtid)

                    if self.tgt2type == 'dist':
                        oifie = ohead + '-pbifie_' + 'dist' + str(self.dist) + '.csv'
                    else:
                        oifie = ohead + '-pbifie_dimer-es-false.csv'

                        if self.addresinfo:
                            for i in range(1, len(self.resname_perfrag)+1):
                                val1 = i
                                val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                                self.pbifdf_filter.I = self.pbifdf_filter.I.replace(val1, val2)
                                self.pbifdf_filter.J = self.pbifdf_filter.J.replace(val1, val2)

                            print(self.pbifdf_filter.I)

                        self.pbifdf_filter.to_csv(path + '/' + oifie, float_format='%.6f')

                        print(path + '/' + oifie, 'was generated.')

                elif self.tgt2type == 'frag':

                    tgt1s = self.tgt1frag
                    count = 0
                    for tgt1 in tgt1s:
                        ohead = head + '-frag' + str(tgt1) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1])
                        oifie = ohead + '-ifie.csv'

                        ifdf_filter = self.ifdf_filters[count]

                        if self.addresinfo:
                            for i in range(1, len(self.resname_perfrag)+1):
                                val1 = i
                                val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                                ifdf_filter.I = ifdf_filter.I.replace(val1, val2)
                                ifdf_filter.J = ifdf_filter.J.replace(val1, val2)

                        ifdf_filter.to_csv(path + '/' + oifie, float_format='%.6f')
                        print(path + '/' + oifie, 'was generated.')
                        count += 1
                    # N:1 sheet
                    for j in range(len(self.ifdf_filters)):
                        if j == 0:
                            ifdf_fil_n1 = self.ifdf_filters[j].fillna(0.0)
                            # pidf_fil_n1 = self.pidf_filters[j]
                        else:
                            ifdf_fil_n1 += self.ifdf_filters[j].fillna(0.0).values

                            print(self.ifdf_filters[j])
                            # pidf_fil_n1 += self.pidf_filters[j].values

                    ohead = head + '-frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1]) + 'n-1sum'
                    oifie = ohead + '-ifie.csv'

                    ifdf_fil_n1 = ifdf_fil_n1.reset_index(drop=True)

                    if self.addresinfo:
                        ifdf_fil_n1["I"] = self.resname_perfrag[int(self.tgt1frag[0])-1] + '(' + str(self.tgt1frag[0]) + ')' + '-' + \
                            self.resname_perfrag[int(self.tgt1frag[-1])-1] + '(' + str(self.tgt1frag[-1]) + ')'
                        ifdf_fil_n1["J"] = self.tgt2frag

                        for i in range(1, len(self.resname_perfrag)+1):
                            val1 = i
                            val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                            ifdf_fil_n1.J = ifdf_fil_n1.J.replace(val1, val2)
                    else:
                        # print(self.tgt1frag)
                        ifdf_fil_n1["I"] = str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1])
                        ifdf_fil_n1["J"] = self.tgt2frag

                    del ifdf_fil_n1['DIMER-ES']
                    del ifdf_fil_n1['DIST']

                    ifdf_fil_n1.to_csv(path + '/' + oifie, float_format='%.6f')

                    print(path + '/' + oifie, 'was generated.')

                    # pb
                    if self.pbflag:
                        count = 0
                        for tgt1 in tgt1s:
                            ohead = head + '-frag' + str(tgt1) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1])
                            oifie = ohead + '-pbifie.csv'

                            pbifdf_filter = self.pbifdf_filters[count]

                            if self.addresinfo:
                                for i in range(1, len(self.resname_perfrag)+1):
                                    val1 = i
                                    val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                                    pbifdf_filter.I = pbifdf_filter.I.replace(val1, val2)
                                    pbifdf_filter.J = pbifdf_filter.J.replace(val1, val2)

                            pbifdf_filter.to_csv(path + '/' + oifie, float_format='%.6f')
                            print(path + '/' + oifie, 'was generated.')
                            count += 1

                        # N:1 sheet
                        for j in range(len(self.pbifdf_filters)):
                            if j == 0:
                                pbifdf_fil_n1 = self.pbifdf_filters[j]
                                # pidf_fil_n1 = self.pidf_filters[j]
                            else:
                                pbifdf_fil_n1 += self.pbifdf_filters[j].values
                                # pidf_fil_n1 += self.pidf_filters[j].values

                        ohead = head + '-frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1]) + 'n-1sum'
                        oifie = ohead + '-ifie.csv'

                        pbifdf_fil_n1 = pbifdf_fil_n1.reset_index(drop=True)

                        if self.addresinfo:
                            pbifdf_fil_n1["I"] = self.resname_perfrag[int(self.tgt1frag[0])-1] + '(' + str(self.tgt1frag[0]) + ')' + '-' + \
                                self.resname_perfrag[int(self.tgt1frag[-1])-1] + '(' + str(self.tgt1frag[-1]) + ')'
                            pbifdf_fil_n1["J"] = self.tgt2frag

                            for i in range(1, len(self.resname_perfrag)+1):
                                val1 = i
                                val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                                pbifdf_fil_n1.J = pbifdf_fil_n1.J.replace(val1, val2)
                        else:
                            pbifdf_fil_n1["I"] = self.tgt1frag[0] + '-' + self.tgt1frag[-1]
                            pbifdf_fil_n1["J"] = self.tgt2frag

                        del pbifdf_fil_n1['DIMER-ES']
                        del pbifdf_fil_n1['DIST']

                        pbifdf_fil_n1.to_csv(path + '/' + oifie, float_format='%.6f')

                        print(path + '/' + oifie, 'was generated.')

        if self.anlmode == 'multi':
            head = self.ilog_head.split('/')[-1]
            if tgt2type == 'frag':

                if self.matrixtype == 'times-frags':

                    if self.is_disp:
                        datadfs = [
                                    self.esdfs,
                                    self.exdfs,
                                    self.ctdfs,
                                    self.hfdfs,
                                    self.mp2corrdfs,
                                    self.prmp2corrdfs,
                                    self.mp2tdfs,
                                    self.prmp2tdfs,
                                    self.distdfs,
                                    self.didfs,
                                    self.erestdfs
                                ]
                        names = [
                                    'ES',
                                    'EX',
                                    'CT',
                                    'HF',
                                    'MP2corr',
                                    'PRMP2corr',
                                    'MP2total',
                                    'PRMP2total',
                                    'distance',
                                    'LRD',
                                    'PRMP2Erest',
                                ]
                    else:
                        datadfs = [
                                    self.esdfs,
                                    self.exdfs,
                                    self.ctdfs,
                                    self.hfdfs,
                                    self.mp2corrdfs,
                                    self.prmp2corrdfs,
                                    self.mp2tdfs,
                                    self.prmp2tdfs,
                                    self.distdfs
                                ]
                        names = [
                                    'ES',
                                    'EX',
                                    'CT',
                                    'HF',
                                    'MP2corr',
                                    'PRMP2corr',
                                    'MP2total',
                                    'PRMP2total',
                                    'distance',
                                ]

                    tgt2str = str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1])
                    if self.depth(self.tgt2frag) >= 2:
                        tgt2str = str(self.tgt2frag[0][0]) + '-end'

                    print('--- out files ---')
                    # i:type, j:tgt1frag
                    for i in range(len(datadfs)):
                        sum2df = pd.DataFrame(index=datadfs[i][0].columns)
                        for j in range(len(datadfs[i])):
                            tgt1frag = self.tgt1frag[j]
                            ocsv = head + '_frag' + str(tgt1frag) + '-frag' + str(tgt2str) + '-' + names[i] + '-tfmatrix.csv'
                            if self.addresinfo:
                                datadfs[i][j].rename(index = lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')', inplace=True)
                            datadfs[i][j].T.to_csv(path + '/' + ocsv, float_format='%.6f')
                            print(path + '/' + ocsv)

                            # gen tgtfrag1 sum matrix
                            if j == 0:
                                sum1df = datadfs[i][j]
                            else:
                                sum1df = sum1df + datadfs[i][j]

                            # gen tgtfrag2 sum matrix
                            sum2dfbuf = datadfs[i][j].sum()
                            sum2dfbuf.name = tgt1frag
                            sum2df = pd.concat([sum2df, sum2dfbuf], axis=1)

                        osum1csv = head + '_frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(tgt2str) + '-' + names[i] + '-sumtfmatrix.csv'
                        sum1df.T.to_csv(path + '/' + osum1csv, float_format='%.6f')
                        print(path + '/' + osum1csv)

                        osum2csv = head + '_frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(tgt2str) + '-' + names[i] + '-sum2matrix.csv'
                        sum2df.to_csv(path + '/' + osum2csv, float_format='%.6f')
                        print(path + '/' + osum2csv)

                else:
                    tgt1frag = self.tgt1frag[0]
                    tgt2frag = self.tgt2frag[0]

                    if self.addresinfo:
                        for i in range(1, len(self.resname_perfrag)+1):
                            val1 = i
                            val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                            self.ifdf_filters.I = self.ifdf_filters.I.replace(val1, val2)
                            self.ifdf_filters.J = self.ifdf_filters.J.replace(val1, val2)

                    oifie = 'frag' + str(tgt1frag) + '-frag' + str(tgt2frag) + '-ifie.csv'
                    self.ifdf_filters.to_csv(path + '/' + oifie, float_format='%.6f')
                    print(path + '/' + oifie, 'was created.')

            if tgt2type in ['dist', 'dimer-es']:

                tgt1frag = self.tgt1frag[0]

                if tgt2type == 'molname':
                    tgt2molname = self.tgt2molname
                    oifie = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-dist' + str(self.dist) + '-ifiesum.csv'
                    oifiedt = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-dist' + str(self.dist) + '-ifiedt.csv'

                if tgt2type == 'dist':
                    tgt2dist = self.dist
                    oifie = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-ifiesum.csv'
                    oifiedt = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-ifiedt.csv'

                if tgt2type == 'dimer-es':
                    oifie = 'frag' + str(tgt1frag) + '-dimer-es-false-ifiesum.csv'
                    oifiedt = 'frag' + str(tgt1frag) + '-dimer-es-false-ifiedt.csv'

                if self.addresinfo:
                    for i in range(1, len(self.resname_perfrag)+1):
                        val1 = i
                        val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                        self.ifdf_filters.I = self.ifdf_filters.I.replace(val1, val2)
                        self.ifdf_filters.J = self.ifdf_filters.J.replace(val1, val2)

                self.ifdfsum.to_csv(path + '/' + oifie, float_format='%.6f')
                self.ifdf_filters.to_csv(path + '/' + oifiedt, float_format='%.6f')

                print(path + '/' + oifie)
                print(path + '/' + oifiedt, 'was created.')

            if tgt2type in ['molname']:

                for j in range(len(self.tgt1frag)):
                    tgt1frag = self.tgt1frag[j]

                    if tgt2type == 'molname':
                        tgt2molname = self.tgt2molname
                        oifie = 'frag' + str(tgt1frag) + '-' + \
                            str(tgt2molname) + '-dist' + \
                            str(self.dist) + '-ifiesum.csv'
                        oifiedt = 'frag' + str(tgt1frag) + \
                            '-' + str(tgt2molname) + '-dist' + \
                            str(self.dist) + '-ifiedt.csv'

                    if self.addresinfo:
                        for i in range(1, len(self.resname_perfrag)+1):
                            val1 = i
                            val2 = self.resname_perfrag[i-1] + '(' + \
                                str(val1) + ')'
                            self.ifdf_filters[j].I = \
                                self.ifdf_filters[j].I.replace(val1, val2)
                            self.ifdf_filters[j].J = \
                                self.ifdf_filters[j].J.replace(val1, val2)

                    self.ifdfsum[j].to_csv(path + '/' + oifie, float_format='%.6f')
                    self.ifdf_filters[j].to_csv(path + '/' + oifiedt, float_format='%.6f')

                    print(path + '/' + oifie)
                    print(path + '/' + oifiedt, 'was created.')

        if self.anlmode == 'mol':
            if head is None:
                head = os.path.splitext(self.tgtlogs)[0].split('/')[-1]

            dist = self.dist
            selecttype = self.selecttype
            # piedamol_mol = self.piedamol_mol
            ifdf_frag_mols = self.ifdf_frag_mols
            ifdfmol_mols = self.ifdfmol_mols

            # pieda_frag_mols = self.pieda_frag_mols
            if selecttype == 'molid':
                tgtid = self.tgtmolid
            else:
                tgtid = self.tgt1frag[0]

            idstr = str(tgtid[0]) + '-' + str(tgtid[-1])
            ilogdtname = path + '/' + head + '_ifie-fragmol-' + selecttype + \
                idstr + 'dist' + str(dist) + '.csv'
            imolname = path + '/' + head + '_ifiemol-mol-' + selecttype + \
                idstr + 'dist' + str(dist) + '.csv'
            isumname = path + '/' + head + '_ifiesummol-mol-' + selecttype + \
                idstr + 'dist' + str(dist) + '.csv'

#             ifdf_frag_molsdt = pd.DataFrame()
#             pd.set_option('display.width', 500)
#             for ifdf_frag_mol in ifdf_frag_mols:
#                 ifdf_frag_molsdt = ifdf_frag_molsdt.append(ifdf_frag_mol)

            if self.addresinfo:
                for i in range(1, len(self.resname_perfrag)+1):
                    val1 = i
                    val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                    ifdf_frag_mols.I = ifdf_frag_mols.I.replace(val1, val2)
                    ifdf_frag_mols.J = ifdf_frag_mols.J.replace(val1, val2)

            # print(ifdf_frag_molsdt, file=ilogdt)
            ifdf_frag_mols.to_csv(ilogdtname, float_format='%.6f')
            ifdfmol_mols.to_csv(imolname, float_format='%.6f')
            self.ifdfmolsums.to_csv(isumname, float_format='%.6f')

            print('---out---')
            print(ilogdtname)
            print(imolname)
            print(isumname)

        if self.anlmode == 'fraginmol':
            if head is None:
                head = os.path.splitext(self.tgtlogs)[0].split('/')[-1]

            ohead = head + '-' 'tgt1frag' + str(self.tgt1_lofrag) + '-mol' + \
                str(self.tgt2molname) + 'frag' + str(self.tgt2_lofrag)

            if self.addresinfo is True:
                for i in range(1, len(self.resname_perfrag)+1):
                    val1 = i
                    val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                    self.ifdf_filters.I = \
                        self.ifdf_filters.I.replace(val1, val2)
                    self.ifdf_filters.J = \
                        self.ifdf_filters.J.replace(val1, val2)
                    # self.ifdf.I = self.ifdf.I.replace(val1, val2)
                    # self.ifdf.J = self.ifdf.J.replace(val1, val2)

            # self.ifdf.to_csv(path + '/' + head + '-ifie.csv')
            oifie = path + '/' + ohead + '-ifie_' + 'dist' + \
                str(self.dist) + '.csv'
            oifiesum = path + '/' + ohead + '-ifiesum_' + 'dist' + \
                str(self.dist) + '.csv'
            self.ifdf_filters.to_csv(oifie, float_format='%.6f')
            self.ifdfsums.to_csv(oifiesum, float_format='%.6f')
            print(oifie, 'was generated.')
            print(oifiesum, 'was generated.')
