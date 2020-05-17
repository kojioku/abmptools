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
import csv
import pdb_io as pdio
import time
from ctypes import *
# import setparam as sp
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
        self.ajf_method = 'HF'
        self.ajf_basis_set = '6-31Gdag'
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
        self.tgtpos =[]
        self.icolumn = ['I', 'J', 'DIST', 'DIMER-ES', 'HF-IFIE', 'MP2-IFIE', 'PR-TYPE1', 'GRIMME', 'JUNG', 'HILL']
        self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']

        self.anlmode= 'frag' #frag, 'mol', 'fraginmol', 'ff-multi'
        self.fragmode = 'auto'  #'hybrid', 'auto', 'manual'
        self.dist = 10.0
        self.tgtfrag = None

        self.rpdbflag = False
        self.pdbname = None   # 'iss2-spg2-ok20200130opt-for_abmp.pdb'

        # -- for mol mode or ff-multi mode--
        self.tgt2type = 'frag' #frag: mol-frag, mol: mol-mol

        # mol - mol mode
        self.selecttype = 'molid'
        self.tgtmolid = None

        # -- fraginmol mode --
        self.tgt1_lofrag = None
        self.tgt2_lofrag = None
        self.tgt2molname = None
        # ------ ff-multi mode ------
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
        self.hyfrag = None #320
        self.hynum = None

        self.pynp = 3

        mydir = os.path.dirname(os.path.abspath(__file__))
        self.f90sofile = mydir + '/f90/bin/readifiepiedalib.so'
        self.f90soflag = True

        pass


    def read_fraginfo(self, fname):
        frags = []
        count = 0
        text = open(fname, "r").readlines()
        flag = False
        for i in range(len(text)):
            itemList = text[i][:-1].split()
            if len(itemList) < 2:
                continue
            if itemList[1] == 'AUTOMATIC' or itemList[1] == 'HYBRID':
                flag = True
                continue
            if itemList[1] == 'MANUAL':
                manflag = True
            if itemList[1] == 'system':
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
        neighbors_list = []
        # print(idxs)
        for idx in idxs:
            tgtdf = df[df['I'] == idx]
            tgtdf = tgtdf.append(df[df['J'] == idx])
            tgtdf_zero = tgtdf[tgtdf['DIST'] == 0.]
            # print(tgtdf_zero)
            neighbor_i = [index for index, row in tgtdf_zero.groupby("I")]
            neighbor_j = [index for index, row in tgtdf_zero.groupby("J")]
            neighbors = set(neighbor_i + neighbor_j)
            # print('connect_idx', neighbors)
            neighbors_list.append(neighbors)
        neighbors_flat = list(itertools.chain.from_iterable(neighbors_list))
        # print(idx)
        # print('neighbors_flat', neighbors_flat)

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
        return nf

    def getifiedf(self, ifie):
        df = pd.DataFrame(ifie, columns=self.icolumn)
        df['I'] = df['I'].astype(int)
        df['J'] = df['J'].astype(int)
        df['DIST'] = df['DIST'].astype(float)
        df['HF-IFIE'] = df['HF-IFIE'].astype(float) * 627.5095
        df['MP2-IFIE'] = df['MP2-IFIE'].astype(float) * 627.5095
        df['PR-TYPE1'] = df['PR-TYPE1'].astype(float) * 627.5095
        df['GRIMME'] = df['GRIMME'].astype(float) * 627.5095
        df['JUNG'] = df['JUNG'].astype(float) * 627.5095
        df['HILL'] = df['HILL'].astype(float) * 627.5095

        return df

    def getpiedadf(self, pieda):
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


    def gettgtdf(self, df):
        print('--- ifie near tgt ', self.dist, 'angstrom ----')
        tgtdf = df[((df['I'] == self.tgtfrag) | (df['J'] == self.tgtfrag))]
        # tgtdf = tgtdf.append(df[df['J'] == self.tgtfrag])
        tgtdf_filter = tgtdf[tgtdf['DIST'] < self.dist]

        # print('tgtdf', tgtdf)
        # print('tgtdf, distfilt', tgtdf_filter)
        # print(tgtdf_filter['J'])
        return tgtdf, tgtdf_filter

    def gettgtdf_ffmulti(self, df, frag1, frag2):
        # print('--- ifie frag ', frag1, frag2, '----')
        frag1 = int(frag1)
        frag2 = int(frag2)
        tgtdf_filter = df[((df['I'] == frag1) & (df['J'] == frag2)) | ((df['I'] == frag2) & (df['J'] == frag1))]
        # tgtdf_filter = df[((df['I'] == frag1) & (df['J'] == frag2)) | ((df['I'] == frag2) & (df['J'] == frag1))]
        return tgtdf_filter


    def getifiesummol(self, df):
        molfrags = self.molfrags
        tgtdf_filters = pd.DataFrame(columns=self.icolumn)
        print(self.dist)
        for i in molfrags:
            tgtdf = df[df['I'] == i]
            tgtdf = tgtdf.append(df[df['J'] == i])
            tgtdf_filter = tgtdf[tgtdf['DIST'] < self.dist]
            print(tgtdf_filter)
            tgtdf_filters = tgtdf_filters.append(tgtdf_filter)

        print(tgtdf_filters)

        neighbor_i = [index for index, row in tgtdf_filters.groupby("I")]
        neighbor_j = [index for index, row in tgtdf_filters.groupby("J")]
        neighbors= list(set(neighbor_i + neighbor_j))
        # print(neighbors)
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

        # print('-- ifie permol --')
        ifie_permols = []
        for contactmolfrag in contactmolfrags:
            ifie_permol = pd.DataFrame(columns=self.icolumn)
            for contact in contactmolfrag:
                for tgtfrag in molfrags:
                    # print(contact, tgtfrag)
                    ifie_permol = ifie_permol.append(df[((df['I'] == contact) & (df['J'] == tgtfrag)) | ((df['I'] == tgtfrag) & (df['J'] == contact))])
            ifie_permols.append(ifie_permol)

        count = 0
        ifiesums = [['contactmolfrag', 'tgtmolfrags', 'HF-IFIE', 'MP2-IFIE']]
        # print('contactmolfrag', 'tgtmolfrags', 'HF-IFIE', 'MP2-IFIE')
        for datadf in ifie_permols:
            ifiesums.append([contactmolfrags[count], molfrags, datadf['HF-IFIE'].sum(), datadf['MP2-IFIE'].sum()])
            count += 1

        return contactmolfrags, ifie_permols, ifiesums


    def getpitgtdf(self, pidf, tgtdf_filter):
        # print('--- pieda near tgt ', self.dist, 'angstrom ----')
        print('--- pieda for tgt frag ----')
        pitgtdf = pidf[((pidf['I'] == self.tgtfrag) |(pidf['J'] == self.tgtfrag))]
        # pitgtdf = pitgtdf.append(pidf[pidf['J'] == self.tgtfrag])

        # --- filter ver.1 dist ----
        # pitgtdf_filter =  pd.DataFrame(columns=self.pcolumn)
        # for i in tgtdf_filter['J']:
        #     # print(pitgtdf[pitgtdf['J'] == i])
        #     aa = pitgtdf[pitgtdf['J'] == i]
        #     pitgtdf_filter = pitgtdf_filter.append(aa)

        # --- filter ver.2 except dimer-es resion ----
        # pitgtdf_filter = pitgtdf[~(pitgtdf['CT-mix'] == 0.0)]
        # pitgtdf_filter = pitgtdf_filter[~(pitgtdf_filter['DI(MP2)'] == 0.0)]

        # print(pitgtdf)
        # print(pitgtdf_filter)

        return pitgtdf  #, pitgtdf_filter


    def setupreadparm(self, item1=None, item2=None):
        print('anlmode:' ,self.anlmode)
        print('fragmode:', self.fragmode)
        print('np =:', self.pynp)

        tgtlogs = []
        tgtpdbs = []
        tgttimes = []

        if self.anlmode == 'ff-multi' and self.tgt2type == 'molname':
            self.rpdbflag = True

        # print(self.ilog_head)
        if self.anlmode == 'ff-multi':
            # self.readpdb = True
            print('tgt2type:' ,self.tgt2type)

            if item1 != None:
                self.tgt1frag = item1
                self.tgtfrag = item1

            for i in range(self.start, self.end+1, self.interval):
                tgttimes.append(str(i))
                tgtlogs.append(self.ilog_head + str(i) + self.ilog_tail)
                if self.rpdbflag == True:
                    tgtpdbs.append(self.pdb_head + str(i) + self.pdb_tail)
            print('tgtlogs', tgtlogs)

            if self.tgt2type == 'frag':
                if item2 != None:
                    self.tgt2frag = item2
                print('tgt1frag, tgt2frag', self.tgt1frag, self.tgt2frag)
            if self.tgt2type == 'molname':
                if item2 != None:
                    self.tgt2molname = item2
                print('tgt1frag, tgt2mol', self.tgt1frag, self.tgt2molname)

            self.tgtlogs = tgtlogs
            self.tgttimes = tgttimes
            self.tgtpdbs = tgtpdbs

        else:
            if item1 != None:
                self.tgtlogs = item1
            if item2 != None:
                if self.anlmode == 'mol' and self.selecttype == 'molid':
                    self.tgtmolid = int(item2)
                else:
                    self.tgtfrag = item2


            # --- print setting ---
        # print('logname:', self.tgtlogs)
        # print('dist:', self.dist)
        # print('tgtfragid:', self.tgtfrag)
        # print('anlmode:' , self.anlmode)
        # print('fragmode:', self.fragmode)
        # print('tgt2type:', self.tgt2type)

        ## read pdb
        if self.anlmode == 'ff-multi' and self.tgt2type == 'molname':
            nfs = []
            molnames_perrec = []
            self.assignmolname = False
            for i in range(len(tgtpdbs)):
                self.readpdb(tgtpdbs[i])
                # self.readpdb('aaaa.pdb')
                nf = self.getlognf(tgtlogs[i], self.fragmode)
                nfs.append(nf)
                molnames_perrec.append(self.resnames)

            self.nfs = nfs
            self.molnames_perrec  = molnames_perrec

        elif  self.anlmode == 'fraginmol' or self.rpdbflag == True:
            self.assignmolname = False
            self.readpdb(self.pdbname)
            print(self.resnames)

        # PIEDA
        if self.abinit_ver == 'rev16' or self.abinit_ver == 'rev17':
            self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'Solv(ES)', 'DI(MP2)', 'q(I=>J)']

        ### read fraginfo section
        frags = []
        if self.fragmode != 'manual':
            frags = self.read_fraginfo(self.tgtlogs)
            # print('frags', frags)

        if self.fragmode == 'hybrid':
            getf = frags.pop(hyfrag-1)
            for i in range(self.hynum):
                frags.append(getf)
            # print('frags', frags)

        self.frags = frags

    def read_ifiepieda(self, fname):
        ifie = []
        count = 0
        pieda = []
        pcount = 0
        pflag = False

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
            if itemList[1] == 'PIEDA':
                flag = False
                pflag = True
                continue
            if flag is True:
                count += 1
            if flag is True and count > 2:
                ifie.append(itemList)
            # print itemList
            if len(itemList) < 2:
                continue

            if itemList[1] == 'PIEDA':
                pflag = True
                # head.append(itemList)
                continue
            if itemList[1] == 'Mulliken':
                # flag = False
                break
            if pflag is True:
                pcount += 1
            if pflag is True and pcount > 2:
                pieda.append(itemList)

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

        return ifie, pieda
        # print ifie



#     def readifies(self, tgtlog):
#         # tgtlog, tgttime = args
#         print('read', tgtlog)
#         ifie = self.read_ifie(tgtlog)
#         ifdfs = self.getifiedf(ifie)
#         return ifdfs

    def read_ifiepiedas(self, tgtlog):
        # tgtlog, tgttime = args
        print('read', tgtlog)
        ifie, pieda = self.read_ifiepieda(tgtlog)
        # print('ifie', ifie[0], 'pieda', pieda[0])
        ifdfs = self.getifiedf(ifie)
        pidfs = self.getpiedadf(pieda)
        # print('ifdfs', ifdfs)
        # print('pidfs', pidfs)

        return [ifdfs, pidfs]

    def read_ifpif90(self, tgtlog):
        print('read', tgtlog)
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
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            np.ctypeslib.ndpointer(dtype=np.int32),
            ]

        instr = tgtlog
        enc_str = instr.encode('utf-8')
        instr = create_string_buffer(enc_str)

        ifi = np.empty(10000000, dtype=np.int32)
        ifj = np.empty(10000000, dtype=np.int32)
        pii = np.empty(10000000, dtype=np.int32)
        pij = np.empty(10000000, dtype=np.int32)
        dist = np.empty(10000000, dtype=np.float64)
        hfifie = np.empty(10000000, dtype=np.float64)
        mp2ifie = np.empty(10000000, dtype=np.float64)
        prtype1 = np.empty(10000000, dtype=np.float64)
        grimme = np.empty(10000000, dtype=np.float64)
        jung = np.empty(10000000, dtype=np.float64)
        hill = np.empty(10000000, dtype=np.float64)
        es = np.empty(10000000, dtype=np.float64)
        ex = np.empty(10000000, dtype=np.float64)
        ct = np.empty(10000000, dtype=np.float64)
        di = np.empty(10000000, dtype=np.float64)
        qval = np.empty(10000000, dtype=np.float64)
        fdimesint = np.empty(10000000, dtype=np.int32)
        ifpair = np.empty(1, dtype=np.int32)
        pipair = np.empty(1, dtype=np.int32)

        f.readifiepieda_(instr, ifi, ifj, pii, pij, dist, hfifie, mp2ifie, prtype1, grimme, jung, hill, es, ex, ct, di, qval, fdimesint, ifpair, pipair)

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

        pidf = pd.DataFrame(columns=self.pcolumn)
        pidf['I'] = copy.deepcopy(pii[:pipair[0]])
        pidf['J'] = copy.deepcopy(pij[:pipair[0]])
        pidf['ES'] = copy.deepcopy(es[:pipair[0]])
        pidf['EX'] = copy.deepcopy(ex[:pipair[0]])
        pidf['CT-mix'] = copy.deepcopy(ct[:pipair[0]])
        pidf['DI(MP2)'] = copy.deepcopy(di[:pipair[0]])
        pidf['q(I=>J)'] = copy.deepcopy(qval[:pipair[0]])

        return [ifdf, pidf]

    def readifiewrap(self, item1=None, item2=None):
        self.setupreadparm(item1, item2)
        if self.anlmode == 'ff-multi':
            print('## read multi mode')
            ifdfs = []
            skips = []
            args = []
            st = time.time()
            p = Pool(self.pynp)
            if self.f90soflag == True:
                print("use fortran library")
                ifpidfs = p.map(self.read_ifpif90, self.tgtlogs)
            else:
                ifpidfs = p.map(self.read_ifiepiedas, self.tgtlogs)
            p.close()
            print('read elapsed', time.time() - st)

            for i in range (len(ifpidfs)):
                # print('len tgt', i, ':', len(ifpidfs[i][0]))
                if len(ifpidfs[i][0]) == 0:
                    print('Warning: time', self.tgttimes[i], 'not converged: skip data')
                    skips.append(i)

            dellist = lambda items, indexes: [item for index, item in enumerate(items) if index not in indexes]
            self.tgttimes = dellist(self.tgttimes,skips)
            # print('tgttimes', self.tgttimes)

            if self.tgt2type == 'molname':
                print('### read frag info ###')
                molfragss_perrec = []
                for i in range(len(self.tgtpdbs)):
                    molfragss = self.getallmolfrags(self.tgtlogs[i], ifpidfs[i][0], self.nfs[i])
                    # print('frags_permol\n', molfragss)
                    molfragss_perrec.append(molfragss)
                # print(molfragss_perrec)
                # print(len(molnames), len(molfragss))
                self.molfragss_perrec = molfragss_perrec

            st = time.time()
            for i in range(len(ifpidfs)):
                if len(ifpidfs[i][0]) != 0:
                    self.ifdfs.append(ifpidfs[i][0])
                    self.pidfs.append(ifpidfs[i][1])
            print('elapsed', time.time() - st)

        else:
            print('## read single mode')
            ifie, pieda = self.read_ifiepieda(self.tgtlogs)
            df = self.getifiedf(ifie)

            self.ifdf = df

            # pieda = self.read_pieda(self.tgtlogs)
            pidf = self.getpiedadf(pieda)

            self.pidf = pidf

            # self.ifdfs = []
            # self.ifdfs.append(df)

        ### read pieda (from log) section
        return self

            # print(df.head())

    ### filter section
    def filterifiewrap(self, dist=None):

        tgt2type = self.tgt2type
        if dist != None:
            self.dist = dist
        # frag mode
        if self.anlmode == 'frag':
            # print(self.ifdf)
            tgtdf, ifdf_filter = self.gettgtdf(self.ifdf)
            # print(ifdf_filter)
            self.ifdf_filter = ifdf_filter

        # mol-mol mode
        if self.anlmode == 'mol':
            df = self.ifdf
            if self.selecttype == 'fragid':
                molfrags = self.getmolfrags(self.tgt1frag, df)
                print('target-frags:', molfrags)
            elif self.selecttype == 'molid':
                nf = self.getlognf(self.tgtlogs, self.fragmode)
                molfragss = self.getallmolfrags(self.tgtlogs, df, nf)
                print('frags_permol\n', molfragss)
                molfrags = molfragss[self.tgtmolid-1]
                self.molfrags = molfrags

                contactmolfrags, ifie_permols, ifiesums = self.getifiesummol(df)

                self.contactmolfrags = contactmolfrags
                self.ifie_permols = ifie_permols
                self.ifiesums =ifiesums

        # fraginmol mode
        if self.anlmode == 'fraginmol':
            df = self.ifdf
            tgt1_lofrag = self.tgt1_lofrag
            tgt2_lofrag = self.tgt2_lofrag
            tgt2molname = self.tgt2molname
            nf = self.getlognf(self.tgtlogs, self.fragmode)
            print(nf)
            molfragss = self.getallmolfrags(self.tgtlogs, df, nf)
            print(molfragss)
            tgtmol = self.tgtmolid - 1

            tgt2_glofrags = []
            tgt1_glofrag = molfragss[tgtmol][tgt1_lofrag - 1]
            print('centermolfrag:', tgt1_glofrag)
            print('tgt2molname', tgt2molname)
            for i in range(len(self.resnames)):
                if self.resnames[i] == tgt2molname:
                    tgt2frag = molfragss[i][tgt2_lofrag - 1]
                    tgt2_glofrags.append(tgt2frag)
            print('tgt2_glofrags', tgt2_glofrags)
            tgtdf = df[df['I'] == tgt1_glofrag]
            tgtdf = tgtdf.append(df[df['J'] == tgt1_glofrag])
            tgtdf = tgtdf[tgtdf['DIST'] < self.dist]

            tgt_new2 = pd.DataFrame()
            for tgt2_glofrag in tgt2_glofrags:
                tgt_new = tgtdf[(tgtdf['I'] == tgt2_glofrag) |(tgtdf['J'] == tgt2_glofrag)]
                tgt_new2 = tgt_new2.append(tgt_new)

            print('tgt_new2\n', tgt_new2)

            self.tgt1_glofrag = tgt1_glofrag
            self.tgt2_glofrags = tgt2_glofrags
            self.ifdf_filters = tgt_new2

        # multi mode
        if self.anlmode == 'ff-multi':

            ifdfs = self.ifdfs
            tgtdf_filters = pd.DataFrame()
            print('tgttimes', self.tgttimes)

            count = 0
            if tgt2type == 'frag':
                for df in ifdfs:
                    # print('read', tgtlogs[count])
                    # print(df)
                    tgtdf_filters = tgtdf_filters.append(self.gettgtdf_ffmulti(df, self.tgt1frag, self.tgt2frag))
                    count += 1
                tgtdf_filters['TIME'] = self.tgttimes
                # print(tgtdf_filters)
                self.ifdf_filters = tgtdf_filters

            if tgt2type == 'molname':
                # get tgt frag id
                tgtmolfrags_perrec = []
                HF_IFIE_sums = []
                MP2_IFIE_sums = []
                PR_TYPE1_sums = []
                GRIMME_sums = []
                JUNG_sums = []
                HILL_sums = []

                molnames_perrec = self.molnames_perrec
                # print(molnames_perrec)
                for i in range(len(molnames_perrec)):
                    tgtmolfrags = []
                    for j in range(len(molnames_perrec[i])):
                        try:
                            if molnames_perrec[i][j] == self.tgt2molname:
                                tgtmolfrags.append(self.molfragss_perrec[i][j])
                        except:
                            continue
                    # print(tgtmolfrags)

                    tgtmolfrags_perrec.append(tgtmolfrags)
                    self.tgtmolfrags_perrec =  tgtmolfrags_perrec
                # print(len(tgtmolfrags_perrec))

                frag1 = self.tgt1frag
                if type(frag1) == int:
                    frag1s = [frag1]
                    print('frag1s', frag1s)
                else:
                    frag1s = copy.deepcopy(frag1)

                self.frag1s =frag1s
                for i in range(len(tgtmolfrags_perrec)):
                    df = ifdfs[i]
                    if df.empty:
                        # print(i, 'empty!!')
                        continue
                    tgtdf_filters = pd.DataFrame()
                    for frag1p in frag1s:
                        print('frag1p', frag1p)
                        for tgtmolfrags in tgtmolfrags_perrec[i]:
                            for tgt2frag in tgtmolfrags:
                                tgtdf_filters = tgtdf_filters.append(self.gettgtdf_ffmulti(df, frag1p, tgt2frag))
                                count += 1
                        print(tgtdf_filters.tail())

                    HF_IFIE_sums.append(tgtdf_filters['HF-IFIE'].sum())
                    MP2_IFIE_sums.append(tgtdf_filters['MP2-IFIE'].sum())
                    PR_TYPE1_sums.append(tgtdf_filters['PR-TYPE1'].sum())
                    GRIMME_sums.append(tgtdf_filters['GRIMME'].sum())
                    JUNG_sums.append(tgtdf_filters['JUNG'].sum())
                    HILL_sums.append(tgtdf_filters['HILL'].sum())

                tgtifdfsum = pd.DataFrame()
                tgtifdfsum['HF-IFIE'] = HF_IFIE_sums
                tgtifdfsum['MP2-IFIE'] = MP2_IFIE_sums
                tgtifdfsum['PR-TYPE1'] = PR_TYPE1_sums
                tgtifdfsum['GRIMME'] = GRIMME_sums
                tgtifdfsum['JUNG'] = JUNG_sums
                tgtifdfsum['HILL'] = HILL_sums
                tgtifdfsum['TIME'] = self.tgttimes
                self.tgtmolfrags_perrec = tgtmolfrags_perrec

                print(tgtifdfsum)
                self.ifdfsum = tgtifdfsum

            if tgt2type == 'dist':
                print('filter-dist:', self.dist)
                # get tgt frag id
                HF_IFIE_sums = []
                MP2_IFIE_sums = []
                PR_TYPE1_sums = []
                GRIMME_sums = []
                JUNG_sums = []
                HILL_sums = []

                tgtdf_filters = pd.DataFrame()
                for ifdf in self.ifdfs:
                    tgtdf, tgtdf_filter = self.gettgtdf(ifdf)
                    # print(tgtdf_filter.head())
                    tgtdf_filters = tgtdf_filters.append(tgtdf_filter)

                    HF_IFIE_sums.append(tgtdf_filter['HF-IFIE'].sum())
                    MP2_IFIE_sums.append(tgtdf_filter['MP2-IFIE'].sum())
                    PR_TYPE1_sums.append(tgtdf_filter['PR-TYPE1'].sum())
                    GRIMME_sums.append(tgtdf_filter['GRIMME'].sum())
                    JUNG_sums.append(tgtdf_filter['JUNG'].sum())
                    HILL_sums.append(tgtdf_filter['HILL'].sum())

                tgtifdfsum = pd.DataFrame()
                tgtifdfsum['HF-IFIE'] = HF_IFIE_sums
                tgtifdfsum['MP2-IFIE'] = MP2_IFIE_sums
                tgtifdfsum['PR-TYPE1'] = PR_TYPE1_sums
                tgtifdfsum['GRIMME'] = GRIMME_sums
                tgtifdfsum['JUNG'] = JUNG_sums
                tgtifdfsum['HILL'] = HILL_sums
                tgtifdfsum['TIME'] = self.tgttimes

                print(tgtifdfsum)
                self.ifdf_filters = tgtdf_filters
                self.ifdfsum = tgtifdfsum

        return self


    def readpiedawrap(self, item1=None, item2=None):
        print('--- pieda ---')

        # self.setupreadparm(item1, item2)
        if self.abinit_ver == 'rev16' or self.abinit_ver == 'rev17':
            self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'Solv(ES)', 'DI(MP2)', 'q(I=>J)']
        else:
            self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']

        ### read fraginfo section
        frags = []
        if self.fragmode != 'manual':
            frags = self.read_fraginfo(self.tgtlogs)
            # print('frags', frags)

        if self.fragmode == 'hybrid':
            getf = frags.pop(hyfrag-1)
            for i in range(self.hynum):
                frags.append(getf)
            # print('frags', frags)

        self.frags = frags
        ### read pieda (from log) section
        if self.anlmode == 'ff-multi':
            pidfs = []
            for logname in self.tgtlogs:
                pieda = self.read_pieda(logname)
                if len(pieda) != 0:
                    pidfs.append(self.getpiedadf(pieda))

            self.pidfs = pidfs
        else:
            pieda = self.read_pieda(self.tgtlogs)
            pidf = self.getpiedadf(pieda)

            self.pidf = pidf
        return self

    def filterpiedawrap(self):
        ### filter
        # frag-frag mode
        frags = self.frags
        if self.anlmode == 'frag':
            pitgtdf = self.getpitgtdf(self.pidf, self.ifdf_filter)
            if self.fragmode != 'manual':
                # print('len_frags', len(frags))
                #assign resname(e.g. Gly6)
                for i in range(1, len(frags) + 1):
                    pitgtdf.I = pitgtdf.I.replace(i, frags[i-1])
                    pitgtdf.J = pitgtdf.J.replace(i, frags[i-1])

            self.pitgtdf = pitgtdf
        # ff-multi mode
        if self.anlmode == 'ff-multi':
            pidfs = self.pidfs
            if self.tgt2type == 'frag':
                pitgtdf_filters = pd.DataFrame(columns=self.pcolumn)
                # print('first', pitgtdf_filters)

                count = 0
                for pidf in pidfs:
                    # print('read', tgtlogs[count])
                    pitgtdf_filter = self.gettgtdf_ffmulti(pidf, self.tgt1frag, self.tgt2frag)
                    if len(pitgtdf_filter) == 0:
                        esdata = self.ifdf_filters.iat[count, 4]
                        pitgtdf_filter.loc[str(count)] = [self.tgt1frag, self.tgt2frag, esdata, "", "", "", ""]
                    pitgtdf_filters = pitgtdf_filters.append(pitgtdf_filter)
                    count += 1
                pitgtdf_filters['TIME'] = self.tgttimes
                self.pidf_filters = pitgtdf_filters

            if self.tgt2type == 'molname':
                ES_sums = []
                EX_sums = []
                CT_sums = []
                DI_sums = []
                q_sums = []
                for i in range(len(self.tgtmolfrags_perrec)):
                    pidf = pidfs[i]
                    if pidf.empty:
                        continue
                    pitgtdf_filters = pd.DataFrame()
                    for frag1p in self.frag1s:
                        for tgtmolfrags in self.tgtmolfrags_perrec[i]:
                            for tgt2frag in tgtmolfrags:
                                pitgtdf_filters = pitgtdf_filters.append(self.gettgtdf_ffmulti(pidf, frag1p, tgt2frag))
                    ES_sums.append(pitgtdf_filters['ES'].sum())
                    EX_sums.append(pitgtdf_filters['EX'].sum())
                    CT_sums.append(pitgtdf_filters['CT-mix'].sum())
                    DI_sums.append(pitgtdf_filters['DI(MP2)'].sum())
                    q_sums.append(pitgtdf_filters['q(I=>J)'].sum())

                pitgtdfsum = pd.DataFrame()
                pitgtdfsum['ES'] = ES_sums
                pitgtdfsum['EX'] = EX_sums
                pitgtdfsum['CT-mix'] = CT_sums
                pitgtdfsum['DI(MP2)'] = DI_sums
                pitgtdfsum['q(I=>J)'] = q_sums
                pitgtdfsum['TIME'] = self.tgttimes

                print(pitgtdfsum)
                self.pidfsum = pitgtdfsum
            if  self.tgt2type == 'dist':
                ES_sums = []
                EX_sums = []
                CT_sums = []
                DI_sums = []
                q_sums = []

                pitgtdfs = pd.DataFrame()
                for i in range(len(self.pidfs)):
                    # print(self.pidfs[i].head())
                    # print('headheadhead', self.ifdf_filters.head())
                    pitgtdf = self.getpitgtdf(self.pidfs[i], self.ifdf_filters)
                    if self.fragmode != 'manual':
                        # print('len_frags', len(frags))
                        #assign resname(e.g. Gly6)
                        for i in range(1, len(frags) + 1):
                            pitgtdf.I = pitgtdf.I.replace(i, frags[i-1])
                            pitgtdf.J = pitgtdf.J.replace(i, frags[i-1])

                    ES_sums.append(pitgtdf['ES'].sum())
                    EX_sums.append(pitgtdf['EX'].sum())
                    CT_sums.append(pitgtdf['CT-mix'].sum())
                    DI_sums.append(pitgtdf['DI(MP2)'].sum())
                    q_sums.append(pitgtdf['q(I=>J)'].sum())

                    pitgtdfs = pitgtdfs.append(pitgtdf)

                pitgtdfsum = pd.DataFrame()
                pitgtdfsum['ES'] = ES_sums
                pitgtdfsum['EX'] = EX_sums
                pitgtdfsum['CT-mix'] = CT_sums
                pitgtdfsum['DI(MP2)'] = DI_sums
                pitgtdfsum['q(I=>J)'] = q_sums
                pitgtdfsum['TIME'] = self.tgttimes

                print(pitgtdfsum)
                self.pidfsum = pitgtdfsum

        # mol-mol mode
        if self.anlmode == 'mol':
            pidf = self.pidf
            pieda_permols = []
            for contactmolfrag in self.contactmolfrags:
                pieda_permol = pd.DataFrame(columns=self.pcolumn)
                for contact in contactmolfrag:
                    for tgtfrag in self.molfrags:
                        # print(contact, tgtfrag)
                        pieda_permol = pieda_permol.append(pidf[((pidf['I'] == contact) & (pidf['J'] == tgtfrag)) | ((pidf['I'] == tgtfrag) & (pidf['J'] == contact))])
                pieda_permols.append(pieda_permol)
                # print(pieda_permol, file=plogdt)

            count = 0
            if self.abinit_ver == 'rev17' or self.abinit_ver == 'rev16':
                piedasums = [['contactmolfrag', 'tgtmolfrags', 'ES', 'EX', 'Solv(ES)', 'CT-mix', 'DI(MP2)', 'a(I=>J)']]
            else:
                piedasums = [['contactmolfrag', 'tgtmolfrags', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'a(I=>J)']]
            for datadf in pieda_permols:
                piedasums.append([self.contactmolfrags[count], self.molfrags, datadf['ES'].sum(), datadf['EX'].sum(), datadf['CT-mix'].sum(), datadf['DI(MP2)'].sum(),  datadf['q(I=>J)'].sum()])
                count += 1

                self.piedasums = piedasums
                self.pieda_permols = pieda_permols

        # fraginmol mode
        if self.anlmode == 'fraginmol':
            pidf = self.pidf
            tgt1_glofrag = self.tgt1_glofrag
            tgt2_glofrags = self.tgt2_glofrags
            print('--- pieda ----')
            pitgtdf = pidf[pidf['I'] == tgt1_glofrag]
            pitgtdf = pitgtdf.append(pidf[pidf['J'] == tgt1_glofrag])

            pitgt_new2 = pd.DataFrame()
            for tgt2_glofrag in tgt2_glofrags:
                pitgt_new = pitgtdf[(pitgtdf['I'] == tgt2_glofrag) |(pitgtdf['J'] == tgt2_glofrag)]
                pitgt_new2 = pitgt_new2.append(pitgt_new)

            if self.fragmode != 'manual':
                print('len_frags', len(frags))

                for i in range(1, len(frags) + 1):
                    pitgt_new2.I = pitgt_new2.I.replace(i, frags[i-1])
                    pitgt_new2.J = pitgt_new2.J.replace(i, frags[i-1])
            self.pidf_filters = pitgt_new2

        return self


    def writecsvwrap(self, head=None):
        # -------------
        # --- write ---
        # -------------
        path = 'csv'
        tgt2type = self.tgt2type
        if os.path.exists('csv') == False:
            os.mkdir('csv')


        if self.anlmode == 'frag':
            if head == None:
                head, ext = os.path.splitext(self.tgtlogs)

            tgtid = self.tgtfrag
            try:
                ohead = head + '-' + str(tgtid) + '-' + frags[tgtid - 1]
            except:
                ohead = head + '-' + str(tgtid)

            # print(self.pitgtdf.head())

            # tgtdf.to_csv(path + '/' + ohead + '-ifie.csv')
            self.ifdf_filter.to_csv(path + '/' + ohead + '-ifie_' + 'dist' + str(self.dist) + '.csv')
            self.pitgtdf.to_csv(path + '/' + ohead + '-pieda.csv')

            print(ohead + '-ifie.csv', ohead + '-pieda.csv generated.')


        if self.anlmode == 'ff-multi':
            tgt1frag = self.tgt1frag
            if tgt2type == 'frag':
                tgt2frag = self.tgt2frag
                oifie = 'frag' + str(tgt1frag) + '-frag' + str(tgt2frag) + '-ifie.csv'
                opieda = 'frag' + str(tgt1frag) + '-frag' + str(tgt2frag) + '-pieda.csv'
                self.ifdf_filters.to_csv(path + '/' + oifie)
                self.pidf_filters.to_csv(path + '/' + opieda)
                print(path + '/' + oifie, path + '/' + opieda, 'was created.')

            if tgt2type == 'molname':
                tgt2molname = self.tgt2molname
                oifie = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-ifiesum.csv'
                opieda = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-piedasum.csv'
                self.ifdfsum.to_csv(path + '/' + oifie)
                self.pidfsum.to_csv(path + '/' + opieda)
                print(path + '/' + oifie, path + '/' + opieda, 'was created.')


            if tgt2type == 'dist':
                tgt2dist = self.dist
                oifie = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-ifiesum.csv'
                opieda = 'frag' + str(tgt1frag) + '-piedasum.csv'
                self.ifdfsum.to_csv(path + '/' + oifie)
                self.pidfsum.to_csv(path + '/' + opieda)
                print(path + '/' + oifie, path + '/' + opieda, 'was created.')


        if self.anlmode == 'mol':
            if head == None:
                head, ext = os.path.splitext(self.tgtlogs)

            dist = self.dist
            selecttype = self.selecttype
            ifiesums= self.ifiesums
            piedasums = self.piedasums
            ifie_permols = self.ifie_permols
            pieda_permols = self.pieda_permols
            if selecttype == 'molid':
                tgtid = self.tgtmolid
            else:
                tgtid = self.tgtfrag

            ilogdtname = path + '/' + head + '_ifie-mol-' +  selecttype + str(tgtid) + 'dist' + str(dist) + '.txt'
            isumname = path + '/' + head + '_ifiesum-mol-' + selecttype + str(tgtid) + 'dist' + str(dist) + '.csv'
            plogdtname = path + '/' + head + '_pieda-mol-' + selecttype + str(tgtid) + 'dist' + str(dist) + '.txt'
            psumname = path + '/' + head + '_piedasum-mol-' + selecttype  + str(tgtid) + 'dist' + str(dist) + '.csv'

            # write section
            ilogdt = open(ilogdtname, 'w')
            for ifie_permol in ifie_permols:
                print(ifie_permol, file=ilogdt)

            plogdt = open(plogdtname, 'w')
            for pieda_permol in pieda_permols:
                print(pieda_permol, file=plogdt)

            with open(isumname, 'w') as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerows(ifiesums)

            with open(psumname, 'w') as f:
                writer = csv.writer(f, lineterminator='\n')
                # writer.writerow(list)
                writer.writerows(piedasums)

            print('---out---')
            print(ilogdtname)
            print(isumname)
            print(plogdtname)
            print(psumname)

        if self.anlmode == 'fraginmol':
            if head == None:
                head, ext = os.path.splitext(self.tgtlogs)

            ohead = head + '-' 'frag' + str(self.tgt1_glofrag) + '-mol' + str(self.tgt2molname) + 'frag' + str(self.tgt2_lofrag)
            # print(lpitgt_new2.head())

            self.ifdf.to_csv(path + '/' + head + '-ifie.csv')
            self.ifdf_filters.to_csv(path + '/' + ohead + '-ifie_'  + 'dist' + str(self.dist) + '.csv')
            self.pidf_filters.to_csv(path + '/' + ohead + '-pieda.csv')
            print(ohead + '-ifie.csv', ohead + '-pieda.csv generated.')



