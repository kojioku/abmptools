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
        self.ifdfsumcolumn = ['HF-IFIE', 'MP2-IFIE', 'PR-TYPE1', 'GRIMME', 'JUNG', 'HILL']
        self.pidfsumcolumn = ['ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']
        self.logMethod = 'MP2'


        self.anlmode= 'frag' #frag, 'mol', 'fraginmol', 'multi'
        self.fragmode = 'auto'  #'hybrid', 'auto', 'manual'
        self.dist = 10.0
        self.tgt1frag = None

        self.rpdbflag = False
        self.pdbname = None   # 'iss2-spg2-ok20200130opt-for_abmp.pdb'

        # -- for mol mode or multi mode--
        self.tgt2type = 'frag' #frag: mol-frag, mol: mol-mol

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
        self.hyfrag = None #320
        self.hynum = None

        self.pynp = 3

        mydir = os.path.dirname(os.path.abspath(__file__))
        self.f90sofile = mydir + '/f90/bin/readifiepiedalib.so'
        self.f90soflag = True

        # for svd
        self.matrixtype = 'normal'
        self.exceptfrag = []
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


    def gettgtpidf_n2ffmatrix(self):
        print('\n--- generate pieda', str(self.tgt1frag), str(self.tgt2frag), 'ffmatrix ---\n')
        esdf = pd.DataFrame(index=self.tgt2frag)
        exdf = pd.DataFrame(index=self.tgt2frag)
        ctdf = pd.DataFrame(index=self.tgt2frag)
        count = 0

        df = self.pidf
        for f1 in self.tgt1frag:
            fragids = []
            tgtdf = df[(df['I'] == f1) | (df['J'] == f1)]
            tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

            fragis = tgtdf_filter['I'].values.tolist()
            fragjs = tgtdf_filter['J'].values.tolist()
            for i in range(len(fragis)):
                if fragis[i] != f1:
                    fragids.append(fragis[i])
                else:
                    fragids.append(fragjs[i])
            # print(fragids)
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

            for i in range(len(self.tgt2frag)):
                tgtid = self.tgt2frag[i]
                if tgtid  == f1:
                    print('error!! target frag1 and target 2 is duplicate!!')
                    sys.exit()
                if tgtid in fragids:
                    es.append(esbuf[fragids.index(tgtid)])
                    ex.append(exbuf[fragids.index(tgtid)])
                    ct.append(ctbuf[fragids.index(tgtid)])
                else:
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


    def gettgtdf_n2ffmatrix(self):
        # gettgtdf_normal to times-frags
        print('\n--- generate ifie', str(self.tgt1frag), str(self.tgt2frag), 'ffmatrix---\n')
        hfdf = pd.DataFrame(index=self.tgt2frag)
        mp2corrdf = pd.DataFrame(index=self.tgt2frag)
        prmp2corrdf = pd.DataFrame(index=self.tgt2frag)
        mp2tdf =  pd.DataFrame(index=self.tgt2frag)
        prmp2tdf = pd.DataFrame(index=self.tgt2frag)

        mp3corrdf = pd.DataFrame(index=self.tgt2frag)
        usermp3corrdf = pd.DataFrame(index=self.tgt2frag)
        mp3tdf =  pd.DataFrame(index=self.tgt2frag)
        usermp3tdf = pd.DataFrame(index=self.tgt2frag)

        df = self.ifdf
        count = 0

        if self.logMethod == 'MP2':
            for f1 in self.tgt1frag:
                fragids = []
                tgtdf = df[(df['I'] == f1) | (df['J'] == f1)]
                tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

                # fragis = tgtdf_filter['I'].values.tolist()
                # fragjs = tgtdf_filter['J'].values.tolist()
                # for i in range(len(fragis)):
                #     if fragis[i] != self.tgt1frag:
                #         fragids.append(fragis[i])
                #     else:
                #         fragids.append(fragjs[i])
                # print(fragids)
                hfifie = 0
                mp2corr = 0
                prmp2corr = 0
                hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
                mp2corr = tgtdf_filter['MP2-IFIE'].values.tolist()
                prmp2corr = tgtdf_filter['PR-TYPE1'].values.tolist()

                mp2total = []
                prmp2total  = []
                for i in range(len(hfifie)):
                    mp2total.append(hfifie[i] + mp2corr[i])
                    prmp2total.append(hfifie[i] + prmp2corr[i])

                # print('hfifie', hfifie)
                # print('tgtfrag', self.tgt2frag)

                hfdf[str(f1)] = hfifie
                mp2corrdf[str(f1)] = mp2corr
                prmp2corrdf[str(f1)] = prmp2corr
                mp2tdf[str(f1)] = mp2total
                prmp2tdf[str(f1)] = prmp2total

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


        elif self.logMethod == 'MP3':
            for f1 in self.tgt1frag:
                fragids = []
                tgtdf = df[(df['I'] == f1) | (df['J'] == f1)]
                tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

                hfifie = 0
                mp3corr = 0
                usermp3corr = 0
                hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
                mp3corr = tgtdf_filter['MP3-IFIE'].values.tolist()
                usermp3corr = tgtdf_filter['USER-MP3'].values.tolist()

                mp3total = []
                usermp3total  = []
                for i in range(len(hfifie)):
                    mp3total.append(hfifie[i] + mp3corr[i])
                    usermp3total.append(hfifie[i] + usermp3corr[i])

                # print('hfifie', hfifie)
                # print('tgtfrag', self.tgt2frag)

                hfdf[str(f1)] = hfifie
                mp3corrdf[str(f1)] = mp3corr
                usermp3corrdf[str(f1)] = usermp3corr
                mp3tdf[str(f1)] = mp3total
                usermp3tdf[str(f1)] = usermp3total

                count += 1

            print ('HF\n', hfdf.head())
            print ('MP3corr\n', mp3corrdf.head())
            print ('USER-MP3corr\n', usermp3corrdf.head())
            print ('MP3total\n', mp3tdf.head())
            print ('USER-MP3total\n', usermp3tdf.head())

            self.hfdf = hfdf
            self.mp3corrdf = mp3corrdf
            self.usermp3corrdf = usermp3corrdf
            self.mp3tdf = mp3tdf
            self.usermp3tdf = usermp3tdf

        return


    def gettgtdf_n2tfmatrix(self, i, df, tgt1frag):
        # gettgtdf_normal to times-frags
        hfdf = pd.DataFrame(index=self.tgt2frag)
        mp2corrdf = pd.DataFrame(index=self.tgt2frag)
        prmp2corrdf = pd.DataFrame(index=self.tgt2frag)
        mp2tdf =  pd.DataFrame(index=self.tgt2frag)
        prmp2tdf = pd.DataFrame(index=self.tgt2frag)


        fragids = []
        tgtdf = df[(df['I'] == tgt1frag) | (df['J'] == tgt1frag)]
        tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

        # fragis = tgtdf_filter['I'].values.tolist()
        # fragjs = tgtdf_filter['J'].values.tolist()
        # for i in range(len(fragis)):
        #     if fragis[i] != self.tgt1frag:
        #         fragids.append(fragis[i])
        #     else:
        #         fragids.append(fragjs[i])
        # print(fragids)
        hfifie = 0
        mp2corr = 0
        prmp2corr = 0
        hfifie = tgtdf_filter['HF-IFIE'].values.tolist()
        mp2corr = tgtdf_filter['MP2-IFIE'].values.tolist()
        prmp2corr = tgtdf_filter['PR-TYPE1'].values.tolist()

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


        # print (hfdf.head())
        # print (mp2corrdf.head())
        # print (prmp2corrdf.head())
        # print (mp2tdf.head())
        # print (prmp2tdf.head())

        return hfdf, mp2corrdf, prmp2corrdf, mp2tdf, prmp2tdf


    def gettgtpidf_n2tfmatrix(self, i, df, hfdf, tgt1frag):
        esdf = pd.DataFrame(index=self.tgt2frag)
        exdf = pd.DataFrame(index=self.tgt2frag)
        ctdf = pd.DataFrame(index=self.tgt2frag)

        fragids = []
        tgtdf = df[(df['I'] == tgt1frag) | (df['J'] == tgt1frag)]
        tgtdf_filter = tgtdf[(tgtdf['I'].isin(self.tgt2frag)) | (tgtdf['J'].isin(self.tgt2frag))]

        fragis = tgtdf_filter['I'].values.tolist()
        fragjs = tgtdf_filter['J'].values.tolist()
        # print('tgtdf_filter', tgtdf_filter)
        # print('fragis', fragis)
        for j in range(len(fragis)):
            if fragis[j] != tgt1frag:
                fragids.append(fragis[j])
            else:
                fragids.append(fragjs[j])
        # print(fragids)
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

        for j in range(len(self.tgt2frag)):
            tgtid = self.tgt2frag[j]
            if tgtid == tgt1frag:
                continue
            if tgtid in fragids:
                es.append(esbuf[fragids.index(tgtid)])
                ex.append(exbuf[fragids.index(tgtid)])
                ct.append(ctbuf[fragids.index(tgtid)])
            else:
                es.append(hfdf.loc[tgtid, str(self.tgttimes[i])])
                ex.append(0.0)
                ct.append(0.0)

        esdf[self.tgttimes[i]] = es
        exdf[self.tgttimes[i]] = ex
        ctdf[self.tgttimes[i]] = ct

#         print (esdf.head())
#         print (exdf.head())
#         print (ctdf.head())

        # print('esdf\n', esdf)
        return esdf, exdf, ctdf

        # return hfdf, mp2corrdf, prmp2corrdf, mp2tdf, prmp2tdf

    def gettgtdf_fd(self, df):
        print('--- ifie near tgt ', self.dist, 'angstrom ----')
        tgtdf = df[(df['I'] == self.tgt1frag[0]) | (df['J'] == self.tgt1frag[0])]
        # tgtdf = tgtdf.append(df[df['J'] == self.tgt1frag])
        tgtdf_filter = tgtdf[tgtdf['DIST'] < self.dist]

        # print('tgtdf', tgtdf)
        # print('tgtdf, distfilt', tgtdf_filter)
        # print(tgtdf_filter['J'])
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
        tgtdf_filter = df[((df['I'] == frag1) & (df['J'].isin(frag2))) | ((df['I'].isin(frag2)) & (df['J'] == frag1))]
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


    def getfiltifpiff(self, i, ifdf, pidf):
        # in (class var: molnames(rec i), tgt2frag, tgt1frag, tgttimes, matrixtype
        #     local var:i,  ifdf, pidf
        # out tgtdf_filter, pitgtdf_filter

        # ifie
        ifdf_filter = pd.DataFrame()
        ifdf_filter = self.gettgtdf_ff(ifdf, self.tgt1frag[0], self.tgt2frag[0])
        ifdf_filter = ifdf_filter.rename(index={ifdf_filter.index[0]:self.tgttimes[i]})
        # print(ifdf_filter)

        # pieda
        pitgtdf_filter = self.gettgtdf_ff(pidf, self.tgt1frag[0], self.tgt2frag[0])
        if len(pitgtdf_filter) == 0:
            esdata = ifdf_filter['HF-IFIE'][0]
            tmpps = pd.Series([self.tgt1frag[0], self.tgt2frag[0], esdata, 0.0, 0.0, 0.0, 0.0], index=self.pcolumn, name =self.tgttimes[i])
            pitgtdf_filter = pitgtdf_filter.append(tmpps)
        else:
            pitgtdf_filter = pitgtdf_filter.rename(index={pitgtdf_filter.index[0]:self.tgttimes[i]})

            # print(pitgtdf_filter)

        return ifdf_filter, pitgtdf_filter


    def getfiltifpifm(self, i, ifdf, pidf):
        # in (class var: log(rec i), nf(rec i), molnames(rec i), tgt2molname, tgt1frag,
        #     local var: ifdf, pidf
        # out tgtdf_filter, tgtifdfsum
        #     pitgtdf_filter, pitgtdfsum

        print('### read frag info ###')

        molfragss = self.getallmolfrags(self.tgtlogs[i], ifdf, self.nfs[i])
        # get tgt frag id
        # IFIE
        molnames_perrec = self.molnames_perrec
        tgtmolfrags = []
        print(molnames_perrec[i])
        for j in range(len(molnames_perrec[i])):
            try:
                if molnames_perrec[i][j] == self.tgt2molname:
                    tgtmolfrags += molfragss[j]
            except:
                continue

        frag1 = self.tgt1frag
        if type(frag1) == int:
            frag1s = [frag1]
            print('tgtfrag1s', frag1s)
        else:
            frag1s = copy.deepcopy(frag1)

        self.frag1s = frag1s
        # if ifdf.empty:
            # continue
        ifdf_filters = pd.DataFrame()
        for frag1p in frag1s:
            for tgt2frag in tgtmolfrags:
                ifdf_filters = ifdf_filters.append(self.gettgtdf_ff(ifdf, frag1p, tgt2frag))
            print('ifdf_filters\n', ifdf_filters)

        HF_IFIE_sum = ifdf_filters['HF-IFIE'].sum()
        MP2_IFIE_sum = ifdf_filters['MP2-IFIE'].sum()
        PR_TYPE1_sum = ifdf_filters['PR-TYPE1'].sum()
        GRIMME_sum = ifdf_filters['GRIMME'].sum()
        JUNG_sum = ifdf_filters['JUNG'].sum()
        HILL_sum = ifdf_filters['HILL'].sum()

        ifdf_filters['TIMES'] = self.tgttimes[i]

        ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum], index=self.ifdfsumcolumn, name=self.tgttimes[i])

        print('ifdfsum\n', ifdfsum)

        # pieda
        # if pidf.empty:
            # continue
        pidf_filters = pd.DataFrame(columns=self.pcolumn)
        for frag1p in self.frag1s:
            for tgt2frag in tgtmolfrags:
                pidf_filter = self.gettgtdf_ff(pidf, frag1p, tgt2frag)
                if len(pidf_filter) == 0:
                    esdata = self.gettgtdf_ff(ifdf, frag1p, tgt2frag)['HF-IFIE'].values[0]
                    tmpps = pd.Series([frag1p, self.tgt2frag, esdata, 0.0, 0.0, 0.0, 0.0], index=self.pcolumn, name=self.tgttimes[i])
                    # print(tmpps.index)
                    pidf_filters = pidf_filters.append(tmpps)
                else:
                    pidf_filters = pidf_filters.append(pidf_filter)

        ES_sum = pidf_filters['ES'].sum()
        EX_sum = pidf_filters['EX'].sum()
        CT_sum = pidf_filters['CT-mix'].sum()
        DI_sum = pidf_filters['DI(MP2)'].sum()
        q_sum = pidf_filters['q(I=>J)'].sum()

        pidf_filters['TIMES'] = self.tgttimes[i]

        pidfsum = pd.Series([ES_sum, EX_sum, CT_sum, DI_sum, q_sum], index=self.pidfsumcolumn, name=self.tgttimes[i])

        print(pidfsum)

        return ifdf_filters, pidf_filters, ifdfsum, pidfsum


    def getfiltifpifd(self, i, ifdf, pidf):
        # in class var: tgttime
        #    local var: i, ifdf,pidf
        # out local var: tgtifdfsum, tgtdf_filter
        #                pitgtdfsum, pitgtedf

        print('filter-dist:', self.dist)
        # get tgt frag id
        ifdf, ifdf_filter = self.gettgtdf_fd(ifdf)

        HF_IFIE_sum = ifdf_filter['HF-IFIE'].sum()
        MP2_IFIE_sum = ifdf_filter['MP2-IFIE'].sum()
        PR_TYPE1_sum = ifdf_filter['PR-TYPE1'].sum()
        GRIMME_sum = ifdf_filter['GRIMME'].sum()
        JUNG_sum = ifdf_filter['JUNG'].sum()
        HILL_sum = ifdf_filter['HILL'].sum()

        ifdf_filter['TIMES'] = self.tgttimes[i]

        ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum], index=self.ifdfsumcolumn, name=self.tgttimes[i])

        print(ifdfsum)

        # pieda
        frags = self.frags
        pidf_filter = pd.DataFrame()

        pidf_filter = self.getpitgtdf(pidf, ifdf_filter)
        if self.fragmode != 'manual':
            # print('len_frags', len(frags))
            #assign resname(e.g. Gly6)
            for j in range(1, len(frags) + 1):
                pidf_filter.I = pidf_filter.I.replace(j, frags[j-1])
                pidf_filter.J = pidf_filter.J.replace(j, frags[j-1])

        ES_sum = pidf_filter['ES'].sum()
        EX_sum = pidf_filter['EX'].sum()
        CT_sum = pidf_filter['CT-mix'].sum()
        DI_sum = pidf_filter['DI(MP2)'].sum()
        q_sum = pidf_filter['q(I=>J)'].sum()

        pidf_filter['TIMES'] = self.tgttimes[i]

        pidfsum = pd.Series([ES_sum, EX_sum, CT_sum, DI_sum, q_sum], index=self.pidfsumcolumn, name=self.tgttimes[i])

        return ifdf_filter, pidf_filter, ifdfsum, pidfsum


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

        pidf = pd.DataFrame(columns=self.pcolumn)
        pidf['I'] = copy.deepcopy(pii[:pipair[0]])
        pidf['J'] = copy.deepcopy(pij[:pipair[0]])
        pidf['ES'] = copy.deepcopy(es[:pipair[0]])
        pidf['EX'] = copy.deepcopy(ex[:pipair[0]])
        pidf['CT-mix'] = copy.deepcopy(ct[:pipair[0]])
        pidf['DI(MP2)'] = copy.deepcopy(di[:pipair[0]])
        pidf['q(I=>J)'] = copy.deepcopy(qval[:pipair[0]])

        return [ifdf, pidf]


    def read_ifpimulti(self, i):

        # read ifpi f90
        if self.f90soflag == True:
            dfs = self.read_ifpif90(self.tgtlogs[i])
        else:
            dfs = self.read_ifiepiedas(self.tgtlogs[i])
        if len(dfs) == 0:
            return None

        ifdf, pidf = dfs
        # skip = []
        # if len(ifdf) == 0:
            # return

        # IFIE pieda filter
        tgt2type = self.tgt2type
        print('tgttimes', self.tgttimes[i])

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
                for tmptgt1 in self.tgt1frag:
                    hfdf, mp2corrdf, prmp2corrdf, mp2tdf, prmp2tdf = self.gettgtdf_n2tfmatrix(i, ifdf, tmptgt1)
                    esdf, exdf, ctdf = self.gettgtpidf_n2tfmatrix(i, pidf, hfdf, tmptgt1)
                    hfdfs.append(hfdf)
                    mp2corrdfs.append(mp2corrdf)
                    prmp2corrdfs.append(prmp2corrdf)
                    mp2tdfs.append(mp2tdf)
                    prmp2tdfs.append(prmp2tdf)
                    esdfs.append(esdf)
                    exdfs.append(exdf)
                    ctdfs.append(ctdf)
                return hfdfs, mp2corrdfs, prmp2corrdfs, mp2tdfs, prmp2tdfs, esdfs, exdfs, ctdfs

            else:
                ifdf_filter, pidf_filter = self.getfiltifpiff(i, ifdf, pidf)
                iftgtdfsum = None
                pitgtdfsum = None
                return ifdf_filter, pidf_filter

        if tgt2type == 'molname':
            ifdf_filter, pidf_filter, ifdfsum, pidfsum = self.getfiltifpifm(i, ifdf, pidf)
            return ifdf_filter, pidf_filter, ifdfsum, pidfsum

        if tgt2type == 'dist':
            ifdf_filter , pidf_filter, ifdfsum, pidfsum = self.getfiltifpifd(i, ifdf, pidf)
            return ifdf_filter , pidf_filter, ifdfsum, pidfsum

        else:
            print('tgt2type error!!')
            return


    def readifiewrap(self, item1=None, item2=None, item3=None):
        self.setupreadparm(item1, item2, item3)
        if self.anlmode == 'multi':
            print('## read multi mode')
            ifdfs = []
            skips = []
            args = []
            st = time.time()
            p = Pool(self.pynp)
            # for i in range(len(self.tgtlogs)):
            # args.append [self.tgtlogs[i], self.tgttimes[i], self.tgtpdbs[i]]
            if self.f90soflag == True:
                print("use fortran library")

            ifpidfs = p.map(self.read_ifpimulti, [i for i in range(len(self.tgtlogs))])
             # i: time j:type, k:tgt1frag
            # [[hfdf[time 1][frag 1, frag2...], , mp2df[time 1], ...], [hfdf[time 2], mp2df[time 2], ...], ...]
            delids = []
            for i in range(len(ifpidfs)):
                if ifpidfs[i] == None:
                    delids.append(i)
            dellist = lambda items, indexes: [item for index, item in enumerate(items) if index not in indexes]
            ifpidfs = dellist(ifpidfs, delids)

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
                        self.hfdfs.append(hfdf)
                        self.mp2corrdfs.append(mp2corrdf)
                        self.prmp2corrdfs.append(prmp2corrdf)
                        self.mp2tdfs.append(mp2tdf)
                        self.prmp2tdfs.append(prmp2tdf)
                        self.esdfs.append(esdf)
                        self.exdfs.append(exdf)
                        self.ctdfs.append(ctdf)

                else:
                    self.ifdf_filters = pd.DataFrame()
                    self.pidf_filters = pd.DataFrame()
                    for dfs in ifpidfs:
                        self.ifdf_filters = self.ifdf_filters.append(dfs[0])
                        self.pidf_filters = self.pidf_filters.append(dfs[1])


            if self.tgt2type == 'molname' or self.tgt2type == 'dist':
                self.ifdf_filters = pd.DataFrame()
                self.pidf_filters = pd.DataFrame()
                self.ifdfsum = pd.DataFrame()
                self.pidfsum = pd.DataFrame()

                for dfs in ifpidfs:
                    print(dfs[0])
                    self.ifdf_filters = self.ifdf_filters.append(dfs[0])
                    self.pidf_filters = self.pidf_filters.append(dfs[1])
                    self.ifdfsum = self.ifdfsum.append(dfs[2])
                    self.pidfsum = self.pidfsum.append(dfs[3])

            p.close()
            print('read elapsed', time.time() - st)



        else:
            print('## read single mode')
            self.logMethod = self.getlogmethod(self.tgtlogs)
            # self.logMethod = 'MP2'
            if self.matrixtype != 'frags-frags' and self.logMethod == 'MP3':
                print('Error: MP3 mode for this mode is unsupported yet.')
                sys.exit()

            if self.f90soflag == True:
                print("use fortran library")
                ifpidfs = self.read_ifpif90(self.tgtlogs)
                self.ifdf = ifpidfs[0]
                self.pidf = ifpidfs[1]

            else:
                ifie, pieda = self.read_ifiepieda(self.tgtlogs)
                df = self.getifiedf(ifie)
                self.ifdf = df

                pidf = self.getpiedadf(pieda)
                self.pidf = pidf

        return self


    def getlogmethod(self, tgtlog):
        f = open(tgtlog, 'r')
        for line in f:
            Items = line.split()
            if len(Items) <= 2:
                continue
            if Items[0:2] == ['Method', '=']:
                print('logMethod =', Items[2])
                return Items[2]


    def setupreadparm(self, item1=None, item2=None, item3=None):
        print('anlmode:' ,self.anlmode)
        print('fragmode:', self.fragmode)
        print('np =:', self.pynp)

        tgtlogs = []
        tgtpdbs = []
        tgttimes = []

        if self.anlmode == 'multi' and self.tgt2type == 'molname':
            self.rpdbflag = True

        # print(self.ilog_head)
        if self.anlmode == 'multi':
            print('tgt2type:' ,self.tgt2type)

            # setup tgttimes, logs, and pdbs
            for i in range(self.start, self.end+1, self.interval):
                tgttimes.append(str(i))
                tgtlogs.append(self.ilog_head + str(i) + self.ilog_tail)
                if self.rpdbflag == True:
                    tgtpdbs.append(self.pdb_head + str(i) + self.pdb_tail)
            print('tgtlogs', tgtlogs)

            # setup tgt1frag
            if item1 != None:
                print('type', type(item1))
                if type(item1) == str:
                    if '-' in item1:
                        tgt = item1.split('-')
                        print('tgt', tgt)
                        self.tgt1frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                        if self.tgt1frag in self.tgt1frag:
                            del self.tgt1frag[self.tgt1frag.index(self.tgt1frag)]
                    else:
                        self.tgt1frag = [eval(item1)]
                else:
                        self.tgt1frag = item1

            # setup tgt2frag
            if self.tgt2type == 'frag':
                if item2 != None:
                    print('type', type(item2))
                    if type(item2) == str:
                        if '-' in item2:
                            tgt = item2.split('-')
                            print('tgt', tgt)
                            self.tgt2frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                            if self.tgt1frag in self.tgt2frag:
                                del self.tgt2frag[self.tgt2frag.index(self.tgt1frag)]
                        else:
                            self.tgt2frag = [eval(item2)]
                    else:
                            self.tgt2frag = [item2]
                    if type(self.tgt2frag) == list:
                        for dfrag in self.exceptfrag:
                            try:
                                del self.tgt2frag[self.tgt2frag.index(dfrag)]
                                print('del frag', dfrag, 'from tgt2')
                            except:
                                pass
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
                    self.tgt1frag = [item2]
            if self.tgt2type == 'frag':
                if item2 != None:
                    print('type', type(item2))
                    if type(item2) == str:
                        if '-' in item2:
                            tgt = item2.split('-')
                            print('tgt1', tgt)
                            self.tgt1frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                    else:
                        self.tgt1frag = [item2]
                if item3 != None:
                    print('type', type(item3))
                    if type(item3) == str:
                        if '-' in item3:
                            tgt = item3.split('-')
                            print('tgt2', tgt)
                            self.tgt2frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                    else:
                        self.tgt2frag = item3

        if self.matrixtype == 'frags-frags':
            print(self.tgt1frag, self.tgt2frag)
            dp = set(self.tgt1frag) & set(self.tgt2frag)
            if len(dp) != 0:
                print('Error! tgt1 and tgt2 is duplicate')
                sys.exit()

        ## read pdb
        if self.anlmode == 'multi' and self.tgt2type == 'molname':
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
            print('resnames', self.resnames)
            print('len resnames', len(self.resnames))

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

    ### filter section
    def filterifiewrap(self, dist=None):

        tgt2type = self.tgt2type
        if dist != None:
            self.dist = dist
        # frag mode
        if self.anlmode == 'frag':
            if self.matrixtype == 'frags-frags':
                self.gettgtdf_n2ffmatrix()

            else:
                if self.tgt2type == 'dist':
                    # print(self.ifdf)
                    tgtdf, ifdf_filter = self.gettgtdf_fd(self.ifdf)
                    # print(ifdf_filter)
                    self.ifdf_filter = ifdf_filter
                elif self.tgt2type == 'frag':
                    self.ifdf_filters = []
                    if type(self.tgt1frag) == int:
                        self.tgt1frag == [self.tgt1frag]
                    if type(self.tgt2frag) == int:
                        self.tgt1frag == [self.tgt2frag]
                    for tgt1 in self.tgt1frag:
                        self.ifdf_filters.append(self.gettgtdf_ffs(self.ifdf, tgt1, self.tgt2frag))

        # mol-mol mode
        if self.anlmode == 'mol':
            df = self.ifdf
            if self.selecttype == 'fragid':
                molfrags = self.getmolfrags(self.tgt1frag[0], df)
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
            print('nf', nf)
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

        return self


    def filterpiedawrap(self):
        ### filter
        # frag-frag mode
        frags = self.frags
        if self.anlmode == 'frag':
            if  self.matrixtype == 'frags-frags':
                self.gettgtpidf_n2ffmatrix()
            else:

                if self.tgt2type == 'dist':
                    pitgtdf = self.getpitgtdf(self.pidf, self.ifdf_filter)
                    if self.fragmode != 'manual':
                        # print('len_frags', len(frags))
                        #assign resname(e.g. Gly6)
                        for i in range(1, len(frags) + 1):
                            pitgtdf.I = pitgtdf.I.replace(i, frags[i-1])
                            pitgtdf.J = pitgtdf.J.replace(i, frags[i-1])

                    self.pidf_filter = pitgtdf

                elif self.tgt2type == 'frag':
                    t1cnt = 0
                    pidf_filters = []
                    count = 0
                    for tgt1 in self.tgt1frag:
                        pidf_filter = pd.DataFrame(columns = self.pcolumn)
                        for tgt2 in self.tgt2frag:
                            pidf1_filter = self.gettgtdf_ff(self.pidf, tgt1, tgt2)

                            df = self.ifdf_filters[t1cnt]
                            if len(pidf1_filter) == 0:
                                esdata = df[((df['I'] == tgt1) & ( df['J'] == tgt2))|((df['J'] == tgt1) & (df['I'] == tgt2))]['HF-IFIE'].values.tolist()[0]
                                pidf_series = pd.Series([tgt1, tgt2, esdata, 0.0, 0.0, 0.0, 0.0], index=self.pcolumn, name=count)
                                pidf_filter = pidf_filter.append(pidf_series)
                                count += 1
                            else:
                                pidf_filter = pidf_filter.append(pidf1_filter)
                        t1cnt +=1
                        pidf_filters.append(pidf_filter)
                    self.pidf_filters = pidf_filters

            # multi mode

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
        if self.anlmode == 'multi':
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
                                self.prmp2tdf
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
                            ]
                elif self.logMethod == 'MP3':
                    datadfs = [
                                self.esdf,
                                self.exdf,
                                self.ctdf,
                                self.hfdf,
                                self.mp3corrdf,
                                self.usermp3corrdf,
                                self.mp3tdf,
                                self.usermp3tdf
                            ]
                    names = [
                                'ES',
                                'EX',
                                'CT',
                                'HF',
                                'MP3corr',
                                'USER-MP3corr',
                                'MP3total',
                                'USER-MP3total',
                            ]

                tgt1str = str(self.tgt1frag[0]) + '-'  + str(self.tgt1frag[-1])
                tgt2str = str(self.tgt2frag[0]) + '-'  + str(self.tgt2frag[-1])
                for i in range(len(datadfs)):
                    ocsv = head + '_frag' + str(tgt1str) + '-frag' + str(tgt2str) + '-' + names[i] + '-ffmatrix.csv'
                    datadfs[i].T.to_csv(path + '/' + ocsv)
                    print(path + '/' + ocsv + ' was created.')

            else:
                if self.tgt2type == 'dist':
                    tgtid = self.tgt1frag[0]
                    try:
                        ohead = head + '-' + str(tgtid) + '-' + frags[tgtid - 1]
                    except:
                        ohead = head + '-' + str(tgtid)

                    # print(self.pidf_filter.head())

                    # tgtdf.to_csv(path + '/' + ohead + '-ifie.csv')
                    oifie = ohead + '-ifie_' + 'dist' + str(self.dist) + '.csv'
                    opieda = ohead + '-pieda_' + 'dist' + str(self.dist) + '.csv'
                    self.ifdf_filter.to_csv(path + '/' + oifie)
                    self.pidf_filter.to_csv(path + '/' + opieda)

                    print(path + '/' + oifie, path + '/' + opieda, 'was generated.')

                elif self.tgt2type == 'frag':

                    tgt1s = self.tgt1frag
                    count = 0
                    for tgt1 in tgt1s:
                        ohead = head + '-frag' + str(tgt1) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1])
                        # print(self.pidf_filter.head())
                        # tgtdf.to_csv(path + '/' + ohead + '-ifie.csv')
                        oifie = ohead + '-ifie.csv'
                        opieda = ohead + '-pieda.csv'
                        self.ifdf_filters[count].to_csv(path + '/' + oifie)
                        self.pidf_filters[count].to_csv(path + '/' + opieda)
                        print(path + '/' + oifie, path + '/' + opieda, 'was generated.')
                        count += 1
                    # N:1 sheet
                    for j in range(len(self.ifdf_filters)):
                        if j == 0:
                            ifdf_fil_n1 = self.ifdf_filters[j]
                            pidf_fil_n1 = self.pidf_filters[j]
                        else:
                            ifdf_fil_n1 += self.ifdf_filters[j].values
                            pidf_fil_n1 += self.pidf_filters[j].values

                    ohead = head + '-frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1]) + 'n-1sum'

                    oifie = ohead + '-ifie.csv'
                    opieda = ohead + '-pieda.csv'

                    ifdf_fil_n1 = ifdf_fil_n1.reset_index(drop=True)
                    pidf_fil_n1 = pidf_fil_n1.reset_index(drop=True)

                    ifdf_fil_n1["I"] = str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1])
                    ifdf_fil_n1["J"] = self.tgt2frag

                    pidf_fil_n1['I'] = str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1])
                    pidf_fil_n1['J'] = self.tgt2frag
# 
                    del ifdf_fil_n1['DIMER-ES']
                    del ifdf_fil_n1['DIST']

                    ifdf_fil_n1.to_csv(path + '/' + oifie)
                    pidf_fil_n1.to_csv(path + '/' + opieda)

                    print(path + '/' + oifie, path + '/' + opieda, 'was generated.')

        if self.anlmode == 'multi':
            head = self.ilog_head
            if tgt2type == 'frag':

                if self.matrixtype == 'times-frags':
                    datadfs = [
                                self.esdfs,
                                self.exdfs,
                                self.ctdfs,
                                self.hfdfs,
                                self.mp2corrdfs,
                                self.prmp2corrdfs,
                                self.mp2tdfs,
                                self.prmp2tdfs
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
                            ]
                    tgt2str = str(self.tgt2frag[0]) + '-'  + str(self.tgt2frag[-1])

                    print('--- out files ---')
                    # i:type, j:tgt1frag
                    for i in range(len(datadfs)):
                        sum2df = pd.DataFrame(index=datadfs[i][0].columns)
                        for j in range(len(datadfs[i])):
                            tgt1frag = self.tgt1frag[j]
                            ocsv = head + '_frag' + str(tgt1frag) + '-frag' + str(tgt2str) + '-' + names[i] + '-tfmatrix.csv'
                            datadfs[i][j].T.to_csv(path + '/' + ocsv)
                            print(path + '/' + ocsv)

                            # gen tgtfrag1 sum matrix
                            if j == 0:
                                sum1df =  datadfs[i][j]
                            else:
                                sum1df = sum1df + datadfs[i][j]

                            # gen tgtfrag2 sum matrix
                            sum2dfbuf = datadfs[i][j].sum()
                            sum2dfbuf.name = tgt1frag
                            sum2df = pd.concat([sum2df, sum2dfbuf], axis=1)

                        osum1csv = head + '_frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(tgt2str) + '-' + names[i] + '-sumtfmatrix.csv'
                        sum1df.T.to_csv(path + '/' + osum1csv)
                        print(path + '/' + osum1csv)

                        osum2csv = head + '_frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(tgt2str) + '-' + names[i] + '-sum2matrix.csv'
                        sum2df.to_csv(path + '/' + osum2csv)
                        print(path + '/' + osum2csv)

                else:
                    tgt1frag = self.tgt1frag[0]
                    tgt2frag = self.tgt2frag[0]
                    oifie = 'frag' + str(tgt1frag) + '-frag' + str(tgt2frag) + '-ifie.csv'
                    opieda = 'frag' + str(tgt1frag) + '-frag' + str(tgt2frag) + '-pieda.csv'
                    self.ifdf_filters.to_csv(path + '/' + oifie)
                    self.pidf_filters.to_csv(path + '/' + opieda)
                    print(path + '/' + oifie, path + '/' + opieda, 'was created.')


            if tgt2type == 'molname' or tgt2type == 'dist':

                tgt1frag = self.tgt1frag[0]

                if tgt2type == 'molname':
                    tgt2molname = self.tgt2molname
                    oifie = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-ifiesum.csv'
                    opieda = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-piedasum.csv'
                    oifiedt = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-ifiedt.csv'
                    opiedadt = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-piedadt.csv'

                if tgt2type == 'dist':
                    tgt2dist = self.dist
                    oifie = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-ifiesum.csv'
                    opieda = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) +'-piedasum.csv'
                    oifiedt = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-ifiedt.csv'
                    opiedadt = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-piedadt.csv'


                self.ifdfsum.to_csv(path + '/' + oifie)
                self.pidfsum.to_csv(path + '/' + opieda)
                self.ifdf_filters.to_csv(path + '/' + oifiedt)
                self.pidf_filters.to_csv(path + '/' + opiedadt)

                print(path + '/' + oifie, path + '/' + opieda)
                print(path + '/' + oifiedt, path + '/' + opiedadt, 'was created.')



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
                tgtid = self.tgt1frag[0]

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






