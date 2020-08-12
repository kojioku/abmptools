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
        self.ifdfsumcolumn = ['HF-IFIE', 'MP2-IFIE', 'PR-TYPE1', 'GRIMME', 'JUNG', 'HILL', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']
        self.logMethod = 'MP2'


        self.anlmode= 'frag' #frag, 'mol', 'fraginmol', 'multi'
        self.fragmode = 'auto'  #'hybrid', 'auto', 'manual'
        self.dist = 1000.0
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
        self.addresinfo = True
        self.writeresnamecsv = True


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

            if items [0:2] == ['START', 'FRAGMENT']:
                break

            ## AUTOMATIC FRAGMENTATION
            if items[0:3] == ['Seq.', 'Frag.', 'Residue']:
                autoreadflag = True
                continue

            if autoreadflag == True and items[0:2] == ['The', 'system']:
                autoreadflag = False
                continue

            if autoreadflag == True:
               print(items)
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
        logabsitems = os.path.abspath(ifile).split('/')
        logabsitems[-1] = logreadGeom
        # print(logabsitems)
        pdbabs = ""
        for logabsitem in logabsitems:
            pdbabs = pdbabs + logabsitem + '/'

        pdbabs = pdbabs[:-1]
        # print(pdbabs)

        # print('Frag Atom number\n', fragdatas)

        resname_perfrag = []
        if self.fragmode == 'manual':
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
                    alreadys.append(fragnos[i])
            # print(resname_perfrag)

        return resname_perfrag, pdbabs


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
        distdf = pd.DataFrame(index=self.tgt2frag)

        df = self.ifdf
        count = 0

        if self.logMethod == 'MP2':
            mp2corrdf = pd.DataFrame(index=self.tgt2frag)
            prmp2corrdf = pd.DataFrame(index=self.tgt2frag)
            mp2tdf =  pd.DataFrame(index=self.tgt2frag)
            prmp2tdf = pd.DataFrame(index=self.tgt2frag)

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
                dist = tgtdf_filter['DIST'].values.tolist()

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
                distdf[str(f1)] = dist

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
            ifdf_frag_mols.append(ifie_frag_mol)

        #pieda
        for i in range(len(ifdf_frag_mols)):
            ifdf_frag_mols[i] =  pd.merge(ifdf_frag_mols[i], self.pidf, on=['I', 'J'], how='left')
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

        ifdf_mol_mol = pd.DataFrame(columns=self.ifdfsumcolumn)
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
            # print('ifdf_filters\n', ifdf_filters)

        ifdf_filters = pd.merge(ifdf_filters, pidf, on=['I', 'J'], how='left')

        HF_IFIE_sum = ifdf_filters['HF-IFIE'].sum()
        MP2_IFIE_sum = ifdf_filters['MP2-IFIE'].sum()
        PR_TYPE1_sum = ifdf_filters['PR-TYPE1'].sum()
        GRIMME_sum = ifdf_filters['GRIMME'].sum()
        JUNG_sum = ifdf_filters['JUNG'].sum()
        HILL_sum = ifdf_filters['HILL'].sum()

        ES_sum = ifdf_filters['ES'].sum()
        EX_sum = ifdf_filters['EX'].sum()
        CT_sum = ifdf_filters['CT-mix'].sum()
        DI_sum = ifdf_filters['DI(MP2)'].sum()
        q_sum = ifdf_filters['q(I=>J)'].sum()

        ifdf_filters['TIMES'] = self.tgttimes[i]

        ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum], index=self.ifdfsumcolumn, name=self.tgttimes[i])

        print('ifdfsum\n', ifdfsum)

        # pieda
        # if pidf.empty:
            # continue
#         pidf_filters = pd.DataFrame(columns=self.pcolumn)
#         for frag1p in self.frag1s:
#             for tgt2frag in tgtmolfrags:
#                 pidf_filter = self.gettgtdf_ff(pidf, frag1p, tgt2frag)
#                 if len(pidf_filter) == 0:
#                     esdata = self.gettgtdf_ff(ifdf, frag1p, tgt2frag)['HF-IFIE'].values[0]
#                     tmpps = pd.Series([frag1p, self.tgt2frag, esdata, 0.0, 0.0, 0.0, 0.0], index=self.pcolumn, name=self.tgttimes[i])
#                     # print(tmpps.index)
#                     pidf_filters = pidf_filters.append(tmpps)
#                 else:
#                     pidf_filters = pidf_filters.append(pidf_filter)


        # pidf_filters['TIMES'] = self.tgttimes[i]

        # pidfsum = pd.Series([ES_sum, EX_sum, CT_sum, DI_sum, q_sum], index=self.pidfsumcolumn, name=self.tgttimes[i])

        return ifdf_filters, ifdfsum


    def getfiltifpifd(self, i, ifdf, pidf):
        # in class var: tgttime
        #    local var: i, ifdf,pidf
        # out local var: tgtifdfsum, tgtdf_filter
        #                pitgtdfsum, pitgtedf

        print('filter-dist:', self.dist)
        # get tgt frag id
        ifdf, ifdf_filter = self.gettgtdf_fd(ifdf)

        ifdf_filter = pd.merge(ifdf_filter, pidf, on=['I', 'J'], how='left')

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

        ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum], index=self.ifdfsumcolumn, name=self.tgttimes[i])
        print(ifdfsum)

        # pieda
#         frags = self.frags
#         if self.fragmode != 'manual':
#             # print('len_frags', len(frags))
#             #assign resname(e.g. Gly6)
#             for j in range(1, len(frags) + 1):
#                 ifdf_filter.I = ifdf_filter.I.replace(j, frags[j-1])
#                 ifdf_filter.J = ifdf_filter.J.replace(j, frags[j-1])

        return ifdf_filter, ifdfsum


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
                ifdf_filter = self.getfiltifpiff(i, ifdf, pidf)
                iftgtdfsum = None
                return ifdf_filter

        if tgt2type == 'molname':
            ifdf_filter, ifdfsum = self.getfiltifpifm(i, ifdf, pidf)
            return ifdf_filter, ifdfsum

        if tgt2type == 'dist':
            ifdf_filter, ifdfsum  = self.getfiltifpifd(i, ifdf, pidf)
            return ifdf_filter, ifdfsum

        else:
            print('tgt2type error!!')
            return


    def readifiewrap(self, item1=None, item2=None, item3=None):
        self.setupreadparm(item1, item2, item3)

        # read fraginfo

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

            pd.set_option('display.width', 500)
            print(len(ifpidfs))

            delids = []
            for i in range(len(ifpidfs)):
                try:
                    if ifpidfs[i] == None:
                        delids.append(i)
                except:
                    pass
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
                    for dfs in ifpidfs:
                        self.ifdf_filters = self.ifdf_filters.append(dfs)


            if self.tgt2type == 'molname' or self.tgt2type == 'dist':
                self.ifdf_filters = pd.DataFrame()
                self.ifdfsum = pd.DataFrame(columns=self.ifdfsumcolumn)

                for dfs in ifpidfs:
                    self.ifdf_filters = self.ifdf_filters.append(dfs[0])
                    self.ifdfsum = self.ifdfsum.append(dfs[1])

            p.close()
            print('read elapsed', time.time() - st)



        else:
            print('## read single mode')
            self.logMethod = self.getlogmethod(self.tgtlogs)
            # self.logMethod = 'MP2'
            if self.matrixtype != 'frags-frags' and (self.logMethod == 'MP3' or self.logMethod == 'CCPT'):
                print('Error: ' + self.logMethod + ' mode for this mode is unsupported yet.')
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

        tgtlogs = []
        tgtpdbs = []
        tgttimes = []

        if self.anlmode == 'multi' and self.tgt2type == 'molname':
            self.rpdbflag = True

        # multi mode
        if self.anlmode == 'multi':
            print('tgt2type:' ,self.tgt2type)

            # setup tgttimes, logs, and pdbs
            for i in range(self.start, self.end+1, self.interval):
                tgttimes.append(str(i))
                tgtlogs.append(self.ilog_head + str(i) + self.ilog_tail)
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
                                print('- Info: del frag', dfrag, 'from tgt2')
                            except:
                                pass
                print('tgt1frag, tgt2frag', self.tgt1frag, self.tgt2frag)

            if self.tgt2type == 'molname':
                if item2 != None:
                    self.tgt2molname = item2
                print('tgt1frag, tgt2mol', self.tgt1frag, self.tgt2molname)

            self.tgtlogs = tgtlogs
            self.tgttimes = tgttimes

            if self.rpdbflag == True:
                nfs = []
                molnames_perrec = []
                self.assignmolname = False
                for i in range(len(tgtlogs)):
                    self.resname_perfrag, tgtpdb = self.getlogorpdbfrag(self.tgtlogs[i])
                    tgtpdbs.append(tgtpdb)
                    if self.fragmode == 'auto':
                        print('Error: auto fragment in mol mode is not suppoted yet.')
                        sys.exit()
                    nf = self.getlognf(tgtlogs[i], self.fragmode)
                    nfs.append(nf)
                    molnames_perrec.append(self.resnames)
                self.nfs = nfs
                self.molnames_perrec  = molnames_perrec
                self.tgtpdbs = tgtpdbs

        # single mode
        else:
            if item1 != None:
                self.tgtlogs = item1
            if item2 != None:
                if self.anlmode == 'mol' and self.selecttype == 'molid':
                    self.tgtmolid = int(item2)
                else:
                    self.tgt1frag = [item2]
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

            if self.tgt2type == 'frag':
                if item2 != None:
                    # print('type tgt2', type(item2))
                    if type(item2) == str:
                        if '-' in item2:
                            tgt = item2.split('-')
                            print('tgt1', tgt)
                            self.tgt1frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                        else:
                            self.tgt1frag = [eval(item2)]

                    else:
                        self.tgt1frag = [item2]
                if item3 != None:
                    # print('type tgt2', type(item3))
                    if type(item3) == str:
                        if '-' in item3:
                            tgt = item3.split('-')
                            print('tgt2', tgt)
                            self.tgt2frag = [ i for i in range(int(tgt[0]), int(tgt[1]) + 1) ]
                        else:
                            self.tgt2frag = [eval(item3)]

                    else:
                        self.tgt2frag = [item3]


        if self.matrixtype == 'frags-frags':
            print(self.tgt1frag, self.tgt2frag)
            dp = set(self.tgt1frag) & set(self.tgt2frag)
            if len(dp) != 0:
                print('Error! tgt1 and tgt2 is duplicate')
                sys.exit()

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
            print('read reference resname')
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


    ### filter section
    def filterifiewrap(self, dist=None):

        tgt2type = self.tgt2type
        if dist != None:
            self.dist = dist
        # frag mode
        if self.anlmode == 'frag':
            if self.matrixtype == 'frags-frags':
                self.gettgtdf_n2ffmatrix()
                self.gettgtpidf_n2ffmatrix()

            else:
                if self.tgt2type == 'dist':
                    # print(self.ifdf)
                    tgtdf, ifdf_filter = self.gettgtdf_fd(self.ifdf)
                    self.ifdf_filter = pd.merge(ifdf_filter, self.pidf, on=['I', 'J'], how='left')
                    print(self.ifdf_filter)
                    # pitgtdf = self.getpitgtdf(self.pidf, self.ifdf_filter)
#                     if self.fragmode != 'manual':
#                         #assign resname(e.g. Gly6)
#                         for i in range(1, len(frags) + 1):
#                             pitgtdf.I = pitgtdf.I.replace(i, frags[i-1])
#                             pitgtdf.J = pitgtdf.J.replace(i, frags[i-1])

                elif self.tgt2type == 'frag':
                    self.ifdf_filters = []
                    if type(self.tgt1frag) == int:
                        self.tgt1frag == [self.tgt1frag]
                    if type(self.tgt2frag) == int:
                        self.tgt1frag == [self.tgt2frag]
                    for tgt1 in self.tgt1frag:
                        ifdf_filter = self.gettgtdf_ffs(self.ifdf, tgt1, self.tgt2frag)
                        # self.ifdf_filters.append(self.gettgtdf_ffs(self.ifdf, tgt1, self.tgt2frag))
                        self.ifdf_filters.append(pd.merge(ifdf_filter, self.pidf, on=['I', 'J'], how='left'))

                    print(self.ifdf_filters[0].head())

        # mol-mol mode
        if self.anlmode == 'mol':
            #ifie
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
            ifdfmolsums = pd.DataFrame(columns=self.ifdfsumcolumn)

            for i in range(len(self.tgtmolfrags)):
                contactmolfrags, ifdf_frag_mol, ifdfmol_mol, ifdfmolsum = self.getifiesummol(df, tgtmolfrags[i], self.tgtmolid[i])

                # self.contactmolfrags = contactmolfrags
                print(ifdf_frag_mol)
                print(ifdfmol_mol)
                print(ifdfmolsum)

                ifdf_frag_mols = ifdf_frag_mols.append(ifdf_frag_mol)
                ifdfmol_mols = ifdfmol_mols.append(ifdfmol_mol)
                ifdfmolsums = ifdfmolsums.append(ifdfmolsum)

            self.ifdf_frag_mols = ifdf_frag_mols
            self.ifdfmol_mols = ifdfmol_mols
            self.ifdfmolsums = ifdfmolsums

        # fraginmol mode
        if self.anlmode == 'fraginmol':
            ifdf_filters = pd.DataFrame()
            ifdfsums = pd.DataFrame(columns=self.ifdfsumcolumn)

            df = self.ifdf
            tgt1_lofrag = self.tgt1_lofrag
            tgt2_lofrag = self.tgt2_lofrag
            tgt2molname = self.tgt2molname
            nf = self.getlognf(self.tgtlogs, self.fragmode)
            print('nf', nf)
            molfragss = self.getallmolfrags(self.tgtlogs, df, nf)
            print(molfragss)

            tgt2_glofrags = []
            # print('resnames', self.resnames)
            for i in range(len(self.resnames)):
                if self.resnames[i] == tgt2molname:
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
                HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum = self.getsumdf(ifdf_filter)

                ifdfsum = pd.Series([HF_IFIE_sum, MP2_IFIE_sum, PR_TYPE1_sum, GRIMME_sum, JUNG_sum, HILL_sum, ES_sum, EX_sum, CT_sum, DI_sum, q_sum], index=self.ifdfsumcolumn, name='mol' + str(tgtmol+1))
                ifdfsums = ifdfsums.append(ifdfsum)


            # self.tgt1_glofrag = tgt1_glofrags
            # self.tgt2_glofrags = tgt2_glofrags
            self.ifdf_filters = ifdf_filters
            self.ifdfsums = ifdfsums

        return self


    def writecsvwrap(self, head=None):
        # -------------
        # --- write ---
        # -------------
        print('## Write Section')
        path = 'csv'
        tgt2type = self.tgt2type
        if os.path.exists('csv') == False:
            os.mkdir('csv')

        if self.anlmode == 'frag':
            if head == None:
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


                tgt1str = str(self.tgt1frag[0]) + '-'  + str(self.tgt1frag[-1])
                tgt2str = str(self.tgt2frag[0]) + '-'  + str(self.tgt2frag[-1])
                for i in range(len(datadfs)):
                    # rename index-columns
                    if self.addresinfo == True:
                        datadfs[i].rename(index = lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')', \
                                          columns = lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')', inplace=True)

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

                    oifie = ohead + '-ifie_' + 'dist' + str(self.dist) + '.csv'

                    if self.addresinfo == True:
                        for i in range(1, len(self.resname_perfrag)+1):
                            val1 = i
                            val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                            self.ifdf_filter.I = self.ifdf_filter.I.replace(val1, val2)
                            self.ifdf_filter.J = self.ifdf_filter.J.replace(val1, val2)

                        print(self.ifdf_filter.I)

                    self.ifdf_filter.to_csv(path + '/' + oifie)

                    print(path + '/' + oifie, 'was generated.')

                elif self.tgt2type == 'frag':

                    tgt1s = self.tgt1frag
                    count = 0
                    for tgt1 in tgt1s:
                        ohead = head + '-frag' + str(tgt1) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1])
                        oifie = ohead + '-ifie.csv'

                        ifdf_filter = self.ifdf_filters[count]

                        if self.addresinfo == True:
                            for i in range(1, len(self.resname_perfrag)+1):
                                val1 = i
                                val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                                ifdf_filter.I = ifdf_filter.I.replace(val1, val2)
                                ifdf_filter.J = ifdf_filter.J.replace(val1, val2)

                        ifdf_filter.to_csv(path + '/' + oifie)
                        print(path + '/' + oifie, 'was generated.')
                        count += 1
                    # N:1 sheet
                    for j in range(len(self.ifdf_filters)):
                        if j == 0:
                            ifdf_fil_n1 = self.ifdf_filters[j]
                            # pidf_fil_n1 = self.pidf_filters[j]
                        else:
                            ifdf_fil_n1 += self.ifdf_filters[j].values
                            # pidf_fil_n1 += self.pidf_filters[j].values

                    ohead = head + '-frag' + str(self.tgt1frag[0]) + '-' + str(self.tgt1frag[-1]) + '-frag' + str(self.tgt2frag[0]) + '-' + str(self.tgt2frag[-1]) + 'n-1sum'
                    oifie = ohead + '-ifie.csv'

                    ifdf_fil_n1 = ifdf_fil_n1.reset_index(drop=True)


                    if self.addresinfo == True:
                        ifdf_fil_n1["I"] = self.resname_perfrag[int(self.tgt1frag[0])-1] + '(' + str(self.tgt1frag[0]) + ')' + '-' + \
                            self.resname_perfrag[int(self.tgt1frag[-1])-1] + '(' + str(self.tgt1frag[-1]) + ')'
                        ifdf_fil_n1["J"] = self.tgt2frag

                        for i in range(1, len(self.resname_perfrag)+1):
                            val1 = i
                            val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                            ifdf_fil_n1.J = ifdf_fil_n1.J.replace(val1, val2)
                    else:
                        ifdf_fil_n1["I"] = self.tgt1frag[0] + '-' +  self.tgt1frag[-1]
                        ifdf_fil_n1["J"] = self.tgt2frag

                    del ifdf_fil_n1['DIMER-ES']
                    del ifdf_fil_n1['DIST']

                    ifdf_fil_n1.to_csv(path + '/' + oifie)

                    print(path + '/' + oifie, 'was generated.')

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
                            if self.addresinfo == True:
                                datadfs[i][j].rename(index = lambda x: self.resname_perfrag[int(x)-1] + '(' + str(x) + ')', inplace=True)
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

                    if self.addresinfo == True:
                        for i in range(1, len(self.resname_perfrag)+1):
                            val1 = i
                            val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                            self.ifdf_filters.I = self.ifdf_filters.I.replace(val1, val2)
                            self.ifdf_filters.J = self.ifdf_filters.J.replace(val1, val2)

                    oifie = 'frag' + str(tgt1frag) + '-frag' + str(tgt2frag) + '-ifie.csv'
                    self.ifdf_filters.to_csv(path + '/' + oifie)
                    print(path + '/' + oifie, 'was created.')


            if tgt2type == 'molname' or tgt2type == 'dist':

                tgt1frag = self.tgt1frag[0]

                if tgt2type == 'molname':
                    tgt2molname = self.tgt2molname
                    oifie = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-ifiesum.csv'
                    oifiedt = 'frag' + str(tgt1frag) + '-' + str(tgt2molname) + '-ifiedt.csv'

                if tgt2type == 'dist':
                    tgt2dist = self.dist
                    oifie = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-ifiesum.csv'
                    oifiedt = 'frag' + str(tgt1frag) + '-dist' + str(tgt2dist) + '-ifiedt.csv'

                if self.addresinfo == True:
                    for i in range(1, len(self.resname_perfrag)+1):
                        val1 = i
                        val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                        self.ifdf_filters.I = self.ifdf_filters.I.replace(val1, val2)
                        self.ifdf_filters.J = self.ifdf_filters.J.replace(val1, val2)

                self.ifdfsum.to_csv(path + '/' + oifie)
                self.ifdf_filters.to_csv(path + '/' + oifiedt)

                print(path + '/' + oifie)
                print(path + '/' + oifiedt, 'was created.')


        if self.anlmode == 'mol':
            if head == None:
                head, ext = os.path.splitext(self.tgtlogs)

            dist = self.dist
            selecttype = self.selecttype
            # piedamol_mol = self.piedamol_mol
            ifdf_frag_mols = self.ifdf_frag_mols
            ifdfmol_mols= self.ifdfmol_mols

            # pieda_frag_mols = self.pieda_frag_mols
            if selecttype == 'molid':
                tgtid = self.tgtmolid
            else:
                tgtid = self.tgt1frag[0]

            idstr = str(tgtid[0]) + '-' + str(tgtid[-1])
            ilogdtname = path + '/' + head + '_ifie-fragmol-' +  selecttype + idstr + 'dist' + str(dist) + '.csv'
            imolname = path + '/' + head + '_ifiemol-mol-' + selecttype + idstr + 'dist' + str(dist) + '.csv'
            isumname = path + '/' + head + '_ifiesummol-mol-' + selecttype + idstr + 'dist' + str(dist) + '.csv'

#             ifdf_frag_molsdt = pd.DataFrame()
#             pd.set_option('display.width', 500)
#             for ifdf_frag_mol in ifdf_frag_mols:
#                 ifdf_frag_molsdt = ifdf_frag_molsdt.append(ifdf_frag_mol)

            if self.addresinfo == True:
                for i in range(1, len(self.resname_perfrag)+1):
                    val1 = i
                    val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                    ifdf_frag_mols.I = ifdf_frag_mols.I.replace(val1, val2)
                    ifdf_frag_mols.J = ifdf_frag_mols.J.replace(val1, val2)

            # print(ifdf_frag_molsdt, file=ilogdt)
            ifdf_frag_mols.to_csv(ilogdtname)
            ifdfmol_mols.to_csv(imolname)
            self.ifdfmolsums.to_csv(isumname)

            print('---out---')
            print(ilogdtname)
            print(imolname)
            print(isumname)

        if self.anlmode == 'fraginmol':
            if head == None:
                head, ext = os.path.splitext(self.tgtlogs)

            ohead = head + '-' 'tgt1frag' + str(self.tgt1_lofrag) + '-mol' + str(self.tgt2molname) + 'frag' + str(self.tgt2_lofrag)

            if self.addresinfo == True:
                for i in range(1, len(self.resname_perfrag)+1):
                    val1 = i
                    val2 = self.resname_perfrag[i-1] + '(' + str(val1) + ')'
                    self.ifdf_filters.I = self.ifdf_filters.I.replace(val1, val2)
                    self.ifdf_filters.J = self.ifdf_filters.J.replace(val1, val2)
                    # self.ifdf.I = self.ifdf.I.replace(val1, val2)
                    # self.ifdf.J = self.ifdf.J.replace(val1, val2)

            # self.ifdf.to_csv(path + '/' + head + '-ifie.csv')
            oifie = path + '/' + ohead + '-ifie_'  + 'dist' + str(self.dist) + '.csv'
            oifiesum = path + '/' + ohead + '-ifiesum_'  + 'dist' + str(self.dist) + '.csv'
            self.ifdf_filters.to_csv(oifie)
            self.ifdfsums.to_csv(oifiesum)
            print(oifie, 'was generated.')
            print(oifiesum, 'was generated.')



#                     t1cnt = 0
#                     pidf_filters = []
#                     count = 0
#                     for tgt1 in self.tgt1frag:
#                         pidf_filter = pd.DataFrame(columns = self.pcolumn)
#                         for tgt2 in self.tgt2frag:
#                             pidf1_filter = self.gettgtdf_ff(self.pidf, tgt1, tgt2)
#
#                             df = self.ifdf_filters[t1cnt]
#                             if len(pidf1_filter) == 0:
#                                 esdata = df[((df['I'] == tgt1) & ( df['J'] == tgt2))|((df['J'] == tgt1) & (df['I'] == tgt2))]['HF-IFIE'].values.tolist()[0]
#                                 pidf_series = pd.Series([tgt1, tgt2, esdata, 0.0, 0.0, 0.0, 0.0], index=self.pcolumn, name=count)
#                                 pidf_filter = pidf_filter.append(pidf_series)
#                                 count += 1
#                             else:
#                                 pidf_filter = pidf_filter.append(pidf1_filter)
#                         t1cnt +=1
#                         pidf_filters.append(pidf_filter)
#                     self.pidf_filters = pidf_filters

#             pidf = self.pidf
#             tgt1_glofrag = self.tgt1_glofrag
#             tgt2_glofrags = self.tgt2_glofrags
#             print('--- pieda ----')
#             pitgtdf = pidf[pidf['I'] == tgt1_glofrag]
#             pitgtdf = pitgtdf.append(pidf[pidf['J'] == tgt1_glofrag])
#
#             pitgt_new2 = pd.DataFrame()
#             for tgt2_glofrag in tgt2_glofrags:
#                 pitgt_new = pitgtdf[(pitgtdf['I'] == tgt2_glofrag) |(pitgtdf['J'] == tgt2_glofrag)]
#                 pitgt_new2 = pitgt_new2.append(pitgt_new)
#
#             if self.fragmode != 'manual':
#                 print('len_frags', len(frags))
#
#                 for i in range(1, len(frags) + 1):
#                     pitgt_new2.I = pitgt_new2.I.replace(i, frags[i-1])
#                     pitgt_new2.J = pitgt_new2.J.replace(i, frags[i-1])
#             self.pidf_filters = pitgt_new2


#     def readpiedawrap(self, item1=None, item2=None):
#         print('--- pieda ---')
#
#         # self.setupreadparm(item1, item2)
#         if self.abinit_ver == 'rev16' or self.abinit_ver == 'rev17':
#             self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'Solv(ES)', 'DI(MP2)', 'q(I=>J)']
#         else:
#             self.pcolumn = ['I', 'J', 'ES', 'EX', 'CT-mix', 'DI(MP2)', 'q(I=>J)']
#
#         ### read fraginfo section
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
#         ### read pieda (from log) section
#         if self.anlmode == 'multi':
#             pidfs = []
#             for logname in self.tgtlogs:
#                 pieda = self.read_pieda(logname)
#                 if len(pieda) != 0:
#                     pidfs.append(self.getpiedadf(pieda))
#
#             self.pidfs = pidfs
#         else:
#             pieda = self.read_pieda(self.tgtlogs)
#             pidf = self.getpiedadf(pieda)
#
#             self.pidf = pidf
#         return self

