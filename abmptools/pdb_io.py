import numpy as np
import sys
import os
scrdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scrdir)
import math
import subprocess
import re
import time
import copy
import abinit_io as fab
import collections
try:
    from UDFManager import *
except:
    pass

class pdb_io(fab.abinit_io):
    def __init__(self):
        super().__init__()
        # print('## load pdb io init')
        self.solutes = []
        self.getmode = 'resnum'
        self.assignresname = False
        self.refreshatmtype = False
        self.refreshresid = False
        self.cellsize = 0
        self.totalRes = []
        self.atmtypeRes = []
        self.resnames = []
        self.gatmlabRes = []
        self.posRes = []
        self.headRes = []
        self.labRes = []
        self.chainRes = []
        self.resnumRes = []
        self.codeRes = []
        self.occRes = []
        self.tempRes = []
        self.amarkRes = []
        self.chargeRes = []
        self.rescount = []

        pass

    def readpdb(self, fname):

        print('--- get pdbinfo ---')
        print('infile:',  fname)
        lines = open(fname,'r').readlines()
        molnames = []
        poss = []
        atypenames = []
        heads = []
        molnames = []
        atypenames = []
        labs = []
        chains = []
        resnums = []
        codes = []
        poss = []
        occs = []
        temps = []
        amarks = []
        charges = []

        posRes = []
        atmtypeRes = []
        headRes = []
        labRes = []
        chainRes = []
        resnumRes = []
        codeRes = []
        occRes = []
        tempRes = []
        amarkRes = []
        chargeRes = []
        totalRes = 0
        resnames = []
        nums = []

        atomcount = 0
        for line in lines:
            itemlist = line.split()
            # print(itemlist)
            if self.getmode == 'TER':
                if len(itemlist) >= 1 and itemlist[0] == 'TER':
                    print(itemlist)
                    posRes.append(poss)
                    atmtypeRes.append(atypenames)
                    headRes.append(heads)
                    labRes.append(labs)
                    chainRes.append(chains)
                    resnumRes.append(resnums)
                    codeRes.append(codes)
                    occRes.append(occs)
                    tempRes.append(temps)
                    amarkRes.append(amarks)
                    chargeRes.append(charges)
                    poss = []
                    atypenames =[]
                    heads = []
                    labs = []
                    chains = []
                    resnums = []
                    codes = []
                    occs = []
                    temps = []
                    amarks = []
                    charges = []
                    resnames.append(totalRes)
                    totalRes += 1
                    continue

            if len(itemlist) < 3:
                continue

            if line[0:6] == 'HETATM' or line[0:4] == 'ATOM':
            # if itemlist[0] == 'ATOM' or itemlist[0] == 'HETATM':
            #     molname = itemlist[3]
            #     molnames.append(molname)
            #     poss.append([float(e) for e in itemlist[5:8]])
            #     # atypenames.append(itemlist[2])
            #     atypenames.append(line[12:15])

            #     atomcount += 1
            # # over HETATM10000
            # elif len(itemlist[0]) > 6 and itemlist[0][0:6] == 'HETATM':
            #     molname = itemlist[2]
            #     molnames.append(molname)
            #     poss.append([float(e) for e in itemlist[4:7]])
            #     # atypenames.append(itemlist[1])
            #     atypenames.append(line[12:15])
                atomcount += 1

                try:
                    num=int(line[6:12])
                except:
                    head=line[0:6]
                    num=line[6:11]
                    atypename=line[12:16]
                    lab=line[16]
                    res=line[17:20].strip()
                    chain=line[21]
                    resnum=line[22:26]
                    code=line[26]
                    pos=[float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                    occ=line[54:60]
                    temp=line[60:66]
                    amark=line[76:78]
                    charge=line[78:80]

                if num < 100000:
                    head=line[0:6]
                    num=line[6:11]
                    atypename=line[12:16]
                    lab=line[16]
                    res=line[17:20].strip()
                    chain=line[21]
                    resnum=line[22:26]
                    code=line[26]
                    pos=[float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                    occ=line[54:60]
                    temp=line[60:66]
                    amark=line[76:78]
                    charge=line[78:80]

                else:
                    head=line[0:6]
                    num=line[6:12]
                    atypename=line[13:17]
                    lab=line[17]
                    res=line[18:21].strip()
                    chain=line[22]
                    resnum=line[23:27]
                    code=line[27]
                    pos=[float(line[31:39].strip()), float(line[39:47].strip()), float(line[47:55].strip())]
                    occ=line[55:61]
                    temp=line[61:67]
                    amark=line[77:79]
                    charge=line[79:81]

                nums.append(num)
                heads.append(head)
                molnames.append(res)
                atypenames.append(atypename)
                labs.append(lab)
                chains.append(chain)
                resnums.append(resnum)
                codes.append(code)
                poss.append(pos)
                occs.append(occ)
                temps.append(temp)
                amarks.append(amark)
                charges.append(charge)

            totalatom = atomcount
        # print(poss)
        if self.getmode == 'rfile':
            molname_set = set(molnames)
            print('totalatom', totalatom)
            print('molname_set', molname_set)
            mol_confs = []
            for molname in molname_set:
                mol_confs.append(self.config_read(molname, 0))
            # print('mol_confs', mol_confs)
            # print('molnames_atom', molnames)
            molnums = {}
            for molconf in mol_confs:
                molnums[molconf['name']] = sum(molconf['atom'])
            print('molnatoms', molnums)

            #search molhead
            i = 0
            resnames = []
            while True:
                resnames.append(molnames[i])
                atomnum = molnums[molnames[i]]
                i += atomnum
                if i + 1> totalatom:
                    break
                # print(i + 1, molnames[i])
            print('molnames_permol', resnames)
            totalRes = len(resnames)

            gatmlabRes = self.getpermol(totalRes, molnums, resnames, nums)
            resnameRes = self.getpermol(totalRes, molnums, resnames, molnames)
            posRes = self.getpermol(totalRes, molnums, resnames, poss)
            atmtypeRes = self.getpermol(totalRes, molnums, resnames, atypenames)
            headRes = self.getpermol(totalRes, molnums, resnames, heads)
            labRes = self.getpermol(totalRes, molnums, resnames, labs)
            chainRes  = self.getpermol(totalRes, molnums, resnames, chains)
            resnumRes  = self.getpermol(totalRes, molnums, resnames, resnums)
            codeRes  = self.getpermol(totalRes, molnums, resnames, codes)
            occRes  = self.getpermol(totalRes, molnums, resnames, occs)
            tempRes  = self.getpermol(totalRes, molnums, resnames, temps)
            amarkRes  = self.getpermol(totalRes, molnums, resnames, amarks)
            chargeRes = self.getpermol(totalRes, molnums, resnames, charges)

        if self.getmode == 'resnum':
            #search molhead
            i = 0
            atomnum = 0
            resnames = [molnames[0]]
            anummols = []
            # print(resnums)
            while True:
                # print(resnums[i])
                atomnum += 1
                if i + 1 == totalatom:
                    anummols.append(atomnum)
                    break
                if resnums[i] != resnums[i+1]:
                    # print(resnums[i], resnums[i+1])
                    # print(i)
                    anummols.append(atomnum)
                    resnames.append(molnames[i+1])
                    atomnum = 0
                i += 1
                # print(i + 1, molnames[i])

            # print(anummols)
            totalRes = len(anummols)
            # print(totalRes)
            gatmlabRes = self.getpermol2(totalRes, anummols, nums)
            resnameRes = self.getpermol2(totalRes, anummols, molnames)
            posRes = self.getpermol2(totalRes, anummols, poss)
            atmtypeRes = self.getpermol2(totalRes, anummols, atypenames)
            headRes = self.getpermol2(totalRes, anummols, heads)
            labRes = self.getpermol2(totalRes, anummols, labs)
            chainRes  = self.getpermol2(totalRes, anummols, chains)
            resnumRes  = self.getpermol2(totalRes, anummols, resnums)
            codeRes  = self.getpermol2(totalRes, anummols, codes)
            occRes  = self.getpermol2(totalRes, anummols, occs)
            tempRes  = self.getpermol2(totalRes, anummols, temps)
            amarkRes  = self.getpermol2(totalRes, anummols, amarks)
            chargeRes = self.getpermol2(totalRes, anummols, charges)

            if self.assignresname == True:
                moldatas = []
                moldatas_el = []
                molidnames = []
                current = 1
                # loop for allmol
                for i in range(totalRes):
                    if anummols[i] not in moldatas:
                        # print(moldatas, anummols[i])
                        moldatas.append(anummols[i])
                        moldatas_el.append(atmtypeRes[i])
                        # print(anummols[i], 'append')
                        molidnames.append('{:0>3}'.format(str(current)))
                        current += 1
                        continue
                    flags = []
                    # print(moldatas)
                    # loop for databasemol
                    for k in range(len(moldatas)):
                        flag = False
                        if moldatas[k] == anummols[i]:
                            # loop for batabase-1molatom
                            for j in range(len(moldatas_el[k])):
                                if moldatas_el[k][j] != atmtypeRes[i][j]:
                                    flag = True
                                    flags.append(flag)
                                    break
                            flags.append(flag)
                            # print(flags)

                            if flag == False:
                                molidnames.append('{:0>3}'.format(str(k + 1)))
                                break

                    if True in flags and False not in flags:
                        # print('mol', i, 'different!')
                        moldatas.append(anummols[i])
                        moldatas_el.append(atmtypeRes[i])
                        molidnames.append('{:0>3}'.format(str(current)))
                        current += 1

                # print(molidnames)
                resnames = copy.deepcopy(molidnames)
                # print(moldatas)

        self.totalRes = totalRes
        self.atmtypeRes = atmtypeRes
        self.resnames = resnames
        self.resnameRes = resnameRes
        self.gatmlabRes = gatmlabRes
        self.posRes = posRes
        self.headRes = headRes
        self.labRes = labRes
        self.chainRes = chainRes
        self.resnumRes = resnumRes
        self.codeRes = codeRes
        self.occRes = occRes
        self.tempRes = tempRes
        self.amarkRes = amarkRes
        self.chargeRes = chargeRes
        self.rescount = collections.Counter(self.resnames)

        return
        # return totalRes, atmtypeRes, resnames, gatmlabRes, posRes, headRes, labRes, chainRes ,resnumRes ,codeRes ,occRes ,tempRes ,amarkRes ,chargeRes


    def getpermol(self, totalRes, molnums, resnames, datas):
        datamols = []
        count = 0
        for i in range(totalRes):
            datamol = []
            for j in range(molnums[resnames[i]]):
                datamol.append(datas[count])

                count += 1
            datamols.append(datamol)
        return datamols

    def getpermol2(self, totalRes, anummols, datas):
        datamols = []
        # print(anummols)
        count = 0
        for i in range(totalRes):
            datamol = []
            for j in range(anummols[i]):
                datamol.append(datas[count])
                count += 1
            datamols.append(datamol)
            # print(datamols)
        return datamols


#        aobj.exportardpdbfull(opath + '/' + self.readgeom, index, posMol, atomnameMol, self.resnames, heads, labs, chains, resnums, codes, occs, temps, amarks, charges)

    def exportardpdbfull(self, out_file, mollist):

        out_file = out_file + '.pdb'
        print('outfile', out_file)

        # # Export position of mol
        # head, ext = os.path.splitext(str(iname))
        # molid = sorted(set(molnames), key=molnames.index)
        # print(molid)

#         molids = []
#         for molname in molnames:
#             for i in range(len(molid)):
#                 if molname == molid[i]:
#                     molids.append(i)
        # print(molnames)

        # header
        # print(out_file)
        f = open(out_file, "w", newline = "\n")
        print("COMPND    " + out_file, file=f)
        print("AUTHOR    " + "GENERATED BY python script in abmptools", file=f)
        if self.cellsize != 0:
            print('CRYST1{0[0]:>9.3f}{0[1]:>9.3f}{0[2]:>9.3f}  90.00  90.00  90.00               1'.format(self.cellsize), file=f)
            # CRYST1   38.376   38.376   38.376  90.00  90.00  90.00               1
        f.close()

        f = open(out_file, "a+", newline = "\n")
        tatomlab = 0
        reslab = 0
        # print(mollist)


        for i in mollist:
            posMol = self.posRes[i]
            reslab += 1
            if reslab >=10000:
                reslab -= 10000

            for j in range(len(posMol)):
                tatomlab += 1
                if tatomlab >= 100000:
                    tatomlab -= 100000

                if self.refreshatmtype == False:
                    atomname = self.atmtypeRes[i][j]
                    form ='{0[0]:<6}{0[1]:>5} {0[2]:>4}{0[3]:>1}{0[4]:>3} {0[5]:>1}{0[6]:>4}{0[7]:>1}   {0[8]:>8}{0[9]:>8}{0[10]:>8}{0[11]:>6}{0[12]:>6}          {0[13]:>2}{0[14]:>2}'
                else:
                    atomname = self.amarkRes[i][j]
                    form ='{0[0]:<6}{0[1]:>5} {0[2]:>2}  {0[3]:>1}{0[4]:>3} {0[5]:>1}{0[6]:>4}{0[7]:>1}   {0[8]:>8}{0[9]:>8}{0[10]:>8}{0[11]:>6}{0[12]:>6}          {0[13]:>2}{0[14]:>2}'

                if self.refreshresid == True:
                    resid = reslab
                else:
                    resid = self.resnumRes[i][j]

                olist = [self.headRes[i][j], str(tatomlab), atomname, self.labRes[i][j], self.resnames[i], self.chainRes[i][j], resid, self.codeRes[i][j], '{:.3f}'.format(posMol[j][0]), '{:.3f}'.format(posMol[j][1]), '{:.3f}'.format(posMol[j][2]), self.occRes[i][j], self.tempRes[i][j], self.amarkRes[i][j], self.chargeRes[i][j]]
                print(form.format(olist), file=f)
            # ATOM      1  H   UNK     1     -12.899  32.293   3.964  1.00  0.00           H

        print("END", file=f)
        f.close()

    #1 1 – 6  HETATM
    #2 7 – 11 原子の通し番号
    # blank 1
    #3 13 – 16    原子名
    #4 17  Alternate location識別子
    #5 18 - 20 残基名
    # blank 1
    #6 22  鎖名
    #7 23 - 26 残基番号
    #8 27  残基の挿入コード
    # blank 3
    #9 31 - 38 原子のX座標の値（Å単位）
    #10 39 - 46 原子のY座標の値（Å単位）
    #11 47 - 54 原子のZ座標の値（Å単位）
    #12 55 - 60 占有率
    #13 61 - 66 温度因子
    #14 77 - 78 元素記号
    #15 79 - 80 原子の電荷

    def exportardpdb(self, out_file, mollist, posRes, nameAtom, molnames_orig):

        ohead, ext = os.path.splitext(out_file)
        out_file = ohead + '.pdb'

        # # Export position of mol
        # head, ext = os.path.splitext(str(iname))

        molids = sorted(set(molnames_orig), key=molnames_orig.index)
        print(molids)

        molnames = []
        for molname in molnames_orig:
            for i in range(len(molids)):
                if molname == molids[i]:
                    molnames.append(i)

        # print(molnames)
        # header
        # print(out_file)
        f = open(out_file, "w", newline = "\n")
        print("COMPND    " + out_file, file=f)
        print("AUTHOR    " + "GENERATED BY python script in abmptools", file=f)
        f.close()

        f = open(out_file, "a+", newline = "\n")
        tatomlab = 0
        reslab = 0
        print(mollist)
        for i in mollist:
            molname = str(molnames[i])
            Molnum = i
            posMol = posRes[i]
            reslab += 1
            if reslab >=10000:
                reslab -= 10000

            for j in range(len(posMol)):
                tatomlab += 1
                if tatomlab >= 100000:
                    tatomlab -= 100000

                list = ["HETATM", str(tatomlab), nameAtom[i][j], molname.zfill(3), str(i + 1), '{:.3f}'.format(posMol[j][0]), '{:.3f}'.format(posMol[j][1]), '{:.3f}'.format(posMol[j][2]), "1.00", "0.00", nameAtom[i][j]]
                print('{0[0]:<6}{0[1]:>5} {0[2]:>2}   {0[3]:>3}  {0[4]:>4}    {0[5]:>8}{0[6]:>8}{0[7]:>8}{0[8]:>6}{0[9]:>6}{0[10]:>12}'.format(list), file=f)
        # ATOM      1  H   UNK     1     -12.899  32.293   3.964  1.00  0.00           H

        print("END", file=f)
        f.close()

    def getcontact_rmapfmopdb(self, path, fname, oname):
        molname = self.molname
        criteria = self.criteria
        tgtpos = self.tgtpos
        solutes = self.solutes
        self.mainpath = '.'
        self.assignresname = False

        # print(path, molname, fname, oname, tgtpos, criteria)
        # (this func)totalRes, self.atmtypeRes, self.resnames, self.posRes, self.headRes, self.labRes, self.chainRes ,self.resnumRes ,self.codeRes ,self.occRes ,self.tempRes ,self.amarkRes ,self.chargeRes = self.getpdbinfo(fname)
        # (class)    totalRes, atmtypeRes, resnames, gatmlabRes, posRes,     headRes,    labRes,    chainRes , resnumRes ,  codeRes ,  occRes ,  tempRes ,  amarkRes ,  chargeRes

        self.readpdb(fname)

        posRes_orig     = copy.deepcopy(self.posRes)
        atmtypeRes_orig = copy.deepcopy(self.atmtypeRes)
        resnames_orig   = copy.deepcopy(self.resnames)
        headRes_orig    = copy.deepcopy(self.headRes)
        labRes_orig     = copy.deepcopy(self.labRes)
        chainRes_orig   = copy.deepcopy(self.chainRes)
        resnumRes_orig  = copy.deepcopy(self.resnumRes)
        codeRes_orig    = copy.deepcopy( self.codeRes)
        occRes_orig     = copy.deepcopy(self.occRes)
        tempRes_orig    = copy.deepcopy(self.tempRes)
        amarkRes_orig   = copy.deepcopy(self.amarkRes)
        chargeRes_orig  = copy.deepcopy(self.chargeRes)


        print("totalmol:",self.totalRes)
        # getmolatomnum(_udf_, totalRes)

        # 1. -- 切り取らない場合 --
        if self.cutmode == 'none':
            print('none mode')
            # -- get neighbor mol --
            neighborindex = []
            for i in range(self.totalRes):
                neighborindex.append(i)

        # 2. -- 範囲内に原子が入っているかで絞る方法 --
        if self.cutmode == 'sphere':
            print('sphere mode')
            # -- get neighbor mol --
            neighborindex = []
            for i in range(len(self.posRes)):
                # check ion
                icflag = False
                if len(self.ionname) != 0:
                    for ion in self.ionname:
                        if  self.resnames[i].strip() == ion.strip():
                            if self.ionmode == 'del':
                                icflag = True
                                break
                            if self.ionmode == 'remain':
                                icflag = True
                                neighborindex.append(i)
                                # print('add mol', i, self.resnames[i])
                                break

                if icflag == True:
                    icflag == False
                    continue

                for j in range(len(self.posRes[i])):
                    dist = self.getdist(np.array(tgtpos), np.array(self.posRes[i][j]))
                    if dist < criteria:
                        neighborindex.append(i)
                        break

        # 3. -- 球状でなくて正方形で絞る方法
        if self.cutmode == 'cube':
            print('cube mode')
            # -- get neighbor mol --
            tgtx = [tgtpos[0] - criteria[0], tgtpos[0] + criteria[0]]
            tgty = [tgtpos[1] - criteria[1], tgtpos[1] + criteria[1]]
            tgtz = [tgtpos[2] - criteria[2], tgtpos[2] + criteria[2]]

            neighborindex = []
            for i in range(len(self.posRes)):
                # check ion
                icflag = False
                if len(self.ionname) != 0:
                    for ion in self.ionname:
                        if  self.resnames[i].strip() == ion.strip():
                            if self.ionmode == 'del':
                                icflag = True
                                break
                            if self.ionmode == 'remain':
                                icflag = True
                                neighborindex.append(i)
                                # print('add mol', i, self.resnames[i])
                                break

                if icflag == True:
                    icflag == False
                    continue

                for j in range(len(self.posRes[i])):
                    if tgtx[0] < self.posRes[i][j][0] < tgtx[1] and \
                       tgty[0] < self.posRes[i][j][1] < tgty[1] and \
                       tgtz[0] < self.posRes[i][j][2] < tgtz[1]:
                        neighborindex.append(i)
                        break

        # 4. -- 溶質からの距離でスクリーニング --
        if self.cutmode == 'around':
            print('around mode')
            # -- get neighbor mol --
            neighborindex = []
            # solutes = [0, 1]
            for i in self.solutes:
                neighborindex.append(i)
                # print('add mol', i, self.resnames[i])
            for i in range(len(self.posRes)):
                # check ion
                icflag = False
                if len(self.ionname) != 0:
                    for ion in self.ionname:
                        if  self.resnames[i].strip() == ion.strip():
                            if self.ionmode == 'del':
                                icflag = True
                                break
                            if self.ionmode == 'remain':
                                icflag = True
                                neighborindex.append(i)
                                # print('add mol', i, self.resnames[i])
                                break

                if icflag == True:
                    icflag == False
                    continue

                if i % 100 == 0:
                    print('check mol', i)
                nextf = False
                if i in solutes:
                   continue
                for j in range(len(self.posRes[i])):
                    if nextf == True:
                        break
                    for k in self.solutes:
                        if nextf == True:
                            break
                        for l in range(len(self.posRes[k])):
                            dist = self.getdist(np.array(self.posRes[k][l]), np.array(self.posRes[i][j]))
                            if dist < criteria:
                                neighborindex.append(i)
                                # print('add mol', i, self.resnames[i])
                                nextf = True
                                break


        print('neighborindex', neighborindex)
        # print vec
        # print("cellsize", cell)
        posRes = []
        atmtypeRes = []
        resnames = []
        headRes = []
        labRes = []
        chainRes = []
        resnumRes = []
        codeRes = []
        occRes = []
        tempRes = []
        amarkRes = []
        chargeRes = []
        # for i in range(totalRes):
        for i in neighborindex:
        # for i in range(1, 2):
            posRes.append(self.posRes[i])
            atmtypeRes.append(self.atmtypeRes[i])
            resnames.append(self.resnames[i])
            headRes.append(self.headRes[i])
            labRes.append(self.labRes[i])
            chainRes.append(self.chainRes[i])
            resnumRes.append(self.resnumRes[i])
            codeRes.append(self.codeRes[i])
            occRes.append(self.occRes[i])
            tempRes.append(self.tempRes[i])
            amarkRes.append(self.amarkRes[i])
            chargeRes.append(self.chargeRes[i])

        self.posRes = posRes
        self.atmtypeRes = atmtypeRes
        self.resnames = resnames
        self.headRes = headRes
        self.labRes = labRes
        self.chainRes = chainRes
        self.resnumRes = resnumRes
        self.codeRes = codeRes
        self.occRes = occRes
        self.tempRes = tempRes
        self.amarkRes = amarkRes
        self.chargeRes = chargeRes

        # get atomnum
        atomnums = []
        tgtmolnames = []
        for i in range(len(molname)):
            for j in range(self.totalRes):
            # for j in range(1):
                # print (molname[i], self.resnames[j])
                if molname[i] == resnames_orig[j]:
                    atomnums.append(len(posRes_orig[j]))
                    tgtmolnames.append(resnames_orig[j])
                    break
        print ('atomnums', atomnums)

        # print (posRes)
        opath = 'for_abmp'
        # oname = "mdout"
        if os.path.exists(opath) is False:
            print(opath)
            subprocess.call(["mkdir", opath])

        # refresh 
        index = [i for i in range(len(posRes))]
        self.exportardpdbfull(opath + '/' + oname, index)

        self.make_abinput_rmap(tgtmolnames, resnames, oname, opath, atomnums)
        # monomer structure

    def movemoltranspdb(self, posVec, transVec):
        # Parallel shift
        posVec = posVec + transVec
        return posVec

    def getpdbcell(self, fname):
        f = open(fname, 'r').readlines()
        cell = []
        for line in f:
            items = line.split()
            if line[0:5] == 'CRYST':
                cell = [float(items[1]), float(items[2]), float(items[3])]
                break
        return cell

    def moveintocellpdb(self, posMol, totalRes, cell):
        # # Move into cell
        posintoMol = []
        cell2 = []
        print(cell)
        for i in range(3):
            cell2.append(cell[i]/2.0)
        print('cell', cell)
        for i in range(totalRes):
            transVec = np.array([0., 0., 0.])  # x,y,z
            centerOfMol = self.getCenter(posMol[i])
            for j in range(3):
                if centerOfMol[j] > cell2[j]:
                    while centerOfMol[j] > cell2[j]:
                        centerOfMol[j] = centerOfMol[j] - cell[j]
                        transVec[j] = transVec[j] - cell[j]
                elif centerOfMol[j] < -cell2[j]:
                    while centerOfMol[j] < -cell2[j]:
                        centerOfMol[j] = centerOfMol[j] + cell[j]
                        transVec[j] = transVec[j] + cell[j]

            posintoMol.append(self.moveMolTrans(posMol[i], transVec))

        print("move_done.")
        return posintoMol

    def getpdbinfowrap(self, fname):
        # get pdbinfo
        self.readpdb(fname)
        print('totalres', self.totalRes)
        # print('atomnameMol', atomnameMol)
        # print('resnames', molnames)
        print('res_count', collections.Counter(self.resnames))
        # print('nummols', nummols)

        cellsize = self.getpdbcell(fname)
        self.cellsize = cellsize
    #
        if len(self.cellsize) == 0:
            self.cellsize = 0
            print('cellinfo: None')
        else:
            print('cellsize:', self.cellsize)
    #
        return self
