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
import fcews.abinit_io as fab
import fcewsmb.udfcreate as ufc
import rmdpd.udfrm_io as rud
try:
    from UDFManager import *
except:
    pass

class pdb_io(fab.abinit_io):
    def __init__(self):
        super().__init__()
        self.solutes = []
        self.getmode = 'resnum'
        self.assignmolname = True
        self.refreshatmtype = False
        self.cellsize = 0
        pass

    def getpdbinfo(self, fname):

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

        posMols = []
        typenameMols = []
        headMols = []
        labMols = []
        chainMols = []
        resnumMols = []
        codeMols = []
        occMols = []
        tempMols = []
        amarkMols = []
        chargeMols = []
        totalMol = 0
        name_permol = []

        atomcount = 0
        for line in lines:
            itemlist = line.split()
            # print(itemlist)
            if self.getmode == 'TER':
                if len(itemlist) >= 1 and itemlist[0] == 'TER':
                    print(itemlist)
                    posMols.append(poss)
                    typenameMols.append(atypenames)
                    headMols.append(heads)
                    labMols.append(labs)
                    chainMols.append(chains)
                    resnumMols.append(resnums)
                    codeMols.append(codes)
                    occMols.append(occs)
                    tempMols.append(temps)
                    amarkMols.append(amarks)
                    chargeMols.append(charges)
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
                    name_permol.append(totalMol)
                    totalMol += 1
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

                if atomcount < 100000:
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
            print('mol_confs', mol_confs)
            print('molnames_atom', molnames)
            molnums = {}
            for molconf in mol_confs:
                molnums[molconf['name']] = sum(molconf['atom'])
            print('molnums', molnums)

            #search molhead
            i = 0
            name_permol = []
            while True:
                name_permol.append(molnames[i])
                atomnum = molnums[molnames[i]]
                i += atomnum
                if i + 1> totalatom:
                    break
                # print(i + 1, molnames[i])
            print('molnames_permol', name_permol)
            totalMol = len(name_permol)

            posMols = self.getpermol(totalMol, molnums, name_permol, poss)
            typenameMols = self.getpermol(totalMol, molnums, name_permol, atypenames)
            headMols = self.getpermol(totalMol, molnums, name_permol, heads)
            labMols = self.getpermol(totalMol, molnums, name_permol, labs)
            chainMols  = self.getpermol(totalMol, molnums, name_permol, chains)
            resnumMols  = self.getpermol(totalMol, molnums, name_permol, resnums)
            codeMols  = self.getpermol(totalMol, molnums, name_permol, codes)
            occMols  = self.getpermol(totalMol, molnums, name_permol, occs)
            tempMols  = self.getpermol(totalMol, molnums, name_permol, temps)
            amarkMols  = self.getpermol(totalMol, molnums, name_permol, amarks)
            chargeMols = self.getpermol(totalMol, molnums, name_permol, charges)

        if self.getmode == 'resnum':
            #search molhead
            i = 0
            atomnum = 0
            name_permol = [molnames[0]]
            anummols = []
            # print(resnums)
            while True:
                # print(resnums[i])
                atomnum += 1
                if i + 1 == totalatom:
                    anummols.append(atomnum)
                    break
                if resnums[i] < resnums[i+1]:
                    # print(resnums[i], resnums[i+1])
                    # print(i)
                    anummols.append(atomnum)
                    name_permol.append(molnames[i+1])
                    atomnum = 0
                i += 1
                # print(i + 1, molnames[i])

            # print(anummols)
            totalMol = len(anummols)
            # print(totalMol)
            posMols = self.getpermol2(totalMol, anummols, poss)
            typenameMols = self.getpermol2(totalMol, anummols, atypenames)
            headMols = self.getpermol2(totalMol, anummols, heads)
            labMols = self.getpermol2(totalMol, anummols, labs)
            chainMols  = self.getpermol2(totalMol, anummols, chains)
            resnumMols  = self.getpermol2(totalMol, anummols, resnums)
            codeMols  = self.getpermol2(totalMol, anummols, codes)
            occMols  = self.getpermol2(totalMol, anummols, occs)
            tempMols  = self.getpermol2(totalMol, anummols, temps)
            amarkMols  = self.getpermol2(totalMol, anummols, amarks)
            chargeMols = self.getpermol2(totalMol, anummols, charges)

            if self.assignmolname == True:
                moldatas = []
                moldatas_el = []
                molidnames = []
                current = 1
                # loop for allmol
                for i in range(totalMol):
                    if anummols[i] not in moldatas:
                        # print(moldatas, anummols[i])
                        moldatas.append(anummols[i])
                        moldatas_el.append(typenameMols[i])
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
                                if moldatas_el[k][j] != typenameMols[i][j]:
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
                        moldatas_el.append(typenameMols[i])
                        molidnames.append('{:0>3}'.format(str(current)))
                        current += 1

                # print(molidnames)
                name_permol = copy.deepcopy(molidnames)
                # print(moldatas)

        return totalMol, typenameMols, name_permol, posMols, headMols, labMols, chainMols ,resnumMols ,codeMols ,occMols ,tempMols ,amarkMols ,chargeMols


    def getpermol(self, totalMol, molnums, name_permol, datas):
        datamols = []
        count = 0
        for i in range(totalMol):
            datamol = []
            for j in range(molnums[name_permol[i]]):
                datamol.append(datas[count])

                count += 1
            datamols.append(datamol)
        return datamols

    def getpermol2(self, totalMol, anummols, datas):
        datamols = []
        # print(anummols)
        count = 0
        for i in range(totalMol):
            datamol = []
            for j in range(anummols[i]):
                datamol.append(datas[count])
                count += 1
            datamols.append(datamol)
            # print(datamols)
        return datamols


    def exportardpdbfull(self, out_file, mollist, posMols, nameAtom, molnames, heads, labs, chains, resnums, codes, occs, temps, amarks, charges):

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
        print("AUTHOR    " + "GENERATED BY python script in FMOrmap", file=f)
        if self.cellsize != 0:
            print('CRYST1{0[0]:>9.3f}{0[1]:>9.3f}{0[2]:>9.3f}  90.00  90.00  90.00               1'.format(self.cellsize), file=f)
            # CRYST1   38.376   38.376   38.376  90.00  90.00  90.00               1
        f.close()

        f = open(out_file, "a+", newline = "\n")
        tatomlab = 0
        # print(mollist)

        if self.refreshatmtype == False:

            for i in mollist:
                posMol = posMols[i]
                for j in range(len(posMol)):
                    tatomlab += 1
                    olist = [heads[i][j], str(tatomlab), nameAtom[i][j], labs[i][j], molnames[i], chains[i][j], str(i + 1), codes[i][j], '{:.3f}'.format(posMol[j][0]), '{:.3f}'.format(posMol[j][1]), '{:.3f}'.format(posMol[j][2]), occs[i][j], temps[i][j], amarks[i][j], charges[i][j]]
                    print('{0[0]:<6}{0[1]:>5} {0[2]:>4}{0[3]:>1}{0[4]:>3} {0[5]:>1}{0[6]:>4}{0[7]:>1}   {0[8]:>8}{0[9]:>8}{0[10]:>8}{0[11]:>6}{0[12]:>6}          {0[13]:>2}{0[14]:>2}'.format(olist), file=f)
            # ATOM      1  H   UNK     1     -12.899  32.293   3.964  1.00  0.00           H


        else:

            for i in mollist:
                posMol = posMols[i]
                for j in range(len(posMol)):
                    tatomlab += 1
                    olist = [heads[i][j], str(tatomlab), amarks[i][j], labs[i][j], molnames[i], chains[i][j], str(i + 1), codes[i][j], '{:.3f}'.format(posMol[j][0]), '{:.3f}'.format(posMol[j][1]), '{:.3f}'.format(posMol[j][2]), occs[i][j], temps[i][j], amarks[i][j], charges[i][j]]
                    print('{0[0]:<6}{0[1]:>5} {0[2]:>2}  {0[3]:>1}{0[4]:>3} {0[5]:>1}{0[6]:>4}{0[7]:>1}   {0[8]:>8}{0[9]:>8}{0[10]:>8}{0[11]:>6}{0[12]:>6}          {0[13]:>2}{0[14]:>2}'.format(olist), file=f)
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

    def exportardpdb(self, out_file, mollist, posMols, nameAtom, molnames_orig):

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
        print("AUTHOR    " + "GENERATED BY python script in FMOrmap", file=f)
        f.close()

        f = open(out_file, "a+", newline = "\n")
        tatomlab = 0
        print(mollist)
        for i in mollist:
            molname = str(molnames[i])
            Molnum = i
            posMol = posMols[i]
            for j in range(len(posMol)):
                tatomlab += 1
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
        self.assignmolname = False

        # print(path, molname, fname, oname, tgtpos, criteria)
        totalMol, atomnameMol_orig, molnamelist_orig, posMol_orig, heads_orig, labs_orig, chains_orig ,resnums_orig ,codes_orig ,occs_orig ,temps_orig ,amarks_orig ,charges_orig = self.getpdbinfo(fname)
        # print(molnamelist_orig)

        print("totalmol:",totalMol)
        # getmolatomnum(_udf_, totalMol)

        # 1. -- 切り取らない場合 --
        if self.cutmode == 'none':
            print('none mode')
            # -- get neighbor mol --
            neighborindex = []
            for i in range(totalMol):
                neighborindex.append(i)

        # 2. -- 範囲内に原子が入っているかで絞る方法 --
        if self.cutmode == 'sphere':
            print('sphere mode')
            # -- get neighbor mol --
            neighborindex = []
            for i in range(len(posMol_orig)):
                # check ion
                icflag = False
                if len(self.ionname) != 0:
                    for ion in self.ionname:
                        if  molnamelist_orig[i].strip() == ion.strip():
                            if self.ionmode == 'del':
                                icflag = True
                                break
                            if self.ionmode == 'remain':
                                icflag = True
                                neighborindex.append(i)
                                # print('add mol', i, molnamelist_orig[i])
                                break

                if icflag == True:
                    icflag == False
                    continue

                for j in range(len(posMol_orig[i])):
                    dist = self.getdist(np.array(tgtpos), np.array(posMol_orig[i][j]))
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
            for i in range(len(posMol_orig)):
                # check ion
                icflag = False
                if len(self.ionname) != 0:
                    for ion in self.ionname:
                        if  molnamelist_orig[i].strip() == ion.strip():
                            if self.ionmode == 'del':
                                icflag = True
                                break
                            if self.ionmode == 'remain':
                                icflag = True
                                neighborindex.append(i)
                                # print('add mol', i, molnamelist_orig[i])
                                break

                if icflag == True:
                    icflag == False
                    continue

                for j in range(len(posMol_orig[i])):
                    if tgtx[0] < posMol_orig[i][j][0] < tgtx[1] and \
                       tgty[0] < posMol_orig[i][j][1] < tgty[1] and \
                       tgtz[0] < posMol_orig[i][j][2] < tgtz[1]:
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
                # print('add mol', i, molnamelist_orig[i])
            for i in range(len(posMol_orig)):
                # check ion
                icflag = False
                if len(self.ionname) != 0:
                    for ion in self.ionname:
                        if  molnamelist_orig[i].strip() == ion.strip():
                            if self.ionmode == 'del':
                                icflag = True
                                break
                            if self.ionmode == 'remain':
                                icflag = True
                                neighborindex.append(i)
                                # print('add mol', i, molnamelist_orig[i])
                                break

                if icflag == True:
                    icflag == False
                    continue

                if i % 100 == 0:
                    print('check mol', i)
                nextf = False
                if i in solutes:
                   continue
                for j in range(len(posMol_orig[i])):
                    if nextf == True:
                        break
                    for k in self.solutes:
                        if nextf == True:
                            break
                        for l in range(len(posMol_orig[k])):
                            dist = self.getdist(np.array(posMol_orig[k][l]), np.array(posMol_orig[i][j]))
                            if dist < criteria:
                                neighborindex.append(i)
                                # print('add mol', i, molnamelist_orig[i])
                                nextf = True
                                break


        print('neighborindex', neighborindex)
        # print vec
        # print("cellsize", cell)
        posMol = []
        atomnameMol = []
        molnamelist = []
        heads = []
        labs = []
        chains = []
        resnums = []
        codes = []
        occs = []
        temps = []
        amarks = []
        charges = []
        # for i in range(totalMol):
        for i in neighborindex:
        # for i in range(1, 2):
            posMol.append(posMol_orig[i])
            atomnameMol.append(atomnameMol_orig[i])
            molnamelist.append(molnamelist_orig[i])
            heads.append(heads_orig[i])
            labs.append(labs_orig[i])
            chains.append(chains_orig[i])
            resnums.append(resnums_orig[i])
            codes.append(codes_orig[i])
            occs.append(occs_orig[i])
            temps.append(temps_orig[i])
            amarks.append(amarks_orig[i])
            charges.append(charges_orig[i])

        # get atomnum
        atomnums = []
        tgtmolnames = []
        for i in range(len(molname)):
            for j in range(totalMol):
            # for j in range(1):
                # print (molname[i], molnamelist_orig[j])
                if molname[i] == molnamelist_orig[j]:
                    atomnums.append(len(posMol_orig[j]))
                    tgtmolnames.append(molnamelist_orig[j])
                    break
        print ('atomnums', atomnums)

        # print (posMol)
        opath = 'for_abmp'
        # oname = "mdout"
        if os.path.exists(opath) is False:
            print(opath)
            subprocess.call(["mkdir", opath])

        # if os.path.exists(path[0] + '/' + path[1]) is False:
            # subprocess.call(["mkdir -p", path[0] + '/' + path[1]])

        index = [i for i in range(len(posMol))]
        self.exportardpdbfull(opath + '/' + oname, index, posMol, atomnameMol, molnamelist, heads, labs, chains, resnums, codes, occs, temps, amarks, charges)

        self.make_abinput_rmap(tgtmolnames, molnamelist, oname, path, atomnums)
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

    def moveintocellpdb(self, posMol, totalMol, cell):
        # # Move into cell
        posintoMol = []
        cell2 = []
        print(cell)
        for i in range(3):
            cell2.append(cell[i]/2.0)
        print('cell', cell)
        for i in range(totalMol):
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

