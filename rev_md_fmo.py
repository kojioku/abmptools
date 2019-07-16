import numpy as np
from UDFManager import *
import sys
import os
import math
import subprocess
import re
import time
import copy
import fcewsmb.udf_io as mbu
import fcews.abinit_io as fab
import fcewsmb.udfcreate as ufc
import rmdpd.udfrm_io as rud
# Matrix operation

class rmap_fmo(fab.abinit_io, ufc.udfcreate, rud.udfrm_io):
    def __init__(self):
        super().__init__()
        self.ajf_method = 'HF'
        self.ajf_basis_set = '6-31Gdag'
        self.cpfflag = True
        self.solv_flag = False  # True -> in water , False -> in vacuum
        self.verflag = True
        self.memory = 3000
        self.npro = 8
        self.para_job = 1
        self.cutmode = 'sphere'
        self.piedaflag = True
        pass

    def setrfmoparam(self, param_rfmo):
        try:
            self.cutmode = param_rfmo['cutmode']
        except KeyError:
            self.cutmode = 'sphere'

        try:
            self.ajf_method = param_rfmo['ajf_method']
        except KeyError:
            self.ajf_method = 'HF'

        try:
            self.ajf_basis_set = param_rfmo['ajf_basis_set']
        except KeyError:
            self.ajf_basis_set = '6-31Gdag'

        try:
            self.cpfflag = param_rfmo['cpfflag']
        except KeyError:
            self.cpfflag = True

        try:
            self.solv_flag = param_rfmo['solv_flag']
        except KeyError:
            self.solv_flag = False

        try:
            self.piedaflag = param_rfmo['piedaflag']
        except KeyError:
            self.piedaflag = True

        try:
            self.verflag = param_rfmo['verflag']
        except KeyError:
            self.verflag = True

        try:
            self.memory = param_rfmo['memory']
        except KeyError:
            self.memory = 3000

        try:
            self.npro = param_rfmo['npro']
        except KeyError:
            self.npro = 8

    def getendatom(self, connect, term):
        flag = False
        for i in range(len(connect)):
            for j in range(2, len(connect[i])):
                if connect[i][j] == term:
                    return connect[i][1]
        return

    def gettermfrag(self, seg_info, endatom):
        print (seg_info)
        print (endatom)
        for i in range(len(seg_info)):
            for j in range(len(seg_info[i])):
                if seg_info[i][j] - 1 == endatom:
                    print ("hit")
                    return i

    def getpolyconf_rmapfmo(self, seg_conf):
        molname = seg_conf['name']
        seg_atom = seg_conf['atom']
        charge = seg_conf['charge']
        connect_num = seg_conf['connect_num']
        connect = seg_conf['connect']
        seg_info = seg_conf['seg_info']
        repeat = seg_conf['repeat'][0]
        terms = seg_conf['term']

        # get term frag info
        connectinfo = self.getconnectdata(molname + ".pdb")
        termfrags = []
        endatoms = []
        for term in terms:
            endatom = self.getendatom(connectinfo, term - 1)
            termfrag = self.gettermfrag(seg_info, endatom)
            print (termfrag)
            endatoms.append(endatom)
            termfrags.append(termfrag)
        print (termfrags)
        print (charge)
        print (seg_atom)
        print (connect_num)
        print (connect)
        print (terms)
        atom_temp = []
        charge_temp = []
        connect_num_temp = []
        for i in range(repeat):
        # first
            if i == 0:
                for j in range(len(seg_atom)):
                    if j == termfrags[1]:
                        atom_temp.append(seg_atom[j] - 1)
                        charge_temp.append(charge[j] + 1)
                        connect_num_temp.append(connect_num[j])
                    else:
                        atom_temp.append(seg_atom[j])
                        charge_temp.append(charge[j])
                        connect_num_temp.append(connect_num[j])

        # end
            elif i == repeat-1:
                for j in range(len(seg_atom)):
                    if j == termfrags[0]:
                        atom_temp.append(seg_atom[j] - 1)
                        charge_temp.append(charge[j] - 1)
                        connect_num_temp.append(connect_num[j] + 1)
                    else:
                        atom_temp.append(seg_atom[j])
                        charge_temp.append(charge[j])
                        connect_num_temp.append(connect_num[j])

        # middle
            else:
                for j in range(len(seg_atom)):
                    if j == termfrags[0]:
                        atom_temp.append(seg_atom[j] - 1)
                        charge_temp.append(charge[j] - 1)
                        connect_num_temp.append(connect_num[j] + 1)
                    elif j == termfrags[1]:
                        atom_temp.append(seg_atom[j] - 1)
                        charge_temp.append(charge[j] + 1)
                        connect_num_temp.append(connect_num[j])
                    else:
                        atom_temp.append(seg_atom[j])
                        charge_temp.append(charge[j])
                        connect_num_temp.append(connect_num[j])

        #for seginfo
        seg_info_temps = []
        for i in range(repeat):
        # first
            if i == 0:
                for j in range(len(seg_atom)):
                    seg_info_temp = []
                    for k in range(len(seg_info[j])):
                        if seg_info[j][k] == terms[1]:
                            continue
                        elif seg_info[j][k] > terms[1]:
                            seg_info_temp.append(seg_info[j][k] - 1)
                        else:
                            seg_info_temp.append(seg_info[j][k])
                    seg_info_temps.append(seg_info_temp)
                current_num = sum([len(segi) for segi in seg_info_temps])

        # end
            elif i == repeat-1:
                for j in range(len(seg_atom)):
                    seg_info_temp = []
                    for k in range(len(seg_info[j])):
                        if seg_info[j][k] == terms[0]:
                            continue
                        elif seg_info[j][k] > terms[0]:
                            seg_info_temp.append(current_num + seg_info[j][k] - 1)
                        else:
                            seg_info_temp.append(current_num + seg_info[j][k])

                    seg_info_temps.append(seg_info_temp)

        # middle
            else:
                for j in range(len(seg_atom)):
                    seg_info_temp = []
                    for k in range(len(seg_info[j])):
                        if seg_info[j][k] == terms[0] or seg_info[j][k] == terms[1]:
                            continue
                        elif seg_info[j][k] > terms[0] and seg_info[j][k] > terms[1] :
                            seg_info_temp.append(current_num + seg_info[j][k] - 2)
                        elif seg_info[j][k] > terms[0] or seg_info[j][k] > terms[1] :
                            seg_info_temp.append(current_num + seg_info[j][k] - 1)
                        else:
                            seg_info_temp.append(current_num + seg_info[j][k])

                    seg_info_temps.append(seg_info_temp)
                current_num = sum([len(segi) for segi in seg_info_temps])

        #for connect
        segatomnum = sum([len(segi) for segi in seg_info])
        print (endatoms)
        connect_temps = []
        for i in range(repeat):
        # first
            if i == 0:
                for j in range(len(connect)):
                    connect_temp = []
                    for k in range(len(connect[j])):
                        if connect[j][k] > terms[1]:
                            connect_temp.append(connect[j][k] - 1)
                        else:
                            connect_temp.append(connect[j][k])
                    connect_temps.append(connect_temp)
                new_conn = [endatoms[1] + 1, endatoms[0] + segatomnum]
                connect_temps.append(new_conn)
                current_num = segatomnum -1

        # end
            elif i == repeat-1:
                for j in range(len(connect)):
                    connect_temp = []
                    for k in range(len(connect[j])):
                        if connect[j][k] > terms[0]:
                            connect_temp.append(current_num + connect[j][k] - 1)
                        else:
                            connect_temp.append(current_num + connect[j][k])

                    connect_temps.append(connect_temp)

        # middle
            else:
                for j in range(len(connect)):
                    connect_temp = []
                    for k in range(len(connect[j])):
                        if connect[j][k] > terms[0] and connect[j][k] > terms[1] :
                            connect_temp.append(current_num + connect[j][k] - 2)
                        elif connect[j][k] > terms[0] or connect[j][k] > terms[1] :
                            connect_temp.append(current_num + connect[j][k] - 1)
                        else:
                            connect_temp.append(current_num + connect[j][k])
                    connect_temps.append(connect_temp)

                new_conn = [aaa + segatomnum - 2 for aaa in new_conn]
                connect_temps.append(new_conn)
                current_num = current_num + segatomnum - 2

        print (connect_temps)
        print (seg_info_temps)
        print ('atom_temp', atom_temp)
        print ('charge_temp', charge_temp)
        print ('connect_temp', connect_num_temp)


        # poly_atom =
        poly_conf = {
                'name': molname,
                'atom': atom_temp,
                'charge': charge_temp,
                'connect_num': connect_num_temp,
                'connect': connect_temps,
                'seg_info': seg_info_temps,
                'nummol_seg' : [1],
                'repeat': [repeat]
                    }


        return poly_conf

    def make_abinput_rmap(self, fname, molnamelist, rec, path, atomnums):
        # print make_input_param

        # fragment configure reading
        frag_atom = []
        frag_charge = []
        frag_connect_num = []
        frag_connect = []
        seg_info = []
        nummol_seg = []

        for i in range(len(fname)):
            mol_conf = self.config_read(fname[i], atomnums[i])

            if mol_conf['repeat'][0] != 1:
                mol_conf = self.getpolyconf_rmapfmo(mol_conf)

            frag_atom.append(mol_conf['atom'])
            frag_charge.append(mol_conf['charge'])
            frag_connect_num.append(mol_conf['connect_num'])
            frag_connect.append(mol_conf['connect'])
            seg_info.append(mol_conf['seg_info'])
            nummol_seg.append(mol_conf['nummol_seg'])

        # frag_atom = [mol1_conf['atom'], mol2_conf['atom']]
        # frag_charge = [mol1_conf['charge'], mol2_conf['charge']]
        # frag_connect_num = [mol1_conf['connect_num'], mol2_conf['connect_num']]
        # frag_connect = [mol1_conf['connect'], mol2_conf['connect']]
        # seg_info = [mol1_conf['seg_info'], mol2_conf['seg_info']]
        # nummol_seg = [mol1_conf['nummol_seg'], mol2_conf['nummol_seg']]

        nameid = []
        print("molnamelist", molnamelist)
        for i in range(len(molnamelist)):
            for j in range(len(fname)):
                if molnamelist[i] == fname[j]:
                    nameid.append(j)
        print("molnameid", nameid)
        # fragment body
        ajf_fragment = self.writemb_frag_section(
            [frag_atom, frag_charge, frag_connect_num, frag_connect, seg_info],
            nameid)
        ajf_fragment = ajf_fragment[:-1]

        # ajf_parameter = [charge, method, basis_set, pdb_filename, cpf_filename,
        # fragment section ]
        ajf_charge = 0
        for i in range(len(nameid)):
            ajf_charge += sum(frag_charge[nameid[i]])

        num_fragment = 0
        for i in range(len(nameid)):
            num_fragment += len(frag_atom[nameid[i]])

        print("molnum")
        for i in range(len(nameid)):
            print("[", nameid.count(i), "]")

        print("num_fragment")
        print(num_fragment)

        Si_flag = 0
        ajf_parameter = [ajf_charge, "", "", ajf_fragment, num_fragment, Si_flag]
        # gen ajf file
        # --for one file--
        name = str(rec)
        ajf_file_name = path[0] + "/" + path[1] + "/" + name + ".ajf"
        ajf_file = open(ajf_file_name, 'w')

        ajf_parameter[1] = "'pdb/" + name + ".pdb'"
        ajf_parameter[2] = "'" + name + '-' + \
            self.ajf_method + '-' + self.ajf_basis_set + ".cpf'"
        ajf_body = self.gen_ajf_body(ajf_parameter)
        print(ajf_body, file=ajf_file)

        return

    def getpdbinfo(self, fname):

        print(fname)
        lines = open(fname,'r').readlines()
        molnames = []
        poss = []
        atypenames = []
        atomcount = 0
        for line in lines:
            itemlist = line.split()
            # print(itemlist)
            if len(itemlist) < 3:
                continue
            if itemlist[0] == 'ATOM' or itemlist[0] == 'HETATM':
                molname = itemlist[3]
                molnames.append(molname)
                poss.append([float(e) for e in itemlist[5:8]])
                atypenames.append(itemlist[2])
                atomcount += 1
            # over HETATM10000
            elif len(itemlist[0]) > 6 and itemlist[0][0:6] == 'HETATM':
                molname = itemlist[2]
                molnames.append(molname)
                poss.append([float(e) for e in itemlist[4:7]])
                atypenames.append(itemlist[1])
                atomcount += 1

        # print(poss)
        molname_set = set(molnames)
        totalatom = atomcount
        print('totalatom', totalatom)
        print(molname_set)
        mol_confs = []
        for molname in molname_set:
            mol_confs.append(self.config_read(molname, 0))
        print(mol_confs)
        # print(molnames)
        molnums = {}
        for molconf in mol_confs:
            molnums[molconf['name']] = sum(molconf['atom'])
        print(molnums)

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

        posMols = []
        count = 0
        print('totalposs', len(poss))
        for i in range(totalMol):
            posMol = []
            for j in range(molnums[name_permol[i]]):
                posMol.append(poss[count])
                count += 1
            posMols.append(posMol)
        # print(posMols)

        typenameMols = []
        count = 0
        for i in range(totalMol):
            typenameMol = []
            for j in range(molnums[name_permol[i]]):
                typenameMol.append(atypenames[count])
                count += 1
            typenameMols.append(typenameMol)
        # print(typenameMols)

        return totalMol, typenameMols, name_permol, posMols


    def getcontact_rmapfmopdb(self, path,  molname, fname, oname, tgtpos, criteria):
        self.mainpath = '.'

        # print(path, molname, fname, oname, tgtpos, criteria)
        totalMol, atomnameMol_orig, molnamelist_orig, posMol_orig = self.getpdbinfo(fname)

        print("totalmol:",totalMol)
        # getmolatomnum(_udf_, totalMol)


        # 2. -- 範囲内に原子が入っているかで絞る方法 --
        if self.cutmode == 'sphere':
            print('sphere mode')
            # -- get neighbor mol --
            neighborindex = []
            for i in range(len(posMol_orig)):
                for j in range(len(posMol_orig[i])):
                    dist = self.getdist(tgtpos, posMol_orig[i][j])
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
                for j in range(len(posMol_orig[i])):
                    if tgtx[0] < posMol_orig[i][j][0] < tgtx[1] and \
                       tgty[0] < posMol_orig[i][j][1] < tgty[1] and \
                       tgtz[0] < posMol_orig[i][j][2] < tgtz[1]:
                        neighborindex.append(i)
                        break

        print('neighborindex', neighborindex)
        # print vec
        # print("cellsize", cell)
        posMol = []
        atomnameMol = []
        molnamelist = []
        # for i in range(totalMol):
        for i in neighborindex:
        # for i in range(1, 2):
            posMol.append(posMol_orig[i])
            atomnameMol.append(atomnameMol_orig[i])
            molnamelist.append(molnamelist_orig[i])

        # get atomnum
        atomnums = []
        for i in range(len(molname)):
            for j in range(totalMol):
            # for j in range(1):
                # print (molname[i], molnamelist_orig[j])
                if molname[i] == molnamelist_orig[j]:
                    atomnums.append(len(posMol_orig[j]))
                    break
        print ('atomnums', atomnums)

        # print (posMol)
        opath = 'pdb'
        # oname = "mdout"
        if os.path.exists(opath) is False:
            print(opath)
            subprocess.call(["mkdir", opath])

        index = [i for i in range(len(posMol))]
        self.exportardpdb(opath + '/' + oname, index, posMol, atomnameMol, molnamelist)

        self.make_abinput_rmap(molname, molnamelist, oname, path, atomnums)
        # monomer structure


    def getcontact_rmapfmo(self, rec, uobj, totalMol, inmol, path,  molname, oname, tgtpos, criteria):

        uobj.jump(rec)
        self.moveintocell_rec(_udf_, rec, totalMol)

        cell = uobj.get("Structure.Unit_Cell.Cell_Size")
        print("totalmol:",totalMol, "rec", rec)
        # getmolatomnum(_udf_, totalMol)

        index = [ i for i in range(totalMol)]
        molnamelist_orig = self.getnamelist(index, uobj, totalMol)
        # print (molnamelist)

        # start
        posMol_orig = []
        typenameMol_orig = []
        atomname_orig = []

        t1 = time.time()
        for i in range(totalMol):
            posMol_orig.append(self.getposmol(uobj, i))
            typenameMol_orig.append(self.getAtomtypename(uobj, i))
            atomname_orig.append(self.getnameAtom(uobj, i))
        t2 = time.time()
        print("elapsed_time:", t2-t1)

        # print (typenameMol_orig)
        # print (posMol_orig)


        # 1. -- 重心位置が範囲内かで絞る方法 --
        # centerOfMol = []
        # for i in range(len(posMol_orig)):
        #     centerOfMol.append(getCenter(posMol_orig[i]))

        # # exportxyz(path[0] + "/" + path[1], centerOfMol, "com")

        # # -- get com dist --
        # distlist = []
        # for i in range(len(posMol_orig)):
        #     distlist.append(getdist(tgtpos, centerOfMol[i]))

        # # -- get neighbor mol --
        # neighborindex = []
        # for i in range(len(posMol_orig)):
        #     if distlist[i] < criteria:
        #         neighborindex.append(i)

        # 2. -- 範囲内に原子が入っているかで絞る方法 --
        if self.cutmode == 'sphere':
            print('sphere mode')
            # -- get neighbor mol --
            neighborindex = []
            for i in range(len(posMol_orig)):
                for j in range(len(posMol_orig[i])):
                    dist = self.getdist(tgtpos, posMol_orig[i][j])
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
                for j in range(len(posMol_orig[i])):
                    if tgtx[0] < posMol_orig[i][j][0] < tgtx[1] and \
                       tgty[0] < posMol_orig[i][j][1] < tgty[1] and \
                       tgtz[0] < posMol_orig[i][j][2] < tgtz[1]:
                        neighborindex.append(i)
                        break

        print(neighborindex)
        # print vec
        print("cellsize", cell)
        posMol = []
        typenameMol = []
        molnamelist = []
        atomnameMol = []
        # for i in range(totalMol):
        for i in neighborindex:
        # for i in range(1, 2):
            posMol.append(posMol_orig[i])
            typenameMol.append(typenameMol_orig[i])
            molnamelist.append(molnamelist_orig[i])
            atomnameMol.append(atomname_orig[i])

        # get atomnum
        atomnums = []
        for i in range(len(molname)):
            for j in range(totalMol):
            # for j in range(1):
                print (molname[i], molnamelist_orig[j])
                if molname[i] == molnamelist_orig[j]:
                    atomnums.append(len(posMol_orig[j]))
                    break
        print ('atomnums', atomnums)

        # print (posMol)
        opath = 'pdb'
        # oname = "mdout"
        index = [i for i in range(len(posMol))]
        self.Exportardpos(opath, oname, index, posMol, atomnameMol)
        self.exportardpdb(opath + '/' + oname, index, posMol, atomnameMol, molnamelist)

        #fmo param
        # self.ajf_method = 'HF'
        # self.ajf_basis_set = '6-31Gdag'
        # self.cpfflag = True
        # self.solv_flag = False  # True -> in water , False -> in vacuum
        # self.verflag = True
        # self.memory = 3000
        # self.npro = 8
        # self.para_job = 1
        # molnamelist: name list for each molecules
        # oname = 'mdout'
        self.make_abinput_rmap(molname, molnamelist, oname, path, atomnums)
        # monomer structure


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


        # aaa = [0.8855]
        # print '{0[0]:.3f}'.format(aaa)
        # print pos[0]
        f = open(out_file, "a+", newline = "\n")
        tatomlab = 0
        print('aaaaaaaaa', mollist)
        for i in mollist:
            molname = str(molnames[i])
            Molnum = i
            posMol = posMols[i]
            for j in range(len(posMol)):
                tatomlab += 1
                list = ["HETATM", str(tatomlab), nameAtom[i][j], molname.zfill(3), str(i), '{:.3f}'.format(posMol[j][0]), '{:.3f}'.format(posMol[j][1]), '{:.3f}'.format(posMol[j][2]), "1.00", "0.00", nameAtom[i][j]]
                print('{0[0]:<6}{0[1]:>5}{0[2]:>4}{0[3]:>5}{0[4]:>6}{0[5]:>12}{0[6]:>8}{0[7]:>8}{0[8]:>6}{0[9]:>6}{0[10]:>12}'.format(list), file=f)
        # ATOM      1  H   UNK     1     -12.899  32.293   3.964  1.00  0.00           H

        print("END", file=f)
        f.close()


if __name__ == "__main__":
    # main
    argvs = sys.argv
    fname = str(argvs[1])
    oname = str(argvs[2])
    _udf_ = UDFManager(fname)

    obj = rmap_fmo()
    totalMol, totalRec = obj.gettotalmol_rec(_udf_)
    path = ['.', '.']

    param_read = {}
    exec(open("input_param", 'r').read(), param_read)
    param_rfmo = param_read['param']
    obj.setrfmoparam(param_rfmo)

    tgtpos = param_rfmo['tgtpos']
    criteria = param_rfmo['criteria']
    molname = param_rfmo['molname']

    obj.getcontact_rmapfmo(
        totalRec-1, _udf_, totalMol, totalMol, path, molname, oname, tgtpos, criteria)
