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
        # super().__init__(mainpath)
        pass

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

        ajf_parameter[1] = """'pdb/""" + name + """.pdb'"""
        ajf_parameter[2] = """'""" + name + '-' + \
            self.ajf_method + '-' + self.ajf_basis_set[1:-1] + """.cpf'"""
        ajf_body = self.gen_ajf_body(ajf_parameter)
        print(ajf_body, file=ajf_file)

        return

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

        t1 = time.time()
        for i in range(totalMol):
            posMol_orig.append(self.getposmol(uobj, i))
            typenameMol_orig.append(self.getAtomtypename(uobj, i))
        t2 = time.time()
        print("elapsed_time:", t2-t1)

        # print (typenameMol_orig)
        # print (posMol_orig)


        # -- 重心位置が範囲内かで絞る方法 --
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

        # -- 範囲内に原子が入っているかで絞る方法 --

        # -- get neighbor mol --
        neighborindex = []
        for i in range(len(posMol_orig)):
            for j in range(len(posMol_orig[i])):
                dist = self.getdist(tgtpos, posMol_orig[i][j])
                if dist < criteria:
                    neighborindex.append(i)
                    break

        print(neighborindex)
        # print vec
        print("cellsize", cell)
        posMol = []
        typenameMol = []
        molnamelist = []
        # for i in range(totalMol):
        for i in neighborindex:
        # for i in range(1, 2):
            posMol.append(posMol_orig[i])
            typenameMol.append(typenameMol_orig[i])
            molnamelist.append(molnamelist_orig[i])

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
        self.Exportardpos(opath, oname, index, posMol, typenameMol)

        #fmo param
        self.ajf_method = 'MP2'
        self.ajf_basis_set = '6-31Gdag'
        self.solv_flag = False  # True -> in water , False -> in vacuum
        self.verflag = True
        self.memory = 3000
        self.npro = 8
        self.para_job = 1
        # molnamelist: name list for each molecules
        # oname = 'mdout'
        self.make_abinput_rmap(molname, molnamelist, oname, path, atomnums)
        # monomer structure


if __name__ == "__main__":
    # main
    argvs = sys.argv
    fname = str(argvs[1])
    oname = str(argvs[2])
    _udf_ = UDFManager(fname)

    obj = rmap_fmo()
    totalMol, totalRec = obj.gettotalmol_rec(_udf_)

    tgtpos = [20.0, 20.0, 20.0]
    criteria = 50.0
    molname = ['nafion', 'Bulk_water']
    path = ['.', '.']

    obj.getcontact_rmapfmo(
        totalRec-1, _udf_, totalMol, totalMol, path, molname, oname, tgtpos, criteria)
