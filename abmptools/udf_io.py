import sys
import os
scrdir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scrdir)
import math
import subprocess
import re
import time
import copy
import molcalc as molc
try:
    from UDFManager import *
except ImportError:
    pass
try:
    import numpy as np
except ImportError:
    pass


class udf_io(molc.molcalc):
    def __init__(self):
        self.verflag = True

    def getposatom(self, uobj, indexatom):
        # UDF operation
        # get the position of a molecule
        i = indexatom
        list = uobj.get("Structure.Position.mol[].atom[]", [i])
        position = np.array(list)

        return position

    def getposmolrec(self, uobj, indexMol, record):
        # UDF operation
        # get the position of a molecule
        uobj.jump(record)
        i = indexMol
        aaa = uobj.get("Structure.Position.mol[" + str(i) + "].atom[]")
        position = np.array(aaa)

        return position

    def getposmol(self, uobj, indexMol):
        # UDF operation
        # get the position of a molecule
        i = indexMol
        list = uobj.get("Structure.Position.mol[" + str(i) + "].atom[]")
        position = np.array(list)

        return position

    def getnameAtom(self, uobj, indexMol):
        # get atom name
        i = indexMol
        atom = uobj.get("Set_of_Molecules.molecule[].atom[].Atom_Name", [i])
        return atom

    def getAtomtypename(self, uobj, indexMol):
        # get atom name
        i = indexMol
        atomtype = uobj.get(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [i])
        return atomtype

    def putPositionsMol(self, uobj, indexMol, position):
        # ## put the position of the molecule to UDF
        i = indexMol
        numAtm = uobj.size("Structure.Position.mol[].atom[]", [i])
        for j in range(numAtm):
            uobj.put(position[j, 0], "Structure.Position.mol[].atom[].x", [i, j])
            uobj.put(position[j, 1], "Structure.Position.mol[].atom[].y", [i, j])
            uobj.put(position[j, 2], "Structure.Position.mol[].atom[].z", [i, j])

    def Exportpos(self, path, Rec, totalMol, uobj, oname):
        # # Export position of mol
        # head, ext = os.path.splitext(str(iname))
        if os.path.exists(path + "/pdb") is False:
            subprocess.call(["mkdir", path + "/pdb"])

        out_file = path + "/pdb/" + str(oname)
        print (out_file)

        numlist = []
        # get total atom
        for i in range(totalMol):
            numAtm = uobj.size("Structure.Position.mol[].atom[]", [i])
            numlist.append(numAtm)
        totalAtm = sum(numlist)
        print(("totalAtom", totalAtm))

        f = open(out_file, "w")
        print(totalAtm, file=f)
        print(out_file, file=f)
        f.close()

        f = open(out_file, "a+")
        uobj.jump(Rec)

        mollist = [i for i in range(totalMol)]
        for i in mollist:
            Molnum = i
            posMol = self.getposmol(uobj, Molnum)
            # print posMol
            nameAtom = self.getnameAtom(uobj, Molnum)
            atom_index = 0
            for j in range(len(posMol)):
                print(nameAtom[atom_index], \
                    posMol[j][0], posMol[j][1], posMol[j][2], file=f)
                atom_index += 1
        f.close()

        # subprocess.call(["babel", "-ixyz", out_file, "-opdb",
        #                  path + "/pdb/mdout_orig.pdb"])
        self.exportpdb(uobj, Rec, out_file, mollist)

    def Exporttgtmolpos(self, path, oname_i, Rec, mollist, uobj):

        if os.path.exists(path) is False:
            print(path)
            subprocess.call(["mkdir", path])

        # print mollist
        # # Export position of mol
        ohead, ext = os.path.splitext(oname_i)
        out_file = path + "/pdb/" + str(ohead) + 'mol' + str(mollist[0]) + '.xyz'
        # print out_file

        numlist = []
        # get total atom
        for i in mollist:
            numAtm = uobj.size("Structure.Position.mol[].atom[]", [i])
            numlist.append(numAtm)
        totalAtm = sum(numlist)
        # print "totalAtom", totalAtm

        f = open(out_file, "w")
        print(totalAtm, file=f)
        print(str(ohead), file=f)
        f.close()

        f = open(out_file, "a+")
        uobj.jump(Rec)

        for i in mollist:
            Molnum = i
            posMol = self.getposmol(uobj, Molnum)
            nameAtom = self.getnameAtom(uobj, Molnum)
            atom_index = 0
            for j in range(len(posMol)):
                print(nameAtom[atom_index], \
                    posMol[j][0], posMol[j][1], posMol[j][2], file=f)
                atom_index += 1
        f.close()

        self.exportpdb(uobj, Rec, out_file, mollist)

        # cmd = "babel -ixyz " +  out_file +  " -opdb ", path + "/" + str(ohead) + ".pdb"
        # subprocess.call(cmd, shell = True)

    def exportpdb(self, uobj, Rec, out_file, mollist):

        ohead, ext = os.path.splitext(out_file)
        out_file = ohead + '.pdb'
        uobj.jump(Rec)

        # # Export position of mol
        # head, ext = os.path.splitext(str(iname))

        molnames_orig = []
        for i in mollist:
            molnames_orig.append(self.getmolname(i, uobj))
        # print(molnames)
        # print(len(molnames_orig))

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
        if self.cell != None:
            print('CRYST1{0[0]:>9.3f}{0[1]:>9.3f}{0[2]:>9.3f}  90.00  90.00  90.00               1'.format(self.cell), file=f)
        f.close()

        # aaa = [0.8855]
        # print '{0[0]:.3f}'.format(aaa)
        # print pos[0]
        f = open(out_file, "a+", newline = "\n")
        tatomlab = 0
        mollab = 0
        for i in mollist:
            molname = str(molnames[mollab])
            mollab += 1
            Molnum = i
            posMol = self.getposmol(uobj, Molnum)
            nameAtom = self.getnameAtom(uobj, Molnum)
            atom_index = 0
            for j in range(len(posMol)):
                tatomlab += 1
                list = ["HETATM", str(tatomlab), nameAtom[atom_index], molname.zfill(3), str(mollab), '{:.3f}'.format(posMol[j][0]), '{:.3f}'.format(posMol[j][1]), '{:.3f}'.format(posMol[j][2]), "1.00", "0.00", nameAtom[atom_index]]
                print('{0[0]:<6}{0[1]:>5}{0[2]:>3}{0[3]:>6}{0[4]:>6}{0[5]:>12}{0[6]:>8}{0[7]:>8}{0[8]:>6}{0[9]:>6}{0[10]:>12}'.format(list), file=f)
        # ATOM      1  H   UNK     1     -12.899  32.293   3.964  1.00  0.00           H
                atom_index += 1

        print("END", file=f)
        f.close()

    def Exportspecificpos(self, path, iname, Rec, mollist, uobj, centermol):

        if os.path.exists(path) is False:
            print(path)
            subprocess.call(["mkdir", path])

        mollist.append(centermol)
        # print mollist
        # # Export position of mol
        out_file = path + "/" + str(iname) + ".xyz"
        # print out_file

        numlist = []
        # get total atom
        for i in mollist:
            numAtm = uobj.size("Structure.Position.mol[].atom[]", [i])
            numlist.append(numAtm)
        totalAtm = sum(numlist)
        # print "totalAtom", totalAtm

        f = open(out_file, "w")
        print(totalAtm, file=f)
        print(str(iname), file=f)
        f.close()

        f = open(out_file, "a+")
        uobj.jump(Rec)

        for i in mollist:
            Molnum = i
            posMol = self.getposmol(uobj, Molnum)
            nameAtom = self.getnameAtom(uobj, Molnum)
            atom_index = 0
            for j in range(len(posMol)):
                print(nameAtom[atom_index], \
                    posMol[j][0], posMol[j][1], posMol[j][2], file=f)
                atom_index += 1
        f.close()

        subprocess.call(["babel", "-ixyz", out_file, "-opdb",
                         path + "/" + str(iname) + ".pdb"])

    def moveintocell(self, uobj, totalRec, totalMol,):
        # # Move into cell
        for k in range(totalRec):
            uobj.jump(k)
            cell = uobj.get("Structure.Unit_Cell.Cell_Size")
            for i in range(totalMol):
                Molnum = i
                transVec = np.array([0., 0., 0.])  # x,y,z
                posMol = self.getposmol(uobj, Molnum)
                centerOfMol = self.getCenter(posMol)
                for j in range(3):
                    if centerOfMol[j] > cell[j]:
                        while centerOfMol[j] > cell[j]:
                            centerOfMol[j] = centerOfMol[j] - cell[j]
                            transVec[j] = transVec[j] - cell[j]
                    elif centerOfMol[j] < 0:
                        while centerOfMol[j] < 0:
                            centerOfMol[j] = centerOfMol[j] + cell[j]
                            transVec[j] = transVec[j] + cell[j]

                posMol = self.moveMolTrans(posMol, transVec)
                self.putPositionsMol(uobj, Molnum, posMol)
        print("move_done.")

    def putnvtnewfile(self, uobj, Rec, iname, addname):
        head, ext = os.path.splitext(str(iname))
        oname = head + addname + ".udf"
        uobj.jump(Rec)
        uobj.write(oname, record=-1, define=1)
        uobj.write(oname, currentRecord, append)
        uobj.put("NVT_Nose_Hoover",
                 "Simulation_Conditions.Solver.Dynamics.Dynamics_Algorithm")
        uobj.write()

    def getudfinfowrap(self, uobj):
        totalMol, totalRec = self.gettotalmol_rec(uobj)
        totalAtm = self.gettotalAtm(uobj)
        cell = self.getcellsize(uobj, totalRec-1)

        print('totalMol:', totalMol)
        print('totalAtm:', totalAtm)
        print('cellinfo:', cell)

        molnamelist = self.getnamelist([i for i in range(totalMol)], uobj, totalMol)
        import collections
        print(collections.Counter(molnamelist))

        return totalMol, totalRec, totalAtm, cell

    def gettotalmol_rec(self, uobj):
        totalMol = uobj.size("Set_of_Molecules.molecule[]")
        totalRec = uobj.totalRecord()
        return totalMol, totalRec

    def gettotalAtm(self, uobj):
        totalMol = uobj.size("Set_of_Molecules.molecule[]")
        numlist = []
        for i in range(totalMol):
            numAtm = uobj.size("Set_of_Molecules.molecule[].atom[]", [i])
            numlist.append(numAtm)
        totalAtm = sum(numlist)
        return totalAtm

    def getcellsize(self, uobj, rec):
        uobj.jump(rec)
        cell = uobj.get("Structure.Unit_Cell.Cell_Size")
        return cell

    def getmolatomnum(self, uobj, totalMol):
        atmnumlist = []
        for i in range(totalMol):
            molnum = uobj.size("Set_of_Molecules.molecule[" + str(i) + "].atom[]")
            atmnumlist.append(molnum)
        # print atmnumlist
        return atmnumlist

#     def getdist(self, p1, p2):
#         dist = math.sqrt(sum((p1 - p2)**2))
#         return dist

    def getinteractionsitetable(self, uobj, indexatom):
        # UDF operation
        # get the position of a molecule
        name = uobj.get("Molecular_Attributes.Interaction_Site_Type[].Name")
        site = uobj.get("Molecular_Attributes.Interaction_Site_Type[].Range")

        return [name, site]

    def calcMMinteraction(self, index, posMol, typenameMol, molnamelist,
                          clist, clu_num, fname, uobj):
        # index: contact mol id total info
        # posmol: all mol position list
        # typenameMol: atom type in mol
        # molnamelist: molname list for contact
        # clist: contact list per mol(renum)
        # clu_num: cluster_num
        print("******calc MM****")
        typenamelist = []
        poslist = []
        chglist = []
        # get atomtype list, position list, charge list
        for i in index:
            typenamelist.append(typenameMol[i])
            poslist.append(posMol[i])

        for i in range(len(molnamelist)):
            chglist.append(self.getchg(
                "monomer/" + molnamelist[i] + ".out", len(typenamelist[i]), 0))

    #    print typenamelist
    #    print poslist
    #    print chglist

        ielist = []
        for i in range(len(clist)):
            for j in range(1, len(clist[i])):
                # get interaction between clist[i][0] and clist[i][k]
                mol1id = clist[i][0]
                mol2id = clist[i][j]
                LJsum = 0
                coulombsum = 0
                # calc interaction energy per mol
                for k in range(len(typenamelist[mol1id])):
                    for l in range(len(typenamelist[mol2id])):
                        atomname1 = typenamelist[mol1id][k]
                        atomname2 = typenamelist[mol2id][l]
                        pos1 = poslist[mol1id][k]
                        pos2 = poslist[mol2id][l]
                        q1 = chglist[mol1id][k]
                        q2 = chglist[mol2id][l]
                        SIparam = self.getsigmaepsilon(atomname1, atomname2, uobj)
                        LJ = self.calcLJPairInteraction(pos1, pos2, SIparam)
                        coulomb = self.calcCoulombInteraction(pos1, pos2, q1, q2, 1.0)
                        LJsum += LJ
                        coulombsum += coulomb
                ielist.append([LJsum, coulombsum])
        ielist = np.array(ielist)
        print("ielist:", ielist)
    #
    #    #get sigma,epsilon from atomtypelist
    #
    #    #calc
    #    calcLJPairInteraction
    #    calcCoulombInteraction

    def getnamelist(self, index, uobj, totalMol):
        # --get used mol infomation--
        molnamelist = []
        for i in index:
            molnamelist.append(
                    uobj.get("Set_of_Molecules.molecule[" +
                             str(i % totalMol) + "].Mol_Name"))
        # print molnamelist

        return molnamelist

    def getcontactstructure(self, rec, uobj, totalMol, inmol, path, molname):
        uobj.jump(rec)
        cell = uobj.get("Structure.Unit_Cell.Cell_Size")
        print("totalmol:", totalMol, "rec", rec)

        posMol_orig = []
        typenameMol_orig = []
        elemMol_orig = []

        for i in range(totalMol):
            posMol_orig.append(self.getposmol(uobj, i))
            typenameMol_orig.append(self.getAtomtypename(uobj, i))
            elemMol_orig.append(self.getnameAtom(uobj, i))

        vec = [[0, 0, 0]]
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    vec.append([i, j, k])
        vec = np.array(vec)

        # print vec
        print("cellsize", cell)

        posMol = []
        typenameMol = []
        elemMol = []
        # print posMol_orig[0]
        for i in range(len(vec)):
            for j in range(totalMol):
                posMol.append(self.moveMolTrans(
                    posMol_orig[j], vec[i] * float(cell[0])))
                typenameMol.append(typenameMol_orig[j])
                elemMol.append(elemMol_orig[j])

        # --definition--
        # posMol: molecular position list (27cell)
        # posMol_orig: molecular position list (original cell)
        # typenameMol: Molecular typename list (27cell)

        # ---get contact_mol---
        isitelist = self.getinteractionsitetable(uobj, 1)
        # print isitelist

        site = []
        # print typenameMol[neighborMol[0][0]]
        for i in range(len(posMol)):
            site.append(self.getatomisite(isitelist, typenameMol[i]))

        # print "vec",vec
        print("pos27molnum", len(posMol), end=' ')
        print("pos_orig molnum", len(posMol_orig))
        print("typenamemol len", len(typenameMol))
        # print posMol

        posfrag_mols, typenamefrag_mols, sitefrag_mols, fragids, infrag = \
            self.getfraginfomb(molname, posMol, typenameMol,
                               site, len(posMol_orig))
        infrag = infrag * inmol
        print("infrag:", infrag)

        # centerOfMol: com of each molecules
        centerOfMol = []
        for i in range(len(posMol)):
            centerOfMol.append(self.getCenter(posMol[i]))

        # export com
        self.exportxyz(path[0] + "/" + path[1], centerOfMol, "com")

        # check
        print("nummol", len(posMol), "numatom_mol0", len(posMol[0]))
        # posMol: num of molecules in cell
        # posMol[0]: num of atoms in Mol[0]

        # --move check--
        aaa = self.moveMolTrans(posMol[0], -centerOfMol[0])
        print(len(aaa))
        print("center", self.getCenter(aaa))

        # --move all mol to origin--
        originmol = []
        for i in range(len(posMol)):
            originmol.append(self.moveMolTrans(posMol[i], -centerOfMol[i]))

        # for i in range(len(movemol)):
        #    exportxyz(path[0] + "/" + path[1],movemol[i],i)

        # --get mol radius --
        radius = []
        for i in range(len(posMol)):
            radius.append(self.getmolradius(originmol[i]))

        # --get com dist--
        distlist = []
        for i in range(inmol):
            dist = []
            for j in range(i+1, len(posMol)):
                dist.append(self.getdist(centerOfMol[i], centerOfMol[j]))
            distlist.append(dist)

        # inmol: mol num need for check contact
        # this area can be parrallel tuning
        # compare com dist and radius
        neighborMol = []
        for i in range(inmol):
            k = 0
            for j in range(i+1, len(posMol)):
                if distlist[i][k] < (radius[i] + radius[j]) * 2:
                    neighborMol.append([i, j])
                k += 1

        # get contact list
        clistall = self.getcontactlist(inmol, posMol, site, neighborMol)
        clistfrag = self.getcontactfrag(clistall, posfrag_mols,
                                        sitefrag_mols, fragids, infrag)

        index = self.getindex(clistall)
        frag_index = self.getindex(clistfrag)
        molnamelist = self.getnamelist(index, uobj, totalMol)

        # --export contact pos--
    #    for i in range(inmol):
    #        opath = path[0] + "/" + path[1] + "/contact"
    #        Exportardpos(opath, rec, clistall[i], posMol,typenameMol)

        # --export whole orig pos--
        # oname = "mdout_orig" + str(rec) + ".xyz"
        # Exportpos(path[0] + "/" + path[1],totalRec-1,totalMol,uobj,oname)

        # --export pos for abinit --
        opath = path[0] + "/" + path[1] + "/pdb"
        if self.mixflag is True or self.clusterflag is True:
            # contact
            oname = rec
        else:
            # plus 1 layer
            oname = "mdout"
        self.Exportardpos(opath, oname, index, posMol, elemMol)

        index_renum, clistall = self.getrenumindex(index, clistall)
        fragindex_renum, clistall = \
            self.getrenumfrag(frag_index, clistfrag, fragids)

        # --export clistall--
        opath = path[0] + "/" + path[1] + "/contactfrag"
        self.exportdata(opath, oname, clistall)

        # calc MM interacion
        # self.calcMMinteraction(index,posMol,typenameMol,molnamelist,clistall,inmol,molname,uobj)

        # bind structure
        # molnamelist: name list for each molecules
        self.make_abinputmb(molname, molnamelist, oname, path)

    def getsigmaepsilon(self, atom1, atom2, uobj):
        # size=uobj.size("Interactions.Pair_Interaction[]")
        aaa = uobj.get("Interactions.Pair_Interaction[]")
        for list in aaa:
            if (list[2] == atom1 and list[3] == atom2) or \
               (list[2] == atom2 and list[3] == atom1):
                # print atom1,atom2 ,list[2],list[3],list[6],"--> ok"
                break
        return list[6]
        # list[6]: [sigma(/2^(1/6)), epsilon(sqrt(e1 * e2))]
