import sys
import os
import math
import subprocess
import re
from multiprocessing import Pool
import udf_io as uio
try:
    from UDFManager import *
except:
    pass
try:
    import numpy as np
except:
    pass


class udfrm_io(uio.udf_io):
    def __init__(self):
        self.molflag = False
        self.cell = None
        pass

    def run_convert(self, args):
        fname, tgtrec, tgtmol, moveflag = args
        _udf_ = UDFManager(fname)

        totalMol, totalRec = self.gettotalmol_rec(_udf_)
        print("totalRec = ", totalRec)
        print("totalMol = ", totalMol)
        if tgtrec == -1:
            tgtrec = totalRec -1

        if tgtmol == -1:
            tgtmol = totalMol

        print("targetRec = ", tgtrec)
        print("moveflag = ", moveflag)

        if moveflag == True:
            self.moveintocell_rec(_udf_, tgtrec, totalMol)

        oname = os.path.splitext(fname)[0].split("/")[-1]
        print(oname)
        if self.molflag == True:
            self.convert_udf_pdb(tgtrec, _udf_, tgtmol, oname)
        else:
            self.convert_udf_pdb(tgtrec, _udf_, totalMol, oname)

    def convert_udf_pdb(self, rec, uobj, totalMol, ohead, writef=True):
        uobj.jump(rec)
        self.cell = uobj.get("Structure.Unit_Cell.Cell_Size")
        print(totalMol, rec)
        # getmolatomnum(_udf_, totalMol)

        posMol = []
        typenameMol = []
        molnamelist = []

        if self.molflag == True:
            tgtmol = totalMol
            print('tgtmol =', tgtmol)
            for i in range(tgtmol, tgtmol + 1):
                posMol.append(self.getposmol(uobj, i))
                typenameMol.append(self.getAtomtypename(uobj, i))
                molnamelist.append(self.getmolname(i, uobj))

            # print (typenameMol)
            # print (posMol)
            # print (molnamelist)

            oname = ohead + ".xyz"
            # self.Exportpos(".", rec, tgtmol, uobj, oname)
            if writef == True:
                self.Exporttgtmolpos(".", oname, rec, [tgtmol], uobj)

        else:
            for i in range(totalMol):
                posMol.append(self.getposmol(uobj, i))
                typenameMol.append(self.getAtomtypename(uobj, i))
                molnamelist.append(self.getmolname(i, uobj))

            # print (typenameMol)
            # print (posMol)
            # print (molnamelist)

            oname = ohead + ".xyz"
            if writef == True:
                self.Exportpos(".", rec, totalMol, uobj, oname)

        return [typenameMol, posMol, molnamelist]


    def moveintocell_rec(self, uobj, Rec, totalMol):
        # # Move into cell
        uobj.jump(Rec)
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


    def getmolname(self, i, uobj):
        # --get used mol infomation--
        molname = uobj.get("Set_of_Molecules.molecule[" +
                           str(i) + "].Mol_Name")
        # print molnamelist

        return molname

