import numpy as np
from UDFManager import *
import sys
import os
import math
import subprocess
import re
# Matrix operation


def moveMolTrans(posVec, transVec):
    # Parallel shift
    posVec = posVec + transVec
    return posVec


def getposmol(uobj, indexMol):
    # UDF operation
    # get the position of a molecule
        i = indexMol
        list = uobj.get("Structure.Position.mol[" + str(i) + "].atom[]")
        position = np.array(list)

        return position


def putPositionsMol(uobj, indexMol, position):
    # ## put the position of the molecule to UDF
    i = indexMol
    numAtm = uobj.size("Structure.Position.mol[].atom[]", [i])
    for j in range(numAtm):
        uobj.put(position[j, 0], "Structure.Position.mol[].atom[].x", [i, j])
        uobj.put(position[j, 1], "Structure.Position.mol[].atom[].y", [i, j])
        uobj.put(position[j, 2], "Structure.Position.mol[].atom[].z", [i, j])


def gettotalmol_rec(uobj):
    totalMol = uobj.size("Set_of_Molecules.molecule[]")
    totalRec = uobj.totalRecord()
    return totalMol, totalRec


if __name__ == "__main__":
    # main
    argvs = sys.argv
    fname = str(argvs[1])
    xx = argvs[2]
    yy = argvs[3]
    zz = argvs[4]
    try:
        molid = int(argvs[5])
    except:
        molid = 'all'

    head, ext = os.path.splitext(fname)
    oname = head + '_move' + ext
    _udf_ = UDFManager(fname)

    totalMol, totalRec = gettotalmol_rec(_udf_)
    transVec = np.array([xx, yy, zz], dtype='float64')  # x,y,z

    for j in range(totalRec):
        print ("Rec", str(j))
        _udf_.jump(j)

        if molid == 'all':
            for i in range(totalMol):
                posMol = getposmol(_udf_, i)
                posmove = moveMolTrans(posMol, transVec)
                putPositionsMol(_udf_, i, posmove)
        else:
            posMol = getposmol(_udf_, molid)
            posmove = moveMolTrans(posMol, transVec)
            putPositionsMol(_udf_, molid, posmove)

    _udf_.write(oname)
