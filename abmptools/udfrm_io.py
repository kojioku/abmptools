from __future__ import annotations

import sys
import os
import math
import subprocess
import re
import logging
from multiprocessing import Pool
from typing import Any

from .udf_io import udf_io as uio

logger = logging.getLogger(__name__)
try:
    from UDFManager import *
except ImportError:
    pass
try:
    import numpy as np
except ImportError:
    pass


class udfrm_io(uio):
    def __init__(self) -> None:
        super().__init__()
        self.molflag: bool = False
        self.cell: list[float] | None = None

    def run_convert(self, args: tuple[str, int, int, bool]) -> None:
        fname, tgtrec, tgtmol, moveflag = args
        _udf_ = UDFManager(fname)

        totalMol, totalRec = self.gettotalmol_rec(_udf_)
        logger.info("totalRec = %s", totalRec)
        logger.info("totalMol = %s", totalMol)
        if tgtrec == -1:
            tgtrec = totalRec -1

        if tgtmol == -1:
            tgtmol = totalMol

        logger.info("targetRec = %s", tgtrec)
        logger.info("moveflag = %s", moveflag)

        if moveflag:
            self.moveintocell_rec(_udf_, tgtrec, totalMol)

        oname = os.path.splitext(fname)[0].split("/")[-1]
        logger.debug("%s", oname)
        if self.molflag:
            self.convert_udf_pdb(tgtrec, _udf_, tgtmol, oname)
        else:
            self.convert_udf_pdb(tgtrec, _udf_, totalMol, oname)

    def convert_udf_pdb(self, rec: int, uobj: Any, totalMol: int, ohead: str, writef: bool = True) -> list[list[Any]]:
        uobj.jump(rec)
        self.cell = uobj.get("Structure.Unit_Cell.Cell_Size")
        logger.debug("totalMol=%s, rec=%s", totalMol, rec)
        # getmolatomnum(_udf_, totalMol)

        posMol = []
        typenameMol = []
        molnamelist = []

        if self.molflag:
            tgtmol = totalMol
            logger.debug("tgtmol = %s", tgtmol)
            for i in range(tgtmol, tgtmol + 1):
                posMol.append(self.getposmol(uobj, i))
                typenameMol.append(self.getAtomtypename(uobj, i))
                molnamelist.append(self.getmolname(i, uobj))

            # print (typenameMol)
            # print (posMol)
            # print (molnamelist)

            oname = ohead + ".xyz"
            # self.Exportpos(".", rec, tgtmol, uobj, oname)
            if writef:
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
            if writef:
                self.Exportpos(".", rec, totalMol, uobj, oname)

        return [typenameMol, posMol, molnamelist]


    def moveintocell_rec(self, uobj: Any, Rec: int, totalMol: int) -> None:
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
        logger.info("move_done.")


    def getmolname(self, i: int, uobj: Any) -> str:
        # --get used mol infomation--
        molname = uobj.get("Set_of_Molecules.molecule[" +
                           str(i) + "].Mol_Name")
        # print molnamelist

        return molname

