from __future__ import annotations
"""UDFファイルからPDB/XYZ形式への変換を行うモジュール。

UDFManager経由で分子座標を読み込み、セル内移動や
PDB/XYZ形式でのエクスポートを提供する。
"""

import os
import logging
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
    """UDFファイルの座標データをPDB/XYZ形式に変換するクラス。

    指定レコード・分子の座標をセル内に移動し、XYZ/PDBファイルとして出力する。
    """

    def __init__(self) -> None:
        super().__init__()
        self.molflag: bool = False
        self.cell: list[float] | None = None

    def run_convert(self, args: tuple[str, int, int, bool]) -> None:
        """UDFファイルを読み込み、指定条件でPDB/XYZ形式に変換する。

        Args:
            args: (ファイル名, 対象レコード番号, 対象分子番号, セル内移動フラグ) のタプル。
                  レコード/分子番号に-1を指定すると最終レコード/全分子を対象とする。
        """
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
        """UDFの指定レコードから分子座標・原子タイプ名を取得し、XYZファイルに出力する。

        Args:
            rec: 対象レコード番号。
            uobj: UDFManagerオブジェクト。
            totalMol: 対象分子数。
            ohead: 出力ファイル名のプレフィックス。
            writef: Trueの場合ファイルに書き出す。

        Returns:
            [原子タイプ名リスト, 座標リスト, 分子名リスト] のリスト。
        """
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
        """指定レコードの全分子をセル内に移動する。

        各分子の重心がセル範囲外にある場合、セルベクトル分だけ平行移動する。

        Args:
            uobj: UDFManagerオブジェクト。
            Rec: 対象レコード番号。
            totalMol: 全分子数。
        """
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
        """UDFから指定インデックスの分子名を取得する。

        Args:
            i: 分子のインデックス番号。
            uobj: UDFManagerオブジェクト。

        Returns:
            str: 分子名。
        """
        # --get used mol infomation--
        molname = uobj.get("Set_of_Molecules.molecule[" +
                           str(i) + "].Mol_Name")
        # print molnamelist

        return molname

