"""UDFファイルの入出力および分子座標操作を提供するモジュール。"""

from __future__ import annotations

import sys
import os
import subprocess
import re
import time
import copy
import logging
from typing import Any
from .molcalc import molcalc as molc

logger = logging.getLogger(__name__)
try:
    from UDFManager import *
except ImportError:
    pass
try:
    import numpy as np
except ImportError:
    pass


class udf_io(molc):
    """UDFファイルの読み書きと分子座標の取得・設定・変換を行うクラス。"""

    def __init__(self) -> None:
        super().__init__()
        self.verflag = True

    def getposatom(self, uobj: Any, indexatom: int) -> np.ndarray:
        """指定原子インデックスの座標を取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            indexatom: 原子のインデックス。

        Returns:
            原子の座標を格納したnumpy配列。
        """
        # UDF operation
        # get the position of a molecule
        i = indexatom
        list = uobj.get("Structure.Position.mol[].atom[]", [i])
        position = np.array(list)

        return position

    def getposmolrec(self, uobj: Any, indexMol: int, record: int) -> np.ndarray:
        """指定レコードにおける分子の原子座標を取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            indexMol: 分子のインデックス。
            record: 対象レコード番号。

        Returns:
            分子内の全原子座標を格納したnumpy配列。
        """
        # UDF operation
        # get the position of a molecule
        uobj.jump(record)
        i = indexMol
        aaa = uobj.get("Structure.Position.mol[" + str(i) + "].atom[]")
        position = np.array(aaa)

        return position

    def getposmol(self, uobj: Any, indexMol: int) -> np.ndarray:
        """指定分子の全原子座標を取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            indexMol: 分子のインデックス。

        Returns:
            分子内の全原子座標を格納したnumpy配列。
        """
        # UDF operation
        # get the position of a molecule
        i = indexMol
        list = uobj.get("Structure.Position.mol[" + str(i) + "].atom[]")
        position = np.array(list)

        return position

    def getnameAtom(self, uobj: Any, indexMol: int) -> list[str]:
        """指定分子内の全原子名を取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            indexMol: 分子のインデックス。

        Returns:
            原子名のリスト。
        """
        # get atom name
        i = indexMol
        atom = uobj.get("Set_of_Molecules.molecule[].atom[].Atom_Name", [i])
        return atom

    def getAtomtypename(self, uobj: Any, indexMol: int) -> list[str]:
        """指定分子内の全原子タイプ名を取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            indexMol: 分子のインデックス。

        Returns:
            原子タイプ名のリスト。
        """
        # get atom name
        i = indexMol
        atomtype = uobj.get(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [i])
        return atomtype

    def putPositionsMol(self, uobj: Any, indexMol: int, position: np.ndarray) -> None:
        """分子の原子座標をUDFに書き込む。

        Args:
            uobj: UDFManagerオブジェクト。
            indexMol: 分子のインデックス。
            position: 書き込む原子座標のnumpy配列。
        """
        # ## put the position of the molecule to UDF
        i = indexMol
        numAtm = uobj.size("Structure.Position.mol[].atom[]", [i])
        for j in range(numAtm):
            uobj.put(position[j, 0], "Structure.Position.mol[].atom[].x", [i, j])
            uobj.put(position[j, 1], "Structure.Position.mol[].atom[].y", [i, j])
            uobj.put(position[j, 2], "Structure.Position.mol[].atom[].z", [i, j])

    def Exportpos(self, path: str, Rec: int, totalMol: int, uobj: Any, oname: str) -> None:
        """全分子の座標をXYZおよびPDB形式でエクスポートする。

        Args:
            path: 出力先ディレクトリパス。
            Rec: 対象レコード番号。
            totalMol: 全分子数。
            uobj: UDFManagerオブジェクト。
            oname: 出力ファイル名。
        """
        # # Export position of mol
        # head, ext = os.path.splitext(str(iname))
        os.makedirs(path + "/pdb", exist_ok=True)

        out_file = path + "/pdb/" + str(oname)
        logger.info(out_file)

        numlist = []
        # get total atom
        for i in range(totalMol):
            numAtm = uobj.size("Structure.Position.mol[].atom[]", [i])
            numlist.append(numAtm)
        totalAtm = sum(numlist)
        logger.info("totalAtom %s", totalAtm)

        with open(out_file, "w") as f:
            print(totalAtm, file=f)
            print(out_file, file=f)

        with open(out_file, "a+") as f:
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

        # subprocess.call(["babel", "-ixyz", out_file, "-opdb",
        #                  path + "/pdb/mdout_orig.pdb"])
        self.exportpdb(uobj, Rec, out_file, mollist)

    def Exporttgtmolpos(self, path: str, oname_i: str, Rec: int, mollist: list[int], uobj: Any) -> None:
        """指定分子リストの座標をXYZおよびPDB形式でエクスポートする。

        Args:
            path: 出力先ディレクトリパス。
            oname_i: 出力ファイル名のベース。
            Rec: 対象レコード番号。
            mollist: 対象分子インデックスのリスト。
            uobj: UDFManagerオブジェクト。
        """

        os.makedirs(path, exist_ok=True)

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

        with open(out_file, "w") as f:
            print(totalAtm, file=f)
            print(str(ohead), file=f)

        with open(out_file, "a+") as f:
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

        self.exportpdb(uobj, Rec, out_file, mollist)

        # cmd = "babel -ixyz " +  out_file +  " -opdb ", path + "/" + str(ohead) + ".pdb"
        # subprocess.call(cmd, shell = True)

    def exportpdb(self, uobj: Any, Rec: int, out_file: str, mollist: list[int]) -> None:
        """指定分子リストの座標をPDB形式で出力する。

        Args:
            uobj: UDFManagerオブジェクト。
            Rec: 対象レコード番号。
            out_file: 出力ファイルパス。
            mollist: 対象分子インデックスのリスト。
        """

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
        logger.debug(molids)
        molnames = []
        for molname in molnames_orig:
            for i in range(len(molids)):
                if molname == molids[i]:
                    molnames.append(i)
        # print(molnames)
        # header
        # print(out_file)
        with open(out_file, "w", newline = "\n") as f:
            print("COMPND    " + out_file, file=f)
            print("AUTHOR    " + "GENERATED BY python script in FMOrmap", file=f)
            if self.cell is not None:
                print('CRYST1{0[0]:>9.3f}{0[1]:>9.3f}{0[2]:>9.3f}  90.00  90.00  90.00               1'.format(self.cell), file=f)

        # aaa = [0.8855]
        # print '{0[0]:.3f}'.format(aaa)
        # print pos[0]
        with open(out_file, "a+", newline = "\n") as f:
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

    def Exportspecificpos(self, path: str, iname: str, Rec: int, mollist: list[int], uobj: Any, centermol: int) -> None:
        """指定分子リストに中心分子を加えた座標をXYZおよびPDB形式でエクスポートする。

        Args:
            path: 出力先ディレクトリパス。
            iname: 出力ファイル名のベース。
            Rec: 対象レコード番号。
            mollist: 対象分子インデックスのリスト。
            uobj: UDFManagerオブジェクト。
            centermol: 中心分子のインデックス。
        """

        os.makedirs(path, exist_ok=True)

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

        with open(out_file, "w") as f:
            print(totalAtm, file=f)
            print(str(iname), file=f)

        with open(out_file, "a+") as f:
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

        subprocess.call(["babel", "-ixyz", out_file, "-opdb",
                         path + "/" + str(iname) + ".pdb"])

    def moveintocell(self, uobj: Any, totalRec: int, totalMol: int) -> None:
        """全レコード・全分子をセル内に移動させる。

        Args:
            uobj: UDFManagerオブジェクト。
            totalRec: 全レコード数。
            totalMol: 全分子数。
        """
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
        logger.info("move_done.")

    def moveintocell_mol(self, uobj: Any, totalRec: int, startmol: int, endmol: int) -> None:
        """指定範囲の分子を全レコードにわたってセル内に移動させる。

        Args:
            uobj: UDFManagerオブジェクト。
            totalRec: 全レコード数。
            startmol: 開始分子インデックス。
            endmol: 終了分子インデックス（含まない）。
        """
        # # Move into cell
        for k in range(totalRec):
            uobj.jump(k)
            cell = uobj.get("Structure.Unit_Cell.Cell_Size")
            for i in range(startmol, endmol):
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

    def putnvtnewfile(self, uobj: Any, Rec: int, iname: str, addname: str) -> None:
        """NVT Nose-Hoover設定を適用した新規UDFファイルを出力する。

        Args:
            uobj: UDFManagerオブジェクト。
            Rec: 対象レコード番号。
            iname: 入力ファイル名。
            addname: 出力ファイル名に付加する文字列。
        """
        head, ext = os.path.splitext(str(iname))
        oname = head + addname + ".udf"
        uobj.jump(Rec)
        uobj.put("NVT_Nose_Hoover",
                 "Simulation_Conditions.Solver.Dynamics.Dynamics_Algorithm")
        uobj.write(oname, record=-1, define=1)
        uobj.write(oname, currentRecord, append)

    def putnvtnewfilemb(self, uobj: Any, Rec: int, iname: str, addname: str) -> None:
        """NVTシミュレーション条件を設定しリスタート用UDFファイルを出力する。

        Args:
            uobj: UDFManagerオブジェクト。
            Rec: 対象レコード番号。
            iname: 入力ファイル名。
            addname: 出力ファイル名に付加する文字列。
        """
        head, ext = os.path.splitext(str(iname))
        logger.info("totalstep: %s", self.nvtstep)
        logger.info("outstep: %s", self.nvtoutstep)
        logger.info("Algorithm: %s", self.nvtalgo[0])
        oname = head + addname + ".udf"

        uobj.jump(-1)   # jump to common record
        uobj.put(self.nvtstep, "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps")
        uobj.put(self.nvtoutstep, "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps")
        uobj.put(self.nvtalgo[0],
                 "Simulation_Conditions.Solver.Dynamics.Dynamics_Algorithm")
        uobj.put("Restart", "Initial_Structure.Generate_Method.Method")
        uobj.put(-1, 'Initial_Structure.Generate_Method.Restart.Record')
        uobj.put(1, 'Initial_Structure.Generate_Method.Restart.Restore_Cell')
        uobj.put(1, 'Initial_Structure.Generate_Method.Restart.Restore_Velocity')

        '''
        record = allRecord(0)/currentRecord(1)/initialRecord(-1)
        mode = overwrite(0)/append(1)：
        define=0/1
        '''

        uobj.jump(Rec)   # jump to target record
        uobj.write(oname, record=-1, define=1)   # save common record
        uobj.write(oname, currentRecord, append) # save target record

    def getudfinfowrap(self, uobj: Any) -> tuple[int, int, int, list[float]]:
        """UDFファイルの基本情報（分子数、レコード数、原子数、セルサイズ）を一括取得する。

        Args:
            uobj: UDFManagerオブジェクト。

        Returns:
            (全分子数, 全レコード数, 全原子数, セルサイズ)のタプル。
        """
        totalMol, totalRec = self.gettotalmol_rec(uobj)
        totalAtm = self.gettotalAtm(uobj)
        cell = self.getcellsize(uobj, totalRec-1)

        logger.info("totalMol: %s", totalMol)
        logger.info("totalAtm: %s", totalAtm)
        logger.info("cellinfo: %s", cell)

        molnamelist = self.getnamelist([i for i in range(totalMol)], uobj, totalMol)
        import collections
        logger.info(collections.Counter(molnamelist))

        return totalMol, totalRec, totalAtm, cell

    def gettotalmol_rec(self, uobj: Any) -> tuple[int, int]:
        """UDFの全分子数と全レコード数を取得する。

        Args:
            uobj: UDFManagerオブジェクト。

        Returns:
            (全分子数, 全レコード数)のタプル。
        """
        totalMol = uobj.size("Set_of_Molecules.molecule[]")
        totalRec = uobj.totalRecord()
        return totalMol, totalRec

    def gettotalAtm(self, uobj: Any) -> int:
        """UDF内の全原子数を取得する。

        Args:
            uobj: UDFManagerオブジェクト。

        Returns:
            全原子数。
        """
        totalMol = uobj.size("Set_of_Molecules.molecule[]")
        numlist = []
        for i in range(totalMol):
            numAtm = uobj.size("Set_of_Molecules.molecule[].atom[]", [i])
            numlist.append(numAtm)
        totalAtm = sum(numlist)
        return totalAtm

    def getcellsize(self, uobj: Any, rec: int) -> list[float]:
        """指定レコードのセルサイズを取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            rec: 対象レコード番号。

        Returns:
            セルサイズのリスト。
        """
        uobj.jump(rec)
        cell = uobj.get("Structure.Unit_Cell.Cell_Size")
        return cell

    def getmolatomnum(self, uobj: Any, totalMol: int) -> list[int]:
        """各分子の原子数リストを取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            totalMol: 全分子数。

        Returns:
            各分子の原子数を格納したリスト。
        """
        atmnumlist = []
        for i in range(totalMol):
            molnum = uobj.size("Set_of_Molecules.molecule[" + str(i) + "].atom[]")
            atmnumlist.append(molnum)
        # print atmnumlist
        return atmnumlist

#     def getdist(self, p1, p2):
#         dist = math.sqrt(sum((p1 - p2)**2))
#         return dist

    def getinteractionsitetable(self, uobj: Any, indexatom: int) -> list[list[Any]]:
        """相互作用サイトの名前と範囲のテーブルを取得する。

        Args:
            uobj: UDFManagerオブジェクト。
            indexatom: 原子のインデックス。

        Returns:
            [名前リスト, 範囲リスト]の2要素リスト。
        """
        # UDF operation
        # get the position of a molecule
        name = uobj.get("Molecular_Attributes.Interaction_Site_Type[].Name")
        site = uobj.get("Molecular_Attributes.Interaction_Site_Type[].Range")

        return [name, site]

    def calcMMinteraction(self, index: list[int], posMol: list[np.ndarray], typenameMol: list[list[str]], molnamelist: list[str],
                          clist: list[list[int]], clu_num: int, fname: str, uobj: Any) -> None:
        """コンタクトリストに基づきLJおよびクーロン分子間相互作用エネルギーを計算する。

        Args:
            index: コンタクト分子のインデックスリスト。
            posMol: 全分子の原子座標リスト。
            typenameMol: 各分子の原子タイプ名リスト。
            molnamelist: コンタクト分子の名前リスト。
            clist: 分子ごとのコンタクトリスト。
            clu_num: クラスタ番号。
            fname: 分子名。
            uobj: UDFManagerオブジェクト。
        """
        # index: contact mol id total info
        # posmol: all mol position list
        # typenameMol: atom type in mol
        # molnamelist: molname list for contact
        # clist: contact list per mol(renum)
        # clu_num: cluster_num
        logger.info("******calc MM****")
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
        logger.debug("ielist: %s", ielist)
    #
    #    #get sigma,epsilon from atomtypelist
    #
    #    #calc
    #    calcLJPairInteraction
    #    calcCoulombInteraction

    def getnamelist(self, index: list[int], uobj: Any, totalMol: int) -> list[str]:
        """指定インデックスの分子名リストを取得する。

        Args:
            index: 分子インデックスのリスト。
            uobj: UDFManagerオブジェクト。
            totalMol: 全分子数（周期境界条件での剰余計算に使用）。

        Returns:
            分子名のリスト。
        """
        # --get used mol infomation--
        molnamelist = []
        for i in index:
            molnamelist.append(
                    uobj.get("Set_of_Molecules.molecule[" +
                             str(i % totalMol) + "].Mol_Name"))
        # print molnamelist

        return molnamelist

    def getcontactstructure(self, rec: int, uobj: Any, totalMol: int, seg1_clunum: int, path: list[str], molname: str,
                            lammpsflag: bool = False, masses: Any = None,
                            atom_molecule_map: Any = None, tgtfile: str | None = None,
                            molecule_list: list[Any] | None = None) -> None:
        '''
        seg1_clunum(for center mode): number of molecules in the first cluster
        seg1_clunum(for whole mode):  number of moleclues in the system
        '''

        if not lammpsflag:
            uobj.jump(rec)
            cell = uobj.get("Structure.Unit_Cell.Cell_Size")
            logger.info("totalmol: %s rec %s", totalMol, rec)

            posMol_orig = []
            typenameMol_orig = []
            elemMol_orig = []

            for i in range(totalMol):
                posMol_orig.append(self.getposmol(uobj, i))
                typenameMol_orig.append(self.getAtomtypename(uobj, i))
                elemMol_orig.append(self.getnameAtom(uobj, i))

        if lammpsflag:
            logger.info("LAMMPS mode")
            # require
            # posMol_orig
            # typenameMol_orig
            # elemMol_orig
            # cell
            # Example usage for trajectory file
            logger.debug(tgtfile)
            trajectory_file_path = tgtfile
            timesteps, box_bounds = \
                self.parse_lammps_trajectory(trajectory_file_path)

            for timestep, atom_data in timesteps.items():
                atom_data.sort(key=lambda x: x[0])  # 追加: 原子IDでソート
                real_coords = self.scale_to_real_coords(
                    atom_data, box_bounds, atom_molecule_map)
                molecule_real_coords = \
                    self.group_real_coords_by_molecule(real_coords)

                for molecule_tag in molecule_real_coords:
                    molecule_real_coords[molecule_tag] = \
                        self.shift_molecule_to_primary_cell(
                        molecule_real_coords[molecule_tag],
                        box_bounds)  # 分子全体を本セル内にシフト
                posMol_orig = list(molecule_real_coords.values())
                # print(f"Molecule real coordinates for timestep {timestep}:")
                # print(posMol_orig)
                # print(len(posMol_orig))

                molecule_radii = self.get_radii_for_molecules(real_coords, masses)
                radii_list = list(molecule_radii.values())
                # print(f"Molecule radii for timestep {timestep}:")
                # print(radii_list)

                # cell
                # elemMol_orig
                elemMol_orig = molecule_list
                # typenameMol_orig
                typenameMol_orig = radii_list

            cell = [box_bounds[1]-box_bounds[0], box_bounds[3]-box_bounds[2],
                    box_bounds[5]-box_bounds[4]]

        vec = [[0, 0, 0]]
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    vec.append([i, j, k])
        vec = np.array(vec)

        # print vec
        logger.info("cellsize %s", cell)

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

        if not lammpsflag:
            isitelist = self.getinteractionsitetable(uobj, 1)
            # print isitelist

            site = []
            # print typenameMol[neighborMol[0][0]]
            for i in range(len(posMol)):
                site.append(self.getatomisite(isitelist, typenameMol[i]))
        else:
            site = copy.deepcopy(typenameMol)

        # print('site\n',site)
        # sys.exit()

        # print "vec",vec
        logger.debug("pos27molnum %s pos_orig molnum %s", len(posMol), len(posMol_orig))
        logger.debug("typenamemol len %s", len(typenameMol))
        # print posMol

        posfrag_mols, typenamefrag_mols, sitefrag_mols, fragids, infrag = \
            self.getfraginfomb(molname, posMol, typenameMol,
                               site, len(posMol_orig), seg1_clunum)
        infrag = infrag * seg1_clunum
        logger.info("seg1_frag: %s", infrag)
        # sys.exit()

        # centerOfMol: com of each molecules
        centerOfMol = []
        for i in range(len(posMol)):
            centerOfMol.append(self.getCenter(posMol[i]))

        # export com
        self.exportxyz(path[0] + "/" + path[1], centerOfMol, "com")

        # check
        logger.debug("nummol %s numatom_mol0 %s", len(posMol), len(posMol[0]))
        # posMol: num of molecules in cell
        # posMol[0]: num of atoms in Mol[0]

        # --move check--
        aaa = self.moveMolTrans(posMol[0], -centerOfMol[0])
        # print(len(aaa))
        logger.debug("center %s", self.getCenter(aaa))

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

        # --- neighbor selection (cutmode-aware) ---
        # contact_cutmode は udfcreate.setudfparam で input_param['cutmode']
        # から設定される。
        #   'contact' (default, legacy): COM 距離 < 2 * (r_i + r_j)
        #   'around': any atom within criteria of any solute (center) atom。
        #            solute = posMol[0..seg1_clunum-1] (center cluster)。
        contact_cutmode = getattr(self, 'contact_cutmode', 'contact')
        neighborMol = []

        if contact_cutmode == 'around':
            contact_criteria = getattr(self, 'contact_criteria', 4.0)
            crit2 = contact_criteria * contact_criteria
            logger.info("neighbor selection: cutmode='around', criteria=%.2f Å",
                        contact_criteria)
            posMol_np = [np.asarray(pm, dtype=float) for pm in posMol]
            # Precompute COM array for fast pairwise filter
            com_np = np.asarray([np.asarray(c, dtype=float)
                                 for c in centerOfMol])
            radius_np = np.asarray(radius, dtype=float)
            # 中心同士 (cluster center 内の他 mol) は同じクラスタの一部として
            # 常に contact 扱い。clusterflag=True (multi-mol center) のみ
            # 意味を持つ。whole mode (seg1_clunum=totalMol) では適用しない。
            cluster_internal = bool(getattr(self, 'clusterflag', False))
            for i in range(seg1_clunum):
                sol_atoms = posMol_np[i]
                ri = radius[i]
                # vectorized COM distance from solute i to ALL candidates
                diff_com = com_np[i+1:] - com_np[i]
                d_com = np.sqrt(np.einsum('ij,ij->i', diff_com, diff_com))
                cutoff = ri + radius_np[i+1:] + contact_criteria
                candidate_mask = d_com <= cutoff
                for offset, hit in enumerate(candidate_mask):
                    j = i + 1 + offset
                    if cluster_internal and j < seg1_clunum:
                        # 同一 cluster center 内 → 距離問わず contact
                        neighborMol.append([i, j])
                        continue
                    if not hit:
                        continue
                    diff = sol_atoms[:, None, :] - posMol_np[j][None, :, :]
                    d2 = np.einsum('ijk,ijk->ij', diff, diff)
                    if d2.min() < crit2:
                        neighborMol.append([i, j])
        else:
            # legacy 'contact' mode: COM dist < 2 * (r_i + r_j)
            distlist = []
            for i in range(seg1_clunum):
                dist = []
                for j in range(i+1, len(posMol)):
                    dist.append(self.getdist(centerOfMol[i], centerOfMol[j]))
                distlist.append(dist)
            for i in range(seg1_clunum):
                k = 0
                for j in range(i+1, len(posMol)):
                    if distlist[i][k] < (radius[i] + radius[j]) * 2:
                        neighborMol.append([i, j])
                    k += 1

        # get contact list
        # around mode: neighborMol は既に criteria 距離で確定済み (cluster
        # center 内 pair も無条件で追加済) なので、 getcontactlist の
        # vdW 距離 (r1+r2)/1.5 による 2nd-pass 再 filter は SKIP する。
        # contact mode (legacy): 従来通り getcontactlist で site-level 距離
        # チェックを行う。
        if contact_cutmode == 'around':
            contactlist = list(neighborMol)
            clistall_raw = []
            for i in range(seg1_clunum):
                clistmol = []
                for j in contactlist:
                    if j[1] == i:
                        clistmol.append(j[0])
                    if j[0] == i:
                        clistmol.append(j[1])
                clistall_raw.append(clistmol)
            for i in range(len(clistall_raw)):
                clistall_raw[i].insert(0, i)
            clistall = clistall_raw
        else:
            clistall = self.getcontactlist(seg1_clunum, posMol, site, neighborMol)
        logger.debug("clistall %s", clistall)
        clistfrag = self.getcontactfrag(clistall, posfrag_mols,
                                        sitefrag_mols, fragids, infrag)
        logger.debug("clistfrag %s", clistfrag)
        # sys.exit()

        index = self.getindex(clistall)
        frag_index = self.getindex(clistfrag)
        molnamelist = self.getnamelist(index, uobj, totalMol)
        # print("molnamelist", molnamelist)

        # --export contact pos--
    #    for i in range(seg1_clunum):
    #        opath = path[0] + "/" + path[1] + "/contact"
    #        Exportardpos(opath, rec, clistall[i], posMol,typenameMol)

        # --export whole orig pos--
        # oname = "mdout_orig" + str(rec) + ".xyz"
        # Exportpos(path[0] + "/" + path[1],totalRec-1,totalMol,uobj,oname)

        # --- minimum-image wrap (BUG fix 2026-05-31) ---
        # 27-cell expansion (posMol) で隣接画像セルから選ばれた neighbor は、
        # 元コードでは posMol[i] = posMol_orig[i%N] + vec*L をそのまま PDB 出力。
        # 結果: center mode で中央 cluster の周りに「+L / -L 方向に飛び離れた
        # phantom cluster」が 4-6 個現れ、可視化上分子が散乱、FMO/SCC でも
        # 本来近接する neighbor が遠方の独立 monomer として処理され収束しない。
        #
        # 修正: center mode のとき、各 neighbor を center (= posMol[0]) の
        # 最近接画像 (minimum-image convention) になるように整数セル分シフト。
        # 結果として PDB の cluster は中央に集約された「shell」形状になる。
        #
        # whole mode はシフトしない: contactfrag に image-cell index (i >= totalMol)
        # が含まれ、FMO は posMol[i] を image cell の geometric position として
        # 解釈するため、原セルに collapse すると duplicate atom と化して
        # 振動の原因になる。
        if self.mixflag or self.clusterflag:
            try:
                L = np.array([float(cell[0]), float(cell[1]), float(cell[2])])
            except (TypeError, IndexError, ValueError):
                Lscalar = float(cell[0]) if hasattr(cell, '__getitem__') else float(cell)
                L = np.array([Lscalar, Lscalar, Lscalar])

            ref_com = np.array(self.getCenter(posMol[0]))
            for i in index:
                if i < seg1_clunum:
                    continue  # center 自身はシフトしない
                com = np.array(self.getCenter(posMol[i]))
                shift = -L * np.round((com - ref_com) / L)
                if np.any(shift != 0):
                    posMol[i] = [[p[0] + shift[0], p[1] + shift[1],
                                  p[2] + shift[2]] for p in posMol[i]]

        # --export pos for abinit --
        opath = path[0] + "/" + path[1] + "/pdb"
        if self.mixflag or self.clusterflag:
            # contact
            oname = rec
        else:
            # plus 1 layer
            oname = "mdout"
        # `molnamelist` を Exportardpos に渡し PDB の residue name を正しく
        # 派生させる (obabel bond perception を介さない direct writer)。
        self.Exportardpos(opath, oname, index, posMol, elemMol,
                          molnames=molnamelist)

        index_renum, clistall = self.getrenumindex(index, clistall)
        fragindex_renum, clistall = \
            self.getrenumfrag(frag_index, clistfrag, fragids)

        # --export clistall--
        opath = path[0] + "/" + path[1] + "/contactfrag"
        self.exportdata(opath, oname, clistall)

        # calc MM interacion
        # self.calcMMinteraction(index,posMol,typenameMol,molnamelist,clistall,seg1_clunum,molname,uobj)

        # bind structure
        # molnamelist: name list for each molecules
        self.make_abinputmb(molname, molnamelist, oname, path)

    def getsigmaepsilon(self, atom1: str, atom2: str, uobj: Any) -> list[float]:
        """指定原子ペアのLJパラメータ（sigma, epsilon）を取得する。

        Args:
            atom1: 第1原子のタイプ名。
            atom2: 第2原子のタイプ名。
            uobj: UDFManagerオブジェクト。

        Returns:
            [sigma, epsilon]のリスト。
        """
        # size=uobj.size("Interactions.Pair_Interaction[]")
        aaa = uobj.get("Interactions.Pair_Interaction[]")
        for list in aaa:
            if (list[2] == atom1 and list[3] == atom2) or \
               (list[2] == atom2 and list[3] == atom1):
                # print atom1,atom2 ,list[2],list[3],list[6],"--> ok"
                break
        return list[6]
        # list[6]: [sigma(/2^(1/6)), epsilon(sqrt(e1 * e2))]

    def getatomtype(self, uobj: Any) -> list[Any]:
        """UDFから全原子タイプ情報を取得する。

        Args:
            uobj: UDFManagerオブジェクト。

        Returns:
            原子タイプ情報のリスト。
        """
        atomtype = uobj.get("Molecular_Attributes.Atom_Type[]")
        return atomtype
