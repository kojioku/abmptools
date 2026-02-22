# -*- coding: utf-8 -*-
"""
top_parser.py
-------------
Parses GROMACS .top / .itp files and returns raw topology data.

Public API
----------
TopParser().parse(top_path) -> TopRawData
    Parses *top_path* (resolving all #include ITP), returning a TopRawData
    that holds all information needed by TopAdapter.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# Output container for raw parsed data
# ---------------------------------------------------------------------------

@dataclass
class TopRawData:
    """
    Raw output of TopParser.parse().

    Field semantics:

    defaults
    --------
    comb_rule   : int           combining rule (2=LB, 3=geometric)
    fudge_lj    : float
    fudge_qq    : float

    Global type definitions (from [atomtypes], [bondtypes] …)
    --------------------------------------------------------
    atomtypes   : List[[name, mass, charge, sigma, epsilon]]
    bondtypes   : List[[a1_name, a2_name, funct, [params]]]
    angletypes  : List[[a1_name, a2_name, a3_name, funct, [params]]]
    torsiontypes: List[[a1_name, a2_name, a3_name, a4_name, funct, [params]]]

    Topology (from [moleculetype] / [atoms] / [bonds] … / [molecules])
    -------------------------------------------------------------------
    mol_types       : List[str]         unique molecule type names in order
    atomlist        : List[List[...]]   atomlist[i] for mol type i
                      each atom: [atom_name, type_name, charge, mass]
    bondlist        : List[List[...]]   bondlist[i] for mol type i
                      each bond: [atom1, atom2, type_index]   (1-based atoms)
    anglelist       : List[List[...]]
                      each angle: [atom1, atom2, atom3, type_index]
    torsionlist     : List[List[...]]
                      each torsion: [atom1, atom2, atom3, atom4,
                                     type_index, funct, improper, nparams]
    mol_instance_list: List[str]        flat list from [molecules] section

    Potential types derived from inline params (set if inline params present)
    -------------------------------------------------------------------------
    bond_types_from_mol    : List[[a1, a2, funct, [params]]]
    angle_types_from_mol   : List[[a1, a2, a3, funct, [params]]]
    torsion_types_from_mol : List[[a1, a2, a3, a4, funct, improper, [params]]]
    """
    # defaults
    comb_rule: int = 2
    fudge_lj: float = 0.5
    fudge_qq: float = 0.8333

    # global defs
    atomtypes: List = field(default_factory=list)
    bondtypes: List = field(default_factory=list)
    angletypes: List = field(default_factory=list)
    torsiontypes: List = field(default_factory=list)

    # per-molecule topology
    mol_types: List[str] = field(default_factory=list)
    atomlist: List = field(default_factory=list)
    bondlist: List = field(default_factory=list)
    anglelist: List = field(default_factory=list)
    torsionlist: List = field(default_factory=list)
    mol_instance_list: List[str] = field(default_factory=list)

    # inline-param-derived potential types
    bond_types_from_mol: List = field(default_factory=list)
    angle_types_from_mol: List = field(default_factory=list)
    torsion_types_from_mol: List = field(default_factory=list)


# ---------------------------------------------------------------------------
# Pure helper functions (module level, same logic as convert_gromacs_udf.py)
# ---------------------------------------------------------------------------

def _is_comment(line: str) -> bool:
    s = line.strip()
    return bool(s) and (s[0] == ";" or s[0] == "#")


def _is_blank(line: str) -> bool:
    return line.strip() == ""


def get_bond_name(n1: str, n2: str) -> str:
    if n1 <= n2:
        return n1 + "-" + n2
    return n2 + "-" + n1


def get_angle_name(n1: str, n2: str, n3: str) -> str:
    if n1 <= n3:
        return n1 + "-" + n2 + "-" + n3
    return n3 + "-" + n2 + "-" + n1


def get_torsion_name(n1: str, n2: str, n3: str, n4: str) -> str:
    if n1 < n4:
        return n1 + "-" + n2 + "-" + n3 + "-" + n4
    elif n1 == n4:
        if n2 <= n3:
            return n1 + "-" + n2 + "-" + n3 + "-" + n4
        return n4 + "-" + n3 + "-" + n2 + "-" + n1
    return n4 + "-" + n3 + "-" + n2 + "-" + n1


def is_improper(torsion_1based: List[int], bondlist_mol: List) -> bool:
    """Return True if the torsion is an improper dihedral (broken chain)."""
    a1, a2, a3, a4 = torsion_1based
    checklist = [(a1, a2), (a2, a3), (a3, a4)]
    for ca, cb in checklist:
        found = False
        for bond in bondlist_mol:
            if (ca == bond[0] and cb == bond[1]) or (ca == bond[1] and cb == bond[0]):
                found = True
                break
        if not found:
            return True
    return False


# ---------------------------------------------------------------------------
# TopParser class
# ---------------------------------------------------------------------------

class TopParser:
    """Parse a GROMACS .top file (with resolved #include ITP) into TopRawData."""

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def parse(self, top_path: str) -> TopRawData:
        """Parse *top_path* and return a :class:`TopRawData`."""
        lines = self._resolve_includes(top_path)
        raw = TopRawData()

        # Defaults
        raw.comb_rule, raw.fudge_lj, raw.fudge_qq = self._parse_defaults(lines)

        # Global type definitions (atomtypes, bondtypes, …)
        at, bt, agt, tt = self._parse_definitions(lines)
        raw.atomtypes = at
        raw.bondtypes = bt
        raw.angletypes = agt
        raw.torsiontypes = tt

        # Per-molecule topology + molecule instances
        (raw.mol_types,
         raw.atomlist,
         raw.bondlist,
         raw.anglelist,
         raw.torsionlist,
         raw.mol_instance_list,
         raw.bond_types_from_mol,
         raw.angle_types_from_mol,
         raw.torsion_types_from_mol) = self._parse_topology(lines)

        return raw

    # ------------------------------------------------------------------
    # Step 1: resolve #include
    # ------------------------------------------------------------------

    def _resolve_includes(self, top_path: str) -> List[str]:
        """Read *top_path* and inline all ``#include "*.itp"`` files."""
        dirpath = os.path.dirname(os.path.abspath(top_path))
        with open(top_path, "r") as f:
            lines = f.readlines()
        return self._inline_includes(lines, dirpath)

    def _inline_includes(self, lines: List[str], dirpath: str) -> List[str]:
        result: List[str] = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#include") and ".itp" in stripped:
                parts = stripped.split()
                if len(parts) >= 2:
                    itp_name = parts[1].strip('"\'')
                    itp_path = os.path.join(dirpath, itp_name)
                    try:
                        with open(itp_path, "r") as f:
                            itp_lines = f.readlines()
                        result.append("\n")
                        result.extend(itp_lines)
                        print("itp was included from", itp_path)
                        continue
                    except OSError as exc:
                        print("Warning: could not include {}: {}".format(itp_path, exc))
            result.append(line)
        return result

    # ------------------------------------------------------------------
    # Step 2: parse [ defaults ]
    # ------------------------------------------------------------------

    def _parse_defaults(self, lines: List[str]) -> Tuple[int, float, float]:
        """Return (comb_rule, fudge_lj, fudge_qq)."""
        in_section = False
        for line in lines:
            stripped = line.strip()
            if "[ defaults ]" in stripped:
                in_section = True
                continue
            if in_section:
                if _is_comment(line) or _is_blank(line):
                    continue
                parts = line.split()
                # nbfunc comb_rule gen-pairs fudgeLJ fudgeQQ
                comb_rule = int(parts[1])
                fudge_lj = float(parts[3])
                fudge_qq = float(parts[4])
                return comb_rule, fudge_lj, fudge_qq
        return 2, 0.5, 0.8333

    # ------------------------------------------------------------------
    # Step 3: parse global type definitions
    # ------------------------------------------------------------------

    def _parse_definitions(
        self, lines: List[str]
    ) -> Tuple[List, List, List, List]:
        """
        Parse [ atomtypes ], [ bondtypes ], [ angletypes ], [ dihedraltypes ].

        Returns (atomtypes, bondtypes, angletypes, torsiontypes).
        """
        flg_atom = False
        flg_bond = False
        flg_angle = False
        flg_torsion = False

        atomtypes: List = []
        bondtypes: List = []
        angletypes: List = []
        torsiontypes: List = []

        for line in lines:
            stripped = line.strip()

            if "[ atomtypes ]" in stripped:
                print(" reading atomtypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = True, False, False, False
                continue
            elif "[ bondtypes ]" in stripped:
                print(" reading bondtypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = False, True, False, False
                continue
            elif "[ angletypes ]" in stripped:
                print(" reading angletypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = False, False, True, False
                continue
            elif "[ dihedraltypes ]" in stripped:
                print(" reading torsiontypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = False, False, False, True
                continue
            elif stripped.startswith("[") and stripped.endswith("]"):
                # any other section → stop current
                flg_atom = flg_bond = flg_angle = flg_torsion = False
                continue

            if _is_comment(line) or _is_blank(line):
                continue

            parts = line.split()

            if flg_atom:
                # name [atnum] mass charge ptype sigma epsilon
                ndata = len(parts)
                if ";" in parts:
                    ndata = parts.index(";")
                name    = parts[0]
                mass    = float(parts[ndata - 5])
                charge  = float(parts[ndata - 4])
                ptype   = parts[ndata - 3]
                sigma   = float(parts[ndata - 2])
                epsilon = float(parts[ndata - 1])
                if ptype == "A":
                    atomtypes.append([name, mass, charge, sigma, epsilon])
                else:
                    print("Warning! particle type {} detected".format(ptype))

            elif flg_bond:
                atom1 = parts[0]
                atom2 = parts[1]
                funct = int(parts[2])
                params: List[float] = []
                for p in parts[3:]:
                    if p == ";":
                        break
                    params.append(float(p))
                bondtypes.append([atom1, atom2, funct, params])

            elif flg_angle:
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                funct = int(parts[3])
                params = []
                for p in parts[4:]:
                    if p == ";":
                        break
                    params.append(float(p))
                angletypes.append([atom1, atom2, atom3, funct, params])

            elif flg_torsion:
                atom1 = parts[0]
                atom2 = parts[1]
                atom3 = parts[2]
                atom4 = parts[3]
                funct = int(parts[4])
                params = []
                for p in parts[5:]:
                    if p == ";":
                        break
                    params.append(float(p))
                torsiontypes.append([atom1, atom2, atom3, atom4, funct, params])

        print("{} atoms, {} bonds, {} angles, {} dihedrals "
              "have been read (global defs)".format(
                  len(atomtypes), len(bondtypes),
                  len(angletypes), len(torsiontypes)))
        return atomtypes, bondtypes, angletypes, torsiontypes

    # ------------------------------------------------------------------
    # Step 4: parse per-molecule topology
    # ------------------------------------------------------------------

    def _parse_topology(self, lines: List[str]):  # noqa: C901
        """
        Parse [ moleculetype ], [ atoms ], [ bonds ], [ angles ],
        [ dihedrals ], and [ molecules ].

        Returns a 9-tuple matching convert_gromacs_udf.py read_top_data():
          (mol_types, atomlist, bondlist, anglelist, torsionlist,
           mol_instance_list,
           bond_types_from_mol, angle_types_from_mol, torsion_types_from_mol)
        """
        flg1 = False  # moleculetype
        flg2 = False  # atoms
        flg3 = False  # bonds
        flg4 = False  # angles
        flg5 = False  # dihedrals
        flg_mol = False  # molecules

        mol_types: List[str] = []
        atomlist: List = []
        bondlist: List = []
        anglelist: List = []
        torsionlist: List = []
        mol_instance_list: List[str] = []

        atomlist_mol: List = []
        bondlist_mol: List = []
        anglelist_mol: List = []
        torsionlist_mol: List = []

        bond_types_from_mol: List = []
        angle_types_from_mol: List = []
        torsion_types_from_mol: List = []

        bondtypes_map: Dict[str, int] = {}
        angletypes_map: Dict[str, int] = {}
        torsiontypes_map: Dict[str, int] = {}

        # torsion accumulation state
        atom1 = atom2 = atom3 = atom4 = 0
        funct = 0
        params: List[float] = []
        stmp: List[str] = []  # last parsed line tokens

        for line in lines:
            stripped = line.strip()

            # --- section headers ---
            if "[ moleculetype ]" in stripped:
                print(" reading moleculetype")
                if mol_types:
                    atomlist.append(atomlist_mol)
                    bondlist.append(bondlist_mol)
                    anglelist.append(anglelist_mol)
                    torsionlist.append(torsionlist_mol)
                    atomlist_mol = []
                    bondlist_mol = []
                    anglelist_mol = []
                    torsionlist_mol = []
                flg1, flg2, flg3, flg4, flg5 = True, False, False, False, False
                continue
            elif "[ atoms ]" in stripped:
                print(" reading atoms")
                flg1, flg2, flg3, flg4, flg5 = False, True, False, False, False
                continue
            elif "[ bonds ]" in stripped:
                print(" reading bonds")
                flg1, flg2, flg3, flg4, flg5 = False, False, True, False, False
                continue
            elif "[ angles ]" in stripped:
                print(" reading angles")
                flg1, flg2, flg3, flg4, flg5 = False, False, False, True, False
                continue
            elif "[ dihedrals ]" in stripped:
                print(" reading dihedrals")
                # flush pending torsion
                if len(stmp) == 8 and len(params) > 0:
                    imp = is_improper([atom1, atom2, atom3, atom4], bondlist_mol)
                    tidx = self._put_torsion_type(
                        atom1, atom2, atom3, atom4, funct, imp, params,
                        atomlist_mol, torsiontypes_map, torsion_types_from_mol)
                    torsionlist_mol.append([atom1, atom2, atom3, atom4,
                                            tidx, funct, imp, len(params)])
                flg1, flg2, flg3, flg4, flg5 = False, False, False, False, True
                atom1 = atom2 = atom3 = atom4 = 0
                params = []
                stmp = []
                continue
            elif "[ molecules ]" in stripped:
                print(" reading molecules")
                # flush pending torsion
                if flg5 and len(stmp) == 8 and len(params) > 0:
                    imp = is_improper([atom1, atom2, atom3, atom4], bondlist_mol)
                    tidx = self._put_torsion_type(
                        atom1, atom2, atom3, atom4, funct, imp, params,
                        atomlist_mol, torsiontypes_map, torsion_types_from_mol)
                    torsionlist_mol.append([atom1, atom2, atom3, atom4,
                                            tidx, funct, imp, len(params)])
                atomlist.append(atomlist_mol)
                bondlist.append(bondlist_mol)
                anglelist.append(anglelist_mol)
                torsionlist.append(torsionlist_mol)
                flg1 = flg2 = flg3 = flg4 = flg5 = False
                flg_mol = True
                continue
            elif stripped.startswith("[") and stripped.endswith("]"):
                # unknown section — stop current flags
                flg1 = flg2 = flg3 = flg4 = flg5 = False
                continue

            # --- data lines ---
            if _is_comment(line):
                continue
            if _is_blank(line):
                if flg1:
                    flg1 = False
                elif flg2:
                    flg2 = False
                elif flg3:
                    flg3 = False
                elif flg4:
                    flg4 = False
                elif flg5:
                    # blank line in dihedrals: flush pending torsion
                    if len(stmp) == 8 and len(params) > 0:
                        imp = is_improper([atom1, atom2, atom3, atom4], bondlist_mol)
                        tidx = self._put_torsion_type(
                            atom1, atom2, atom3, atom4, funct, imp, params,
                            atomlist_mol, torsiontypes_map, torsion_types_from_mol)
                        torsionlist_mol.append([atom1, atom2, atom3, atom4,
                                                tidx, funct, imp, len(params)])
                        atom1 = atom2 = atom3 = atom4 = 0
                        params = []
                        stmp = []
                    flg5 = False
                elif flg_mol:
                    flg_mol = False
                continue

            stmp = line.split()

            if flg1:
                mol_types.append(stmp[0])

            elif flg2:
                # nr type resnr residu atom cgnr charge [mass]
                atype  = stmp[1]
                atom_n = stmp[4]
                charge = float(stmp[6])
                mass   = float(stmp[7]) if len(stmp) > 7 else 0.0
                atomlist_mol.append([atom_n, atype, charge, mass])

            elif flg3:
                if len(stmp) >= 5:
                    # inline params: ai aj funct b0 kb
                    a1 = int(stmp[0])
                    a2 = int(stmp[1])
                    fn = int(stmp[2])
                    p1 = stmp[3]
                    p2 = stmp[4]
                    tidx = self._put_bond_type(
                        a1, a2, fn, p1, p2, atomlist_mol,
                        bondtypes_map, bond_types_from_mol)
                    bondlist_mol.append([a1, a2, tidx])
                else:
                    a1 = int(stmp[0])
                    a2 = int(stmp[1])
                    btype = int(stmp[2])
                    bondlist_mol.append([a1, a2, btype])

            elif flg4:
                if len(stmp) >= 6:
                    a1 = int(stmp[0])
                    a2 = int(stmp[1])
                    a3 = int(stmp[2])
                    fn = int(stmp[3])
                    p1 = stmp[4]
                    p2 = stmp[5]
                    tidx = self._put_angle_type(
                        a1, a2, a3, fn, p1, p2, atomlist_mol,
                        angletypes_map, angle_types_from_mol)
                    anglelist_mol.append([a1, a2, a3, tidx])
                else:
                    a1 = int(stmp[0])
                    a2 = int(stmp[1])
                    a3 = int(stmp[2])
                    ptype = int(stmp[3])
                    anglelist_mol.append([a1, a2, a3, ptype])

            elif flg5:
                if len(stmp) == 8:
                    # Amber-type: ai aj ak al funct phase k mult
                    t1 = int(stmp[0])
                    t2 = int(stmp[1])
                    t3 = int(stmp[2])
                    t4 = int(stmp[3])
                    tf = int(stmp[4])
                    if (t1 == atom1 and t2 == atom2
                            and t3 == atom3 and t4 == atom4
                            and tf == funct):
                        # multi-term: accumulate
                        params += [float(stmp[5]), float(stmp[6]),
                                   int(stmp[7])]
                        continue
                    # flush previous
                    if len(params) > 0:
                        imp = is_improper([atom1, atom2, atom3, atom4],
                                          bondlist_mol)
                        tidx = self._put_torsion_type(
                            atom1, atom2, atom3, atom4, funct, imp, params,
                            atomlist_mol, torsiontypes_map,
                            torsion_types_from_mol)
                        torsionlist_mol.append([atom1, atom2, atom3, atom4,
                                                tidx, funct, imp, len(params)])
                    atom1, atom2, atom3, atom4, funct = t1, t2, t3, t4, tf
                    params = [float(stmp[5]), float(stmp[6]), int(stmp[7])]

                elif len(stmp) == 11:
                    # Cosine_Polynomial: ai aj ak al funct C0..C5
                    atom1 = int(stmp[0])
                    atom2 = int(stmp[1])
                    atom3 = int(stmp[2])
                    atom4 = int(stmp[3])
                    funct = int(stmp[4])
                    ps = [float(stmp[i]) for i in range(5, 11)]
                    imp = is_improper([atom1, atom2, atom3, atom4], bondlist_mol)
                    tidx = self._put_torsion_type(
                        atom1, atom2, atom3, atom4, funct, imp, ps,
                        atomlist_mol, torsiontypes_map, torsion_types_from_mol)
                    torsionlist_mol.append([atom1, atom2, atom3, atom4,
                                            tidx, funct, imp, len(ps)])
                    params = []
                    stmp = []

                else:
                    # type-reference only
                    atom1 = int(stmp[0])
                    atom2 = int(stmp[1])
                    atom3 = int(stmp[2])
                    atom4 = int(stmp[3])
                    ptype = int(stmp[4])
                    torsionlist_mol.append([atom1, atom2, atom3, atom4, ptype])

            elif flg_mol:
                moltype = stmp[0]
                nmols   = int(stmp[1])
                for _ in range(nmols):
                    mol_instance_list.append(moltype)

        return (mol_types, atomlist, bondlist, anglelist, torsionlist,
                mol_instance_list,
                bond_types_from_mol, angle_types_from_mol,
                torsion_types_from_mol)

    # ------------------------------------------------------------------
    # Internal helpers: put*Type — deduplicate and register type entries
    # ------------------------------------------------------------------

    @staticmethod
    def _put_bond_type(atom1, atom2, funct, para1, para2,
                       atomlist_mol, bondtypes_map, bondtypes):
        a1 = atomlist_mol[atom1 - 1][1]
        a2 = atomlist_mol[atom2 - 1][1]
        key = (get_bond_name(a1, a2) + ","
               + str(funct) + "," + str(para1) + "," + str(para2))
        if key not in bondtypes_map:
            bondtypes_map[key] = len(bondtypes_map)
            bondtypes.append([a1, a2, funct, [float(para1), float(para2)]])
        return bondtypes_map[key]

    @staticmethod
    def _put_angle_type(atom1, atom2, atom3, funct, para1, para2,
                        atomlist_mol, angletypes_map, angletypes):
        a1 = atomlist_mol[atom1 - 1][1]
        a2 = atomlist_mol[atom2 - 1][1]
        a3 = atomlist_mol[atom3 - 1][1]
        key = (get_angle_name(a1, a2, a3) + ","
               + str(funct) + "," + str(para1) + "," + str(para2))
        if key not in angletypes_map:
            angletypes_map[key] = len(angletypes_map)
            angletypes.append([a1, a2, a3, funct,
                                [float(para1), float(para2)]])
        return angletypes_map[key]

    @staticmethod
    def _put_torsion_type(atom1, atom2, atom3, atom4, funct, improper, params,
                          atomlist_mol, torsiontypes_map, torsiontypes):
        a1 = atomlist_mol[atom1 - 1][1]
        a2 = atomlist_mol[atom2 - 1][1]
        a3 = atomlist_mol[atom3 - 1][1]
        a4 = atomlist_mol[atom4 - 1][1]
        key = (get_torsion_name(a1, a2, a3, a4)
               + "," + str(improper) + "," + str(params))
        if key not in torsiontypes_map:
            torsiontypes_map[key] = len(torsiontypes_map)
            torsiontypes.append([a1, a2, a3, a4, funct, improper, params])
        return torsiontypes_map[key]
