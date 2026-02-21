# -*- coding: utf-8 -*-
"""
udf_adapter.py
--------------
Adapter layer: reads a COGNAC UDF file via UDFManager and builds a SystemModel.

The public entry point is::

    model = UdfAdapter(udf_manager).build()

All UDFManager calls are contained here; no UDF access occurs in the Writers.
"""
from __future__ import annotations

import math
import os
from typing import List, Optional, Tuple

from ..core.system_model import (
    AtomType, AtomRecord, BondRecord, PairRecord, AngleRecord, DihedralRecord,
    MoleculeTopology, AtomPosition, CellGeometry, NdxData,
    SimulationParams, SystemModel,
)

# Charge conversion constant (same as original script)
_CHARGE_UNIT = 18.224159264


# ---------------------------------------------------------------------------
# Helpers (ported from udf2gro.py, but now pure functions / instance methods)
# ---------------------------------------------------------------------------

def _shorten_molname(ss: str) -> str:
    if len(ss) > 5:
        ss = ss[0] + ss[1] + ss[len(ss)-3] + ss[len(ss)-2] + ss[len(ss)-1]
    return ss


def _float2str(value: float, ndigit: int) -> str:
    return "{:.15f}".format(value)[:ndigit]


def _is_rectangular(cell_raw, thres: float = 1e-5) -> bool:
    for check in [cell_raw[3], cell_raw[4], cell_raw[5]]:
        if abs(check - 90.0) > thres:
            return False
    return True


def _get_proper_dihedral_params(udf, j: int, n: int, kk: float):
    """Port of get_proper_dihedral_params() from udf2gro.py."""
    loc = "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.p[]"
    gro_k = 0.0
    gro_phi = 0.0
    gro_mult = 0
    if n == 1:
        gro_mult = 0
        A0 = udf.get(loc, [j, 0])
        gro_phi = 0.0 if A0 != 0 else 180.0
        gro_k = A0 / 2 * kk
    elif n == 2:
        gro_mult = 1
        A0 = udf.get(loc, [j, 0])
        A1 = udf.get(loc, [j, 1])
        A01 = A0 * A1
        gro_phi = 0.0 if A01 < 0 else 180.0
        gro_k = A0 * kk
    elif n == 3:
        gro_mult = 2
        A0 = udf.get(loc, [j, 0])
        A2 = udf.get(loc, [j, 2])
        if A0 == 0:
            gro_phi = 0.0
            gro_k = A2 / 2 * kk
        else:
            gro_phi = 180.0
            gro_k = A0 / 2 * kk
    elif n == 4:
        gro_mult = 3
        A0 = udf.get(loc, [j, 0])
        A1 = udf.get(loc, [j, 1])
        A01 = A0 * A1
        gro_phi = 0.0 if A01 >= 0 else 180.0
        gro_k = A0 * kk
    elif n == 5:
        gro_mult = 4
        A0 = udf.get(loc, [j, 0])
        A2 = udf.get(loc, [j, 2])
        if A0 != 0:
            gro_phi = 0.0
            gro_k = A0 / 2 * kk
        else:
            gro_phi = 180.0
            gro_k = A2 / 8 * kk
    elif n == 6:
        gro_mult = 5
        A0 = udf.get(loc, [j, 0])
        A1 = udf.get(loc, [j, 1])
        A01 = A0 * A1
        gro_phi = 0.0 if A01 < 0 else 180.0
        gro_k = A0 * kk
    elif n == 7:
        gro_mult = 6
        A0 = udf.get(loc, [j, 0])
        A2 = udf.get(loc, [j, 2])
        if A0 == 0:
            gro_phi = 0.0
            gro_k = A2 / 18 * kk
        else:
            gro_phi = 180.0
            gro_k = A0 / 2 * kk
    return gro_phi, gro_k, gro_mult


# ---------------------------------------------------------------------------
# Main adapter class
# ---------------------------------------------------------------------------

class UdfAdapter:
    """Reads a UDFManager object and produces a SystemModel."""

    def __init__(self, udf):
        self._udf = udf

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def build(self) -> SystemModel:
        udf = self._udf
        udf.jump(udf.totalRecord() - 1)
        print(" Data output : Record number = ", udf.totalRecord() - 1)

        udf_path = str(udf.udfFile()).strip()
        title = udf_path.split("/")[-1] if "/" in udf_path else udf_path

        calcQQ = udf.get("Simulation_Conditions.Calc_Potential_Flags.Electrostatic")
        print(" Electrostatic is {}.".format("ON" if calcQQ == 1 else "OFF"))

        # --- counts ---
        atm_type_num = udf.size(udf.rlocation("Molecular_Attributes.Atom_Type[]", []))
        mol_num      = udf.size(udf.rlocation("Set_of_Molecules.molecule[]",       []))
        int_num      = udf.size(udf.rlocation("Interactions.Pair_Interaction[]",   []))
        bnd_type_num = udf.size(udf.rlocation("Molecular_Attributes.Bond_Potential[]",     []))
        agl_type_num = udf.size(udf.rlocation("Molecular_Attributes.Angle_Potential[]",    []))
        tor_type_num = udf.size(udf.rlocation("Molecular_Attributes.Torsion_Potential[]",  []))

        all_atm_num = sum(
            udf.size(udf.rlocation("Set_of_Molecules.molecule[].atom[]", [i]))
            for i in range(mol_num)
        )

        # --- OPLS detection ---
        totallyOPLS, partiallyOPLS = self._detect_opls()
        comb_rule = 3 if totallyOPLS else 2
        if totallyOPLS:
            print("OPLS FF is detected. LJ mixing-rule -> geometric")
        elif partiallyOPLS:
            print("Warning!! LJ mixing-rule cannot be identified!")

        # --- fudge factors ---
        flags14 = udf.get("Simulation_Conditions.Calc_Potential_Flags.Non_Bonding_1_4")
        fudgeLJ = udf.getArray("Interactions.Pair_Interaction[].Scale_1_4_Pair", [0])
        fudgeQQ = 0.0
        if calcQQ == 1:
            fudgeQQ = udf.getArray("Interactions.Electrostatic_Interaction[].Scale_1_4_Pair", [0])

        # --- atomname_in_gro map ---
        atomname_in_gro = self._build_atomname_map()

        # --- mol name mappings ---
        mol_name_list, molname_map, mol_name = self._build_mol_name_mappings(mol_num)

        mol_type_num = len(mol_name)
        print(" Molecular type's num = ", mol_type_num)
        print(" Molecular's num      = ", mol_num)
        print(" Atom's num   = ", all_atm_num)
        print(" Atomtype num = ", atm_type_num)
        print(" molecule name map (udf -> top)")
        for name in molname_map.keys():
            print(" {:15s} -> {:5s}".format(name, molname_map[name]))

        # --- bond name -> potential type map ---
        bond_name_potmap = self._build_bond_name_potmap()

        # --- atom types with LJ params ---
        atom_types, lj_cutoff = self._build_atom_types(
            atm_type_num, int_num, totallyOPLS
        )

        # --- mol topologies ---
        mol_topologies = self._build_mol_topologies(
            mol_type_num, mol_name, mol_name_list, mol_num, molname_map,
            atomname_in_gro, calcQQ, bond_name_potmap,
            bnd_type_num, agl_type_num, tor_type_num,
            totallyOPLS, int_num
        )

        # --- mol sequence for [ molecules ] ---
        mol_sequence = self._build_mol_sequence(mol_num, mol_name_list, molname_map)

        # --- atom positions (gro structure) ---
        atom_positions, cell, vel_gen, udf_gro = self._build_atom_positions(
            mol_num, atomname_in_gro
        )

        # --- NDX data (optional) ---
        ndx_data = self._build_ndx_data(mol_num, mol_name_list, molname_map)

        # --- simulation params ---
        sim_params = self._build_sim_params(
            title, calcQQ, lj_cutoff, all_atm_num,
            vel_gen, cell, udf_gro
        )

        return SystemModel(
            title=title,
            udf_path=udf_path,
            comb_rule=comb_rule,
            flags14=flags14,
            fudgeLJ=fudgeLJ,
            fudgeQQ=fudgeQQ,
            calcQQ=calcQQ,
            atom_types=atom_types,
            mol_topologies=mol_topologies,
            mol_sequence=mol_sequence,
            atom_positions=atom_positions,
            cell=cell,
            sim_params=sim_params,
            ndx_data=ndx_data,
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _detect_opls(self):
        atypes = self._udf.get("Molecular_Attributes.Atom_Type[].Name")
        totallyOPLS = True
        partiallyOPLS = False
        for atype in atypes:
            if '&' in atype:
                partiallyOPLS = True
            else:
                totallyOPLS = False
        return totallyOPLS, partiallyOPLS

    def _build_atomname_map(self):
        udf = self._udf
        atomtypedata = udf.get("Molecular_Attributes.Atom_Type[].Name")
        atomname_in_gro = {}
        for i, name in enumerate(atomtypedata):
            atomname_in_gro[name] = name[:2] + hex(i)[2:]
        return atomname_in_gro

    def _build_mol_name_mappings(self, mol_num: int):
        udf = self._udf
        mol_name_list = list(udf.get("Set_of_Molecules.molecule[].Mol_Name"))
        map_molname_atypes = {}
        ncount_diff = 1

        for j, molname in enumerate(mol_name_list):
            atypelist = udf.get(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [j]
            )
            if molname in map_molname_atypes:
                if atypelist == map_molname_atypes[molname]:
                    pass
                else:
                    molname_new = self._search_molname_same_topol(
                        atypelist, map_molname_atypes
                    )
                    if molname_new is None:
                        molname_new = "".join(a[0] for a in atypelist)
                        map_molname_atypes[molname_new] = atypelist
                        ncount_diff += 1
                    print("mol[{}], name={} is renamed to {}".format(
                        j, molname, molname_new
                    ))
                    mol_name_list[j] = molname_new
            else:
                map_molname_atypes[molname] = atypelist

        molname_map = {}
        mol_name = []
        for j, molname in enumerate(mol_name_list):
            if molname in mol_name:
                continue
            molname_gromacs = "M{:04d}".format(j)
            molname_map[molname] = molname_gromacs
            mol_name.append(molname)

        return mol_name_list, molname_map, mol_name

    @staticmethod
    def _search_molname_same_topol(atypes, map_name_atypes):
        for name, atypes_check in map_name_atypes.items():
            if atypes == atypes_check:
                return name
        return None

    def _build_bond_name_potmap(self):
        udf = self._udf
        location = "Molecular_Attributes.Bond_Potential[]"
        ntypes = udf.size(location)
        result = {}
        for itype in range(ntypes):
            name    = udf.get(location + ".Name",           [itype])
            pottype = udf.get(location + ".Potential_Type", [itype])
            result[name] = pottype
        return result

    def _build_atom_types(self, atm_type_num: int, int_num: int,
                          totallyOPLS: bool):
        udf = self._udf
        atom_types = []
        lj_cutoff = 0.0

        for i in range(atm_type_num):
            name = udf.getArray("Molecular_Attributes.Atom_Type[].Name", [i])
            mass = udf.getArray("Molecular_Attributes.Atom_Type[].Mass", [i])

            sigma   = 0.0
            epsilon = 0.0
            found   = False
            for j in range(int_num):
                site1   = udf.getArray("Interactions.Pair_Interaction[].Site1_Name", [j])
                site2   = udf.getArray("Interactions.Pair_Interaction[].Site2_Name", [j])
                pottype = udf.get("Interactions.Pair_Interaction[].Potential_Type", [j])
                if pottype != "Lennard_Jones":
                    print("Error!! non LJ(12,6) type potential is detected")
                    raise RuntimeError("non-LJ potential")
                if site1 == name and site1 == site2:
                    sigma   = udf.get("Interactions.Pair_Interaction[].Lennard_Jones.sigma",   [j], "[nm]")
                    epsilon = udf.get("Interactions.Pair_Interaction[].Lennard_Jones.epsilon",  [j], "[kJ/mol]")
                    cc = udf.get("Interactions.Pair_Interaction[].Cutoff", [j], "[nm]")
                    if cc > lj_cutoff:
                        lj_cutoff = cc
                    found = True
                    break

            atom_types.append(AtomType(
                name=name, mass=mass, sigma=sigma, epsilon=epsilon
            ))

        return atom_types, lj_cutoff

    # ------------------------------------------------------------------
    # Molecule topology building
    # ------------------------------------------------------------------

    def _build_mol_topologies(
        self, mol_type_num, mol_name, mol_name_list, mol_num, molname_map,
        atomname_in_gro, calcQQ, bond_name_potmap,
        bnd_type_num, agl_type_num, tor_type_num,
        totallyOPLS, int_num
    ):
        udf = self._udf
        topologies = []

        for m in range(mol_type_num):
            udf_name = mol_name[m]
            gro_raw  = molname_map[udf_name]
            gro_name = _shorten_molname(gro_raw)
            print(" moltype %d %s" % (m, gro_name))

            # find representative molecule index
            mol_no = next(
                i for i in range(mol_num)
                if mol_name_list[i] == udf_name
            )

            topo = MoleculeTopology(
                udf_name=udf_name,
                gro_name=gro_name,
                nrexcl=3,
            )

            # --- atoms ---
            topo.atoms = self._build_atoms(
                mol_no, gro_name, atomname_in_gro, calcQQ
            )

            # --- bonds ---
            topo.bonds = self._build_bonds(
                mol_no, bond_name_potmap, bnd_type_num
            )

            # --- pairs (1-4) ---
            topo.pairs = self._build_pairs(mol_no, int_num)

            # --- angles ---
            topo.angles = self._build_angles(mol_no, agl_type_num)

            # --- dihedrals ---
            topo.dihedrals = self._build_dihedrals(mol_no, tor_type_num)

            topologies.append(topo)

        return topologies

    def _build_atoms(self, mol_no, gro_name, atomname_in_gro, calcQQ):
        udf = self._udf
        atm_num = udf.size("Set_of_Molecules.molecule[].atom[]", [mol_no])
        ele_num = udf.size("Set_of_Molecules.molecule[].electrostatic_Site[]", [mol_no])
        print("  atom_num     = %d" % atm_num)

        atoms = []
        for i in range(atm_num):
            atype = udf.get(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [mol_no, i]
            )
            gro_aname = atomname_in_gro[atype]

            if ele_num == 0 or calcQQ == 0:
                charge = 0.0
            else:
                q_raw = udf.get(
                    "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element",
                    [mol_no, i]
                )
                charge = round(q_raw / _CHARGE_UNIT, 8)

            atoms.append(AtomRecord(
                index=i + 1,
                type_name=atype,
                gro_name=gro_aname,
                charge=charge,
            ))
        return atoms

    def _build_bonds(self, mol_no, bond_name_potmap, bnd_type_num):
        udf = self._udf
        bnd_num = udf.size(udf.rlocation("Set_of_Molecules.molecule[].bond[]", [mol_no]))
        print("  bond_num     = %d" % bnd_num)

        bonds = []
        for i in range(bnd_num):
            a1 = udf.getArray("Set_of_Molecules.molecule[].bond[].atom1", [mol_no, i]) + 1
            a2 = udf.getArray("Set_of_Molecules.molecule[].bond[].atom2", [mol_no, i]) + 1
            pot_name = udf.getArray(
                "Set_of_Molecules.molecule[].bond[].Potential_Name", [mol_no, i]
            )
            pottype = bond_name_potmap[pot_name]

            funct = "1"
            r0 = 0.0
            kb = 0.0
            for j in range(bnd_type_num):
                name = str(udf.getArray("Molecular_Attributes.Bond_Potential[].Name", [j]))
                if name == pot_name:
                    if pottype == "FENE_LJ":
                        funct = "7"
                        r0 = udf.get(
                            "Molecular_Attributes.Bond_Potential[].FENE_LJ.R_max",
                            [j], "[sigma]"
                        )
                        kb = udf.get(
                            "Molecular_Attributes.Bond_Potential[].FENE_LJ.K",
                            [j], "[epsilon/(sigma^2)]"
                        )
                    elif pottype == "Harmonic":
                        funct = "1"
                        r0 = udf.get("Molecular_Attributes.Bond_Potential[].R0",
                                     [j], "[nm]")
                        kb = udf.get("Molecular_Attributes.Bond_Potential[].Harmonic.K",
                                     [j], "[kJ/mol/nm^2]")
                    break

            bonds.append(BondRecord(atom1=a1, atom2=a2, funct=funct, r0=r0, kb=kb))
        return bonds

    def _build_pairs(self, mol_no: int, int_num: int) -> List[PairRecord]:
        """Build 1-4 pair list using the same algorithm as get14() in original."""
        udf = self._udf
        atom_num = udf.size("Set_of_Molecules.molecule[].atom[]", [mol_no])

        bond_map  = [[] for _ in range(atom_num)]
        angle_map = [[] for _ in range(atom_num)]

        bond_num = udf.size("Set_of_Molecules.molecule[].bond[]", [mol_no])
        for j in range(bond_num):
            a1 = udf.get("Set_of_Molecules.molecule[].bond[].atom1", [mol_no, j])
            a2 = udf.get("Set_of_Molecules.molecule[].bond[].atom2", [mol_no, j])
            bond_map[a1].append(a2)
            bond_map[a2].append(a1)

        angle_num = udf.size("Set_of_Molecules.molecule[].angle[]", [mol_no])
        for j in range(angle_num):
            a1 = udf.get("Set_of_Molecules.molecule[].angle[].atom1", [mol_no, j])
            a3 = udf.get("Set_of_Molecules.molecule[].angle[].atom3", [mol_no, j])
            angle_map[a1].append(a3)
            angle_map[a3].append(a1)

        pair1: List[int] = []
        pair2: List[int] = []

        tors_num = udf.size("Set_of_Molecules.molecule[].torsion[]", [mol_no])
        for j in range(tors_num):
            a1 = udf.get("Set_of_Molecules.molecule[].torsion[].atom1", [mol_no, j])
            a3 = udf.get("Set_of_Molecules.molecule[].torsion[].atom3", [mol_no, j])
            a4 = udf.get("Set_of_Molecules.molecule[].torsion[].atom4", [mol_no, j])
            if a3 not in bond_map[a1] and a4 not in angle_map[a1]:
                chk = True
                for m in range(len(pair1)):
                    if ((pair1[m] == a1 and pair2[m] == a4) or
                            (pair2[m] == a1 and pair1[m] == a4)):
                        chk = False
                        break
                if chk:
                    pair1.append(a1)
                    pair2.append(a4)

        del bond_map
        print("  1-4 pair_num =", len(pair1))

        site1names = udf.get("Interactions.Pair_Interaction[].Site1_Name")
        site2names = udf.get("Interactions.Pair_Interaction[].Site2_Name")

        pairs = []
        for i in range(len(pair1)):
            name1 = udf.getArray(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",
                [mol_no, pair1[i]]
            )
            name2 = udf.getArray(
                "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",
                [mol_no, pair2[i]]
            )
            found = False
            for j in range(int_num):
                s1 = site1names[j]
                s2 = site2names[j]
                if (s1 == name1 and s2 == name2) or (s1 == name2 and s2 == name1):
                    found = True
                    break
            # Always add pair (found or not) — matches original behaviour
            pairs.append(PairRecord(atom1=pair1[i] + 1, atom2=pair2[i] + 1))

        return pairs

    def _build_angles(self, mol_no: int, agl_type_num: int) -> List[AngleRecord]:
        udf = self._udf
        agl_num = udf.size(udf.rlocation("Set_of_Molecules.molecule[].angle[]", [mol_no]))
        print("  angle_num    = %d" % agl_num)

        angles = []
        for i in range(agl_num):
            a1 = udf.getArray("Set_of_Molecules.molecule[].angle[].atom1", [mol_no, i]) + 1
            a2 = udf.getArray("Set_of_Molecules.molecule[].angle[].atom2", [mol_no, i]) + 1
            a3 = udf.getArray("Set_of_Molecules.molecule[].angle[].atom3", [mol_no, i]) + 1
            pot_name = udf.getArray(
                "Set_of_Molecules.molecule[].angle[].Potential_Name", [mol_no, i]
            )
            theta0 = 0.0
            k = 0.0
            for j in range(agl_type_num):
                name = str(udf.getArray("Molecular_Attributes.Angle_Potential[].Name", [j]))
                if name == pot_name:
                    theta0 = 180.0 - float(udf.get(
                        "Molecular_Attributes.Angle_Potential[].theta0", [j], "[degree]"
                    ))
                    k = udf.get(
                        "Molecular_Attributes.Angle_Potential[].Theta.K",
                        [j], "[kJ/mol/rad^2]"
                    )
                    break
            angles.append(AngleRecord(atom1=a1, atom2=a2, atom3=a3,
                                       theta0=theta0, k=k))
        return angles

    def _build_dihedrals(self, mol_no: int, tor_type_num: int) -> List[DihedralRecord]:
        udf = self._udf
        tor_num = udf.size(udf.rlocation("Set_of_Molecules.molecule[].torsion[]", [mol_no]))
        print("  torsion_num  = %d" % tor_num)

        tors_pot_type_list = udf.get("Molecular_Attributes.Torsion_Potential[].Potential_Type")
        use_amber_header = "Amber" in tors_pot_type_list

        dihedrals = []
        for i in range(tor_num):
            a1 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom1", [mol_no, i]) + 1
            a2 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom2", [mol_no, i]) + 1
            a3 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom3", [mol_no, i]) + 1
            a4 = udf.getArray("Set_of_Molecules.molecule[].torsion[].atom4", [mol_no, i]) + 1
            pot_name = udf.getArray(
                "Set_of_Molecules.molecule[].torsion[].Potential_Name", [mol_no, i]
            )

            funct = ""
            params: List[float] = []

            for j in range(tor_type_num):
                name    = str(udf.getArray("Molecular_Attributes.Torsion_Potential[].Name",          [j]))
                pottype = str(udf.getArray("Molecular_Attributes.Torsion_Potential[].Potential_Type", [j]))
                if name != pot_name:
                    continue

                if pottype == "Cosine_Polynomial":
                    kk = float(udf.get(
                        "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.K",
                        [j], "[kJ/mol]"
                    ))
                    n = udf.getArray(
                        "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.N", [j]
                    )
                    if "oopa" in name:
                        funct = "4"
                        gro_phi, gro_k, gro_mult = _get_proper_dihedral_params(udf, j, n, kk)
                        params = [gro_phi, gro_k, gro_mult]
                    elif n <= 6:
                        funct = "3"
                        raw_params = udf.get(
                            "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.p[]",
                            [j]
                        )
                        for k_idx in range(6):
                            if k_idx < len(raw_params):
                                params.append(kk * float(raw_params[k_idx]))
                            else:
                                params.append(0.0)
                    else:
                        funct = "1"
                        gro_phi, gro_k, gro_mult = _get_proper_dihedral_params(udf, j, n, kk)
                        params = [gro_phi, gro_k, gro_mult]
                    break

                elif pottype == "Amber":
                    if ":" in name:
                        funct = "9"
                    else:
                        funct = "1"
                    pk    = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.PK",    [j], "[kJ/mol]")
                    idivf = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.IDIVF", [j])
                    pn    = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.PN",    [j])
                    phase = udf.get("Molecular_Attributes.Torsion_Potential[].Amber.PHASE", [j])
                    k_val = pk / idivf
                    params = [phase, k_val, pn]
                    break

            dihedrals.append(DihedralRecord(
                atom1=a1, atom2=a2, atom3=a3, atom4=a4,
                funct=funct, params=params
            ))
        return dihedrals

    # ------------------------------------------------------------------
    # Mol sequence ([ molecules ] section)
    # ------------------------------------------------------------------

    def _build_mol_sequence(self, mol_num, mol_name_list, molname_map):
        sequence = []
        if mol_num > 1:
            ncount = 1
            for i in range(1, mol_num):
                if mol_name_list[i] == mol_name_list[i - 1]:
                    ncount += 1
                    if i == mol_num - 1:
                        gn = _shorten_molname(molname_map[mol_name_list[i]])
                        sequence.append((gn, ncount))
                else:
                    gn = _shorten_molname(molname_map[mol_name_list[i - 1]])
                    sequence.append((gn, ncount))
                    ncount = 1
                    if i == mol_num - 1:
                        gn = _shorten_molname(molname_map[mol_name_list[i]])
                        sequence.append((gn, ncount))
        else:
            gn = _shorten_molname(molname_map[mol_name_list[0]])
            sequence.append((gn, 1))
        return sequence

    # ------------------------------------------------------------------
    # Atom positions (GRO structure)
    # ------------------------------------------------------------------

    def _build_atom_positions(self, mol_num, atomname_in_gro):
        udf = self._udf
        vel_gen = False
        udf_gro = udf

        if udf.get("Initial_Structure.Generate_Method.Method") == "Restart":
            rest_udfname = udf.get("Initial_Structure.Generate_Method.Restart.UDF_Name")
            rest_record  = udf.get("Initial_Structure.Generate_Method.Restart.Record")
            if len(rest_udfname) > 0:
                if len(os.path.dirname(rest_udfname)) == 0:
                    udfdir = udf.udfDirectory()
                    rest_udfpath = os.path.join(udfdir, rest_udfname)
                else:
                    rest_udfpath = rest_udfname
                print("Restart geometry is read from record {} in '{}'".format(
                    rest_record, rest_udfname
                ))
                if not os.path.exists(rest_udfpath):
                    raise RuntimeError(
                        "Error! UDF file for restart does not exist!\n"
                        "       Check UDF path 'Initial_Structure.Generate_Method.Restart.UDF_Name'"
                    )
                from UDFManager import UDFManager
                udf_gro = UDFManager(rest_udfpath)
                nrecs = udf_gro.totalRecord()
                if rest_record == -1:
                    udf_gro.jump(nrecs - 1)
                else:
                    udf_gro.jump(rest_record)

        positions = []
        for i in range(mol_num):
            atm_num = udf_gro.size(udf_gro.rlocation("Structure.Position.mol[].atom[]", [i]))
            mol_id = min(i + 1, 99999)
            ss = udf_gro.get("Set_of_Molecules.molecule[].Mol_Name", [i])
            if len(ss) > 5:
                ss = ss[0] + ss[1] + ss[2] + ss[len(ss)-2] + ss[len(ss)-1]

            for j in range(atm_num):
                x  = udf_gro.get("Structure.Position.mol[].atom[].x",  [i, j], "[nm]")
                y  = udf_gro.get("Structure.Position.mol[].atom[].y",  [i, j], "[nm]")
                z  = udf_gro.get("Structure.Position.mol[].atom[].z",  [i, j], "[nm]")
                vx = udf_gro.get("Structure.Velocity.mol[].atom[].x",  [i, j], "[nm/ps]")
                vy = udf_gro.get("Structure.Velocity.mol[].atom[].y",  [i, j], "[nm/ps]")
                vz = udf_gro.get("Structure.Velocity.mol[].atom[].z",  [i, j], "[nm/ps]")
                if vx is None or vy is None or vz is None:
                    vel_gen = True
                    vx = vy = vz = 0.0

                atom_type_name = udf_gro.get(
                    "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [i, j]
                )
                atom_gro_name = atomname_in_gro[atom_type_name]
                atom_id = min(
                    udf_gro.get("Set_of_Molecules.molecule[].atom[].Atom_ID", [i, j]) + 1,
                    99999
                )

                positions.append(AtomPosition(
                    mol_id=mol_id,
                    mol_name_short=ss,
                    atom_gro_name=atom_gro_name,
                    atom_id=atom_id,
                    x=x, y=y, z=z,
                    vx=vx, vy=vy, vz=vz,
                ))

        # Cell geometry
        cell_raw = udf_gro.get("Structure.Unit_Cell.Cell_Size")
        if _is_rectangular(cell_raw):
            a = udf_gro.get("Structure.Unit_Cell.Cell_Size.a", "[nm]")
            b = udf_gro.get("Structure.Unit_Cell.Cell_Size.b", "[nm]")
            c = udf_gro.get("Structure.Unit_Cell.Cell_Size.c", "[nm]")
            cell = CellGeometry(a=a, b=b, c=c)
        else:
            cell = CellGeometry(
                a=cell_raw[0] * 0.1,
                b=cell_raw[1] * 0.1,
                c=cell_raw[2] * 0.1,
                alpha=cell_raw[3],
                beta=cell_raw[4],
                gamma=cell_raw[5],
            )

        return positions, cell, vel_gen, udf_gro

    # ------------------------------------------------------------------
    # NDX data
    # ------------------------------------------------------------------

    def _build_ndx_data(self, mol_num, mol_name_list, molname_map) -> Optional[NdxData]:
        udf = self._udf
        loc_constr = "Simulation_Conditions.Constraint_Conditions.Constraint_Atom[]"
        ndata_constr = udf.size(loc_constr)
        if ndata_constr == 0:
            return None

        constr_axis_old = udf.get(loc_constr + ".Constraint_Axis", [0])
        constr_aid_list = []
        constr_axis = None

        for idat in range(ndata_constr):
            idx_mol_atom  = udf.get(loc_constr + ".Index",           [idat])
            constr_axis   = udf.get(loc_constr + ".Constraint_Axis", [idat])
            constr_method = udf.get(loc_constr + ".Method",          [idat])
            if constr_method != "Steady":
                raise RuntimeError(
                    "Error! constraint method {} is not supported!".format(constr_method)
                )
            constr_veloc = udf.get(loc_constr + ".Steady.Velocity", [idat])
            molidx  = idx_mol_atom[0]
            atomidx = idx_mol_atom[1]
            atom_id = udf.get(
                "Set_of_Molecules.molecule[].atom[].Atom_ID", [molidx, atomidx]
            ) + 1
            if constr_veloc[0] == 0.0 and constr_veloc[1] == 0.0 and constr_veloc[2] == 0.0:
                if constr_axis == constr_axis_old:
                    constr_axis_old = constr_axis
                    constr_aid_list.append(atom_id)
                else:
                    raise RuntimeError("Error! Not supported constraint conditions!")
            else:
                raise RuntimeError("Error! constraint velocity other than 0 is not supported!")

        atm_num_last = udf.size("Structure.Position.mol[].atom[]", [mol_num - 1])
        atom_id_max = udf.get(
            "Set_of_Molecules.molecule[].atom[].Atom_ID",
            [mol_num - 1, atm_num_last - 1]
        ) + 1

        molnames = set(mol_name_list)
        mol_groups = {}
        for molname in molnames:
            gn = _shorten_molname(molname_map[molname])
            aid_list = []
            for imol in range(mol_num):
                if udf.get("Set_of_Molecules.molecule[].Mol_Name", [imol]) == molname:
                    natoms = udf.size("Structure.Position.mol[].atom[]", [imol])
                    for iatom in range(natoms):
                        aid = udf.get(
                            "Set_of_Molecules.molecule[].atom[].Atom_ID", [imol, iatom]
                        ) + 1
                        aid_list.append(aid)
            mol_groups[gn] = aid_list

        return NdxData(
            atom_id_max=atom_id_max,
            mol_groups=mol_groups,
            constraint_atom_ids=constr_aid_list,
            constr_axis=constr_axis,
        )

    # ------------------------------------------------------------------
    # Simulation parameters
    # ------------------------------------------------------------------

    def _build_sim_params(
        self, title, calcQQ, lj_cutoff, all_atm_num,
        vel_gen, cell, udf_gro
    ) -> SimulationParams:
        udf = self._udf
        UDFVer = udf.getEngineVersion()

        algorithm = udf.get("Simulation_Conditions.Solver.Dynamics.Dynamics_Algorithm")
        if len(algorithm) == 0:
            raise RuntimeError("Error! Specify MD algorithm!")

        tail_correction = udf.get("Simulation_Conditions.Calc_Potential_Flags.Tail_Correction")
        pbc_a = udf.get("Simulation_Conditions.Boundary_Conditions.a_axis")
        pbc_b = udf.get("Simulation_Conditions.Boundary_Conditions.b_axis")
        pbc_c = udf.get("Simulation_Conditions.Boundary_Conditions.c_axis")

        fix_cell  = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Parrinello_Rahman_Nose_Hoover.Fix_Cell_Length")
        fix_angle = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Parrinello_Rahman_Nose_Hoover.Fix_Angle")

        deform, deform_npt, deform_vel = self._extract_deform(algorithm, UDFVer, cell)

        # --- vel_gen override from restart ---
        gen_method = udf.get("Initial_Structure.Generate_Method.Method")
        print(gen_method)
        if gen_method == "Restart":
            restore_vel = udf.get("Initial_Structure.Generate_Method.Restart.Restore_Velocity")
            print("restore_vel", restore_vel)
            vel_gen = (restore_vel == 0)

        # --- integrator ---
        integrator, ld_seed = self._determine_integrator(algorithm, deform_npt)

        # --- output intervals ---
        outputinterval = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps"
        )
        outputinterval2 = outputinterval
        if outputinterval >= 10000:
            outputinterval2 = int(outputinterval / 10)

        nsteps = udf.get("Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps")
        dt     = udf.get("Simulation_Conditions.Dynamics_Conditions.Time.delta_T", "[ps]")

        # --- constraints ---
        rattle_bond  = bool(udf.get("Simulation_Conditions.Dynamics_Conditions.RATTLE.Bond"))
        rattle_angle = bool(udf.get("Simulation_Conditions.Dynamics_Conditions.RATTLE.Angle"))

        # --- electrostatics ---
        qq_algorithm = ""
        if calcQQ == 1:
            qq_algorithm = udf.getArray("Interactions.Electrostatic_Interaction[].Algorithm", [0])

        # --- cutoffs ---
        cutoff_l = lj_cutoff
        cutoff_cl = cutoff_l
        if calcQQ == 1:
            cutoff_c = udf.get(
                "Interactions.Electrostatic_Interaction[].Cutoff_Coulomb.cutoff", [0], "[nm]"
            )
            if cutoff_c <= 0.0:
                cutoff_c = udf.get(
                    "Interactions.Electrostatic_Interaction[].Ewald.R_cutoff", [0], "[nm]"
                )
                if cutoff_c > 0.0:
                    print(" Ewald.R_cutoff was substituted for rcoulomb.")
            if cutoff_c <= cutoff_l:
                cutoff_c = cutoff_l
            cutoff_cl = max(cutoff_c, cutoff_l)
        else:
            cutoff_cl = cutoff_l

        # --- temperature coupling ---
        t_coupl_str, tau_t, ref_t = self._extract_temperature(algorithm, all_atm_num)

        # --- pressure coupling ---
        p_coupl_str, pcoupltype, tau_p, ref_p, ref_p_tensor, comp, comp_tensor = \
            self._extract_pressure(algorithm, fix_cell, fix_angle, deform_npt, deform_vel, cell)

        # --- pbc ---
        if pbc_a == "NONE" and pbc_b == "NONE" and pbc_c == "NONE":
            pbc = "no"
        elif pbc_a == "PERIODIC" and pbc_b == "PERIODIC" and pbc_c == "PERIODIC":
            pbc = "xyz"
        else:
            print(" Error: Please Check Your Boundary Conditions")
            pbc = "xyz"

        periodic_mol = bool(udf.get("Simulation_Conditions.Boundary_Conditions.Periodic_Bond"))

        # --- constraint freeze ---
        freeze_grps = None
        freeze_dim  = None
        loc_constr  = "Simulation_Conditions.Constraint_Conditions.Constraint_Atom[]"
        ndata_constr = udf.size(loc_constr)
        if ndata_constr > 0:
            constr_axis = udf.get(loc_constr + ".Constraint_Axis", [0])
            freeze_grps = "Constraint"
            parts = []
            for each_axis in constr_axis:
                parts.append("Y" if each_axis == "YES" else "N")
            freeze_dim = " ".join(parts)

        # --- vel_gen temp ---
        gen_temp = None
        if vel_gen:
            gen_temp = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Temperature.Temperature", "[K]"
            )

        return SimulationParams(
            title=title,
            algorithm=algorithm,
            nsteps=nsteps,
            dt=dt,
            outputinterval=outputinterval,
            outputinterval2=outputinterval2,
            integrator=integrator,
            ld_seed=ld_seed,
            vel_gen=vel_gen,
            gen_temp=gen_temp,
            rattle_bond=rattle_bond,
            rattle_angle=rattle_angle,
            calcQQ=calcQQ,
            qq_algorithm=qq_algorithm,
            lj_cutoff=cutoff_cl,
            coulomb_cutoff=cutoff_cl,
            t_coupl=t_coupl_str,
            tau_t=tau_t,
            ref_t=ref_t,
            p_coupl=p_coupl_str,
            pcoupltype=pcoupltype,
            tau_p=tau_p,
            ref_p=ref_p,
            ref_p_tensor=ref_p_tensor,
            compressibility=comp,
            compressibility_tensor=comp_tensor,
            tail_correction=tail_correction,
            pbc=pbc,
            periodic_mol=periodic_mol,
            deform_vel=deform_vel,
            freeze_grps=freeze_grps,
            freeze_dim=freeze_dim,
        )

    def _determine_integrator(self, algorithm: str, deform_npt: bool):
        ld_seed = None
        if "NPT" in algorithm:
            if deform_npt:
                integrator = "md"
            elif "NPT_Andersen" in algorithm:
                integrator = "md-vv"
            else:
                integrator = "md"
        elif "Kremer_Grest" in algorithm:
            integrator = "sd"
            ld_seed = 1993
        else:
            integrator = "md-vv"
        return integrator, ld_seed

    def _extract_deform(self, algorithm, UDFVer, cell):
        udf = self._udf
        deform = udf.get("Simulation_Conditions.Dynamics_Conditions.Deformation.Method")
        list_supported_deform = ["Cell_Deformation"]
        list_supported_deform_method = ["Simple_Elongation", "Deformation_Rate"]
        deform_npt = False
        deform_vel = None

        if len(deform) == 0:
            return deform, deform_npt, deform_vel

        cellsize = udf.get("Structure.Unit_Cell.Cell_Size")

        if deform in list_supported_deform:
            if "NPT" in algorithm:
                deform_npt = True
        else:
            raise RuntimeError(
                "Error!! deformation type '{}' is not supported.".format(deform)
            )

        deform_vel = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        deform_method = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Deformation.Cell_Deformation.Method"
        )
        if deform_method not in list_supported_deform_method:
            raise RuntimeError(
                "Error!! deformation method {} is not supported.".format(deform_method)
            )

        if deform_method == "Simple_Elongation":
            deform_axis = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Deformation"
                ".Cell_Deformation.Simple_Elongation.Axis"
            )
            poisson_ratio = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Deformation"
                ".Cell_Deformation.Simple_Elongation.Poisson_Ratio"
            )
            if UDFVer == "Ver101":
                input_method = udf.get(
                    "Simulation_Conditions.Dynamics_Conditions.Deformation"
                    ".Cell_Deformation.Simple_Elongation.Input_Method"
                )
                if input_method == "Deformation_Speed":
                    deform_rate = udf.get(
                        "Simulation_Conditions.Dynamics_Conditions.Deformation"
                        ".Cell_Deformation.Simple_Elongation.Deformation_Speed.Speed",
                        "[m/s]"
                    )
                elif input_method == "Initial_Strain_Rate":
                    if deform_axis == "z":
                        lz = cellsize[2] * 1e-10
                        rate = udf.get(
                            "Simulation_Conditions.Dynamics_Conditions.Deformation"
                            ".Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate",
                            "[1/s]"
                        )
                        deform_rate = lz * rate
            else:
                deform_rate = udf.get(
                    "Simulation_Conditions.Dynamics_Conditions.Deformation"
                    ".Cell_Deformation.Simple_Elongation.Elongation_Rate",
                    "[m/s]"
                )

            if poisson_ratio == 0.0:
                if deform_axis == "z":
                    deform_vel[2] = deform_rate * 0.001
                elif deform_axis == "xy":
                    deform_vel[0] = deform_rate * 0.001
                    deform_vel[1] = deform_rate * 0.001
                else:
                    raise RuntimeError(
                        'Error!! deform axis "{}" is unknown!'.format(deform_axis)
                    )
            else:
                raise RuntimeError("Error!! Poisson ratio other than 0.0 is not supported!")

        if deform_method == "Deformation_Rate":
            loc = "Simulation_Conditions.Dynamics_Conditions.Deformation.Cell_Deformation.Deformation_Rate"
            La = cellsize[0] * 1e-10
            Lb = cellsize[1] * 1e-10
            Lc = cellsize[2] * 1e-10
            deform_vel[0] = La * udf.get(loc + ".xx", "[Hz]") * 1e-3
            deform_vel[1] = Lb * udf.get(loc + ".yy", "[Hz]") * 1e-3
            deform_vel[2] = Lc * udf.get(loc + ".zz", "[Hz]") * 1e-3
            deform_vel[3] = Lb * udf.get(loc + ".xy", "[Hz]") * 1e-3
            deform_vel[4] = Lc * udf.get(loc + ".zx", "[Hz]") * 1e-3
            deform_vel[5] = Lc * udf.get(loc + ".yz", "[Hz]") * 1e-3

            if deform_npt:
                if deform_vel[3] != 0 or deform_vel[4] != 0 or deform_vel[5] != 0:
                    print("Error!! shear deformation with NPT ensemble is not supported!")

        return deform, deform_npt, deform_vel

    def _extract_temperature(self, algorithm, all_atm_num):
        udf = self._udf
        t_coupl_map = {
            "NVE":                              "no",
            "NVT_Nose_Hoover":                  "nose-hoover",
            "NPT_Parrinello_Rahman_Nose_Hoover": "nose-hoover",
            "NPT_Andersen_Nose_Hoover":          "nose-hoover",
            "NVT_Berendsen":                    "berendsen",
            "NPT_Berendsen":                    "berendsen",
            "NVT_Kremer_Grest":                 "no",
        }
        t_coupl_str = t_coupl_map.get(algorithm, "no")
        if t_coupl_str == "no" and algorithm not in t_coupl_map:
            print("algorithm %s is not supported !!" % algorithm)

        if t_coupl_str == "no":
            return t_coupl_str, 0.1, 300.0

        T_d = udf.get("Simulation_Conditions.Dynamics_Conditions.Temperature.Temperature", "[K]")
        unit_L    = udf.get("Unit_Parameter.Length", "[nm]")
        unit_Mass = udf.get("Unit_Parameter.Mass",   "[amu]")

        tau_t = 0.1
        if algorithm == "NVT_Nose_Hoover":
            Q   = udf.get("Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q")
            Q_d = Q * unit_Mass * unit_L * unit_L
            tau_t = math.sqrt((Q_d / T_d) * 4.0 * math.pi * math.pi)
        elif algorithm == "NPT_Parrinello_Rahman_Nose_Hoover":
            Q   = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Parrinello_Rahman_Nose_Hoover.Q")
            Q_d = Q * unit_Mass * unit_L * unit_L
            tau_t = math.sqrt((Q_d / T_d) * 4.0 * math.pi * math.pi)
        elif algorithm == "NPT_Andersen_Nose_Hoover":
            Q   = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Andersen_Nose_Hoover.Q")
            Q_d = Q * unit_Mass * unit_L * unit_L
            tau_t = math.sqrt((Q_d / T_d) * 4.0 * math.pi * math.pi)
        elif algorithm == "NVT_Berendsen":
            tau_t = udf.get("Simulation_Conditions.Solver.Dynamics.NVT_Berendsen.tau_T", "[ps]")
        elif algorithm == "NPT_Berendsen":
            tau_t = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_T", "[ps]")

        ref_t = T_d
        return t_coupl_str, tau_t, ref_t

    def _extract_pressure(self, algorithm, fix_cell, fix_angle, deform_npt, deform_vel, cell):
        udf = self._udf
        commp = 0.000045  # bar^-1

        p_coupl_map = {
            "NVE":                              "no",
            "NPT_Parrinello_Rahman_Nose_Hoover": None,  # handled below
            "NPT_Andersen_Nose_Hoover":          "MTTK",
            "NPT_Berendsen":                    "berendsen",
        }

        if algorithm == "NPT_Parrinello_Rahman_Nose_Hoover":
            if deform_npt:
                p_coupl_str = "Parrinello-Rahman"
            elif "NPT_Andersen" in algorithm:
                p_coupl_str = "MTTK"
            else:
                p_coupl_str = "Parrinello-Rahman"
        elif algorithm in p_coupl_map:
            p_coupl_str = p_coupl_map[algorithm] or "no"
        else:
            p_coupl_str = "no"

        if p_coupl_str == "no":
            return "no", "isotropic", 2.0, 1.0, None, commp, None

        # pcoupltype
        if algorithm == "NPT_Parrinello_Rahman_Nose_Hoover":
            if fix_cell in ("", "xy", "yz", "zx", "x", "y", "z") or deform_npt:
                pcoupltype = "anisotropic"
            else:
                print(" Error: Fix_Cell_Length {} is not supported!".format(fix_cell))
                pcoupltype = "isotropic"
        else:
            pcoupltype = "isotropic"

        unit_L    = udf.get("Unit_Parameter.Length", "[nm]")
        unit_Mass = udf.get("Unit_Parameter.Mass",   "[amu]")

        # cell dims for tau_p
        cell_x = cell.a
        cell_y = cell.b
        cell_z = cell.c
        max_L = max(cell_x, cell_y, cell_z)

        tau_p = 2.0
        if algorithm == "NPT_Parrinello_Rahman_Nose_Hoover":
            W   = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Parrinello_Rahman_Nose_Hoover.Cell_Mass")
            W_d = W * unit_Mass * unit_L * unit_L
            tau_p = math.sqrt((W_d * 4.0 * math.pi * math.pi * commp) / (3.0 * max_L))
            if tau_p < 2.0:
                tau_p = 2.0
        elif algorithm == "NPT_Andersen_Nose_Hoover":
            W   = udf.get("Simulation_Conditions.Solver.Dynamics.NPT_Andersen_Nose_Hoover.Cell_Mass")
            W   = W / 3.0  # Andersen -> Parrinello-Rahman
            W_d = W * unit_Mass * unit_L * unit_L
            tau_p = math.sqrt((W_d * 4.0 * math.pi * math.pi * commp) / (3.0 * max_L))
            if tau_p < 2.0:
                tau_p = 2.0
        elif algorithm == "NPT_Berendsen":
            unit_P = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure", "[bar]"
            ) / udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure", "[P]"
            )
            tau_p_raw = udf.get(
                "Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_P", "[P*ps]"
            )
            tau_p = tau_p_raw * unit_P * commp

        # ref_p
        pressure = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure", "[bar]"
        )
        strtensor = []
        for elem in ["xx", "yy", "zz", "yz", "zx", "xy"]:
            s = udf.get(
                "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Stress." + elem,
                "[bar]"
            )
            strtensor.append(s)

        ref_p_tensor = None
        if algorithm == "NPT_Parrinello_Rahman_Nose_Hoover":
            offdiag_idx_map = [0, 1, 2, 5, 4, 3]
            ref_p_tensor = []
            for i in offdiag_idx_map:
                if i < 3:
                    ref_p_tensor.append(pressure - strtensor[i])
                else:
                    ref_p_tensor.append(-strtensor[i])
            ref_p = pressure
        else:
            ref_p = pressure

        # compressibility
        comp_tensor = None
        if algorithm == "NPT_Parrinello_Rahman_Nose_Hoover":
            if fix_angle == 0:
                comp_offdiag = [commp, commp, commp]
            elif fix_angle == 1:
                comp_offdiag = [0.0, 0.0, 0.0]
            else:
                raise RuntimeError("Error: Please Check Your Fix_Cell_Angle")

            if fix_cell == "":
                comp_tensor = [commp, commp, commp] + comp_offdiag
            elif fix_cell == "xy":
                comp_tensor = [0.0, 0.0, commp] + comp_offdiag
            elif fix_cell == "yz":
                comp_tensor = [commp, 0.0, 0.0] + comp_offdiag
            elif fix_cell == "zx":
                comp_tensor = [0.0, commp, 0.0] + comp_offdiag
            elif fix_cell == "x":
                comp_tensor = [0.0, commp, commp] + comp_offdiag
            elif fix_cell == "y":
                comp_tensor = [commp, 0.0, commp] + comp_offdiag
            elif fix_cell == "z":
                comp_tensor = [commp, commp, 0.0] + comp_offdiag
            elif deform_npt and deform_vel is not None:
                if fix_angle == 1:
                    diag = []
                    for dv in deform_vel[:3]:
                        diag.append(0.0 if dv != 0.0 else commp)
                    comp_tensor = diag + [0.0, 0.0, 0.0]
                else:
                    raise RuntimeError("Error: Please Check Your Fix_Cell_Angle")
            else:
                raise RuntimeError("Error: Please Check Your Fix_Cell_Length")
            comp = commp
        else:
            comp = commp

        return p_coupl_str, pcoupltype, tau_p, ref_p, ref_p_tensor, comp, comp_tensor
