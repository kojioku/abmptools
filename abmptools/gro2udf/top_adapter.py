# -*- coding: utf-8 -*-
"""
top_adapter.py
--------------
Converts raw parsed data (:class:`TopRawData`) into the typed intermediate
representation (:class:`TopModel`).

Responsibilities:
- Apply mass-based element-name mapping (mass_to_element).
- Build BondTypeSpec / AngleTypeSpec / TorsionTypeSpec from raw lists.
- Assign global_atom_id to each atom across all molecule instances.
- Read GRO frames and map them to the flat coord_list / cell structure
  that TopExporter expects.
"""
from __future__ import annotations

from typing import List, Optional

from .gro_parser import GROParser
from .top_model import (
    AtomTypeSpec,
    BondTypeSpec,
    AngleTypeSpec,
    GROFrameData,
    MolAngleSpec,
    MolAtomSpec,
    MolBondSpec,
    MolSpec,
    MolTorsionSpec,
    TopModel,
    TorsionTypeSpec,
    mass_to_element,
)
from .top_parser import TopRawData, get_bond_name, get_angle_name, get_torsion_name


class TopAdapter:
    """Build a :class:`TopModel` from a :class:`TopRawData` and a GRO path."""

    def build(
        self,
        raw: TopRawData,
        gro_path: str,
    ) -> TopModel:
        """
        Convert *raw* (from TopParser) + GRO frames into a TopModel.

        Parameters
        ----------
        raw       : TopRawData returned by TopParser.parse()
        gro_path  : path to the GROMACS .gro file
        """
        # --- mass dictionary: type_name -> mass ---
        mass_dict = self._build_mass_dict(raw)

        # --- atom type specs ---
        atom_type_specs = [
            AtomTypeSpec(
                name=at[0],
                mass=at[1],
                sigma=at[3],
                epsilon=at[4],
            )
            for at in raw.atomtypes
        ]

        # --- bond / angle / torsion type specs (from inline mol params) ---
        bond_type_specs = self._build_bond_type_specs(raw.bond_types_from_mol)
        angle_type_specs = self._build_angle_type_specs(raw.angle_types_from_mol)
        torsion_type_specs = self._build_torsion_type_specs(
            raw.torsion_types_from_mol)

        # --- molecule specs (one per unique type) ---
        mol_type_names = list(raw.mol_types)
        mol_specs = self._build_mol_specs(
            raw, mass_dict, bond_type_specs, angle_type_specs,
            torsion_type_specs)

        # --- assign global_atom_id across all instances ---
        self._assign_global_ids(mol_specs, mol_type_names,
                                raw.mol_instance_list)

        # --- GRO frames ---
        frames = self._read_gro_frames(gro_path)

        return TopModel(
            comb_rule=raw.comb_rule,
            fudge_lj=raw.fudge_lj,
            fudge_qq=raw.fudge_qq,
            atom_type_specs=atom_type_specs,
            bond_type_specs=bond_type_specs,
            angle_type_specs=angle_type_specs,
            torsion_type_specs=torsion_type_specs,
            mass_dict=mass_dict,
            mol_type_names=mol_type_names,
            mol_specs=mol_specs,
            mol_instance_list=list(raw.mol_instance_list),
            frames=frames,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _build_mass_dict(raw: TopRawData):
        """Build type_name → mass dict, filling gaps from atomlist."""
        mass_dict = {}
        for at in raw.atomtypes:
            mass_dict[at[0]] = at[1]
        # Fill missing entries from per-molecule atom data
        for atomlist_mol in raw.atomlist:
            for atom in atomlist_mol:
                atype = atom[1]
                if atype not in mass_dict or mass_dict[atype] == 0.0:
                    m = atom[3]
                    if m > 0.0:
                        mass_dict[atype] = m
        return mass_dict

    @staticmethod
    def _build_bond_type_specs(bond_types_from_mol: list) -> List[BondTypeSpec]:
        specs = []
        for j, bt in enumerate(bond_types_from_mol):
            a1, a2, funct, params = bt
            name = get_bond_name(a1, a2) + "-" + str(j)
            specs.append(BondTypeSpec(
                name=name,
                name1=a1,
                name2=a2,
                funct=funct,
                r0=params[0],
                kb=params[1],
            ))
        return specs

    @staticmethod
    def _build_angle_type_specs(angle_types_from_mol: list) -> List[AngleTypeSpec]:
        specs = []
        for j, at in enumerate(angle_types_from_mol):
            a1, a2, a3, funct, params = at
            name = get_angle_name(a1, a2, a3) + "-" + str(j)
            specs.append(AngleTypeSpec(
                name=name,
                name1=a1,
                name2=a2,
                name3=a3,
                funct=funct,
                theta0=params[0],
                k=params[1],
            ))
        return specs

    @staticmethod
    def _build_torsion_type_specs(torsion_types_from_mol: list) -> List[TorsionTypeSpec]:
        specs = []
        for jj, tt in enumerate(torsion_types_from_mol):
            a1, a2, a3, a4, funct, improper, params = tt
            basename = get_torsion_name(a1, a2, a3, a4) + "-" + str(jj)
            if funct in (2, 4):
                improper = True
                basename += "-oopa"
            elif improper:
                basename += "-oopa"

            if funct in (1, 9, 4):
                # Amber: params = [phase, k, mult, phase2, k2, mult2, ...]
                # Store one triplet per TorsionTypeSpec entry.
                nmulti = len(params) // 3
                for k in range(nmulti):
                    name = basename if nmulti == 1 else basename + ":" + str(k)
                    specs.append(TorsionTypeSpec(
                        name=name,
                        name1=a1, name2=a2, name3=a3, name4=a4,
                        funct=funct,
                        improper=improper,
                        params=params[k * 3:(k + 1) * 3],  # one [phase, PK, PN] slice
                    ))
            else:
                specs.append(TorsionTypeSpec(
                    name=basename,
                    name1=a1, name2=a2, name3=a3, name4=a4,
                    funct=funct,
                    improper=improper,
                    params=params,
                ))
        return specs

    def _build_mol_specs(
        self,
        raw: TopRawData,
        mass_dict: dict,
        bond_type_specs: List[BondTypeSpec],
        angle_type_specs: List[AngleTypeSpec],
        torsion_type_specs: List[TorsionTypeSpec],
    ) -> List[MolSpec]:
        specs = []
        for i, mol_name in enumerate(raw.mol_types):
            atomlist_mol = raw.atomlist[i]
            bondlist_mol = raw.bondlist[i]
            anglelist_mol = raw.anglelist[i]
            torsionlist_mol = raw.torsionlist[i]

            # atoms
            atoms: List[MolAtomSpec] = []
            for idx, atom in enumerate(atomlist_mol):
                atom_name, atype, charge, atom_mass = atom
                mass = mass_dict.get(atype, atom_mass)
                element = mass_to_element(mass)
                atoms.append(MolAtomSpec(
                    index_1based=idx + 1,
                    atom_name=atom_name,
                    element=element,
                    type_name=atype,
                    charge=charge,
                    global_atom_id=0,  # filled in _assign_global_ids
                ))

            # bonds
            bonds: List[MolBondSpec] = []
            for bd in bondlist_mol:
                a1, a2, tidx = bd
                pot_name = bond_type_specs[tidx].name if tidx < len(bond_type_specs) else "bond-{}".format(tidx)
                bonds.append(MolBondSpec(
                    atom1=a1,
                    atom2=a2,
                    type_index=tidx,
                    potential_name=pot_name,
                ))

            # angles
            angles: List[MolAngleSpec] = []
            for ag in anglelist_mol:
                a1, a2, a3, tidx = ag
                pot_name = angle_type_specs[tidx].name if tidx < len(angle_type_specs) else "angle-{}".format(tidx)
                angles.append(MolAngleSpec(
                    atom1=a1,
                    atom2=a2,
                    atom3=a3,
                    type_index=tidx,
                    potential_name=pot_name,
                ))

            # torsions → expand multi-term Amber
            torsions: List[MolTorsionSpec] = []
            j_tors = 0  # index into torsion_type_specs
            for td in torsionlist_mol:
                if len(td) < 6:
                    # type-reference only (no inline params)
                    a1, a2, a3, a4, ptype = td
                    pot_name = "torsion-{}".format(ptype)
                    torsions.append(MolTorsionSpec(
                        atom1=a1, atom2=a2, atom3=a3, atom4=a4,
                        funct=1, improper=False, n_params=3,
                        potential_name=pot_name,
                    ))
                    continue

                a1, a2, a3, a4, tidx, tfunc, timp, nparams = td
                base_name = (torsion_type_specs[tidx].name
                             if tidx < len(torsion_type_specs)
                             else "torsion-{}".format(tidx))

                if tfunc in (1, 9, 4):
                    nmulti = nparams // 3
                    for k in range(nmulti):
                        if nmulti == 1:
                            pot_name = base_name
                        else:
                            # strip any existing ":n" suffix and re-add
                            if ":" in base_name:
                                base_name_clean = base_name.rsplit(":", 1)[0]
                            else:
                                base_name_clean = base_name
                            pot_name = base_name_clean + ":" + str(k)
                        torsions.append(MolTorsionSpec(
                            atom1=a1, atom2=a2, atom3=a3, atom4=a4,
                            funct=tfunc, improper=timp, n_params=nparams,
                            potential_name=pot_name,
                        ))
                        j_tors += 1
                else:
                    torsions.append(MolTorsionSpec(
                        atom1=a1, atom2=a2, atom3=a3, atom4=a4,
                        funct=tfunc, improper=timp, n_params=nparams,
                        potential_name=base_name,
                    ))
                    j_tors += 1

            specs.append(MolSpec(
                name=mol_name,
                atoms=atoms,
                bonds=bonds,
                angles=angles,
                torsions=torsions,
            ))
        return specs

    @staticmethod
    def _assign_global_ids(
        mol_specs: List[MolSpec],
        mol_type_names: List[str],
        mol_instance_list: List[str],
    ) -> None:
        """Assign sequential 0-based global_atom_id to each atom entry."""
        gid = 0
        for mol_name in mol_instance_list:
            idx = mol_type_names.index(mol_name)
            for atom in mol_specs[idx].atoms:
                atom.global_atom_id = gid
                gid += 1

    @staticmethod
    def _read_gro_frames(gro_path: str) -> List[GROFrameData]:
        """Read all frames from *gro_path* using GROParser."""
        parser = GROParser()
        frames: List[GROFrameData] = []
        for frame in parser.parse_frames(gro_path):
            coord_list = [[a.x, a.y, a.z] for a in frame.atoms]
            cell = frame.box_vals[:3]  # use only first 3 (x,y,z lengths)
            frames.append(GROFrameData(
                step=frame.step,
                time=frame.time,
                coord_list=coord_list,
                cell=cell,
            ))
        return frames
