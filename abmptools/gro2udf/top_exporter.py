# -*- coding: utf-8 -*-
"""
top_exporter.py
---------------
Writes a :class:`TopModel` to a COGNAC UDF file using UDFManager.

This module ports the following functions from convert_gromacs_udf.py:
- add_set_of_molecules_byTop  → _write_set_of_molecules
- add_molecular_attributes_byTop → _write_molecular_attributes
- add_interactions            → _write_interactions
- append_structure            → _append_structure
- set_default_condition       → _set_default_condition
"""
from __future__ import annotations

import logging
import os
import shutil
from contextlib import contextmanager
from pathlib import Path
from typing import List, Optional

from .top_model import KB_AMU_A2_PS2_K, GROFrameData, TopModel
from .top_parser import TopParser
from .top_adapter import TopAdapter

logger = logging.getLogger(__name__)


class UDFExportError(RuntimeError):
    """Raised when writing a UDF field via UDFManager fails.

    The message carries the *section* of the export pipeline that was running
    (e.g. ``"Set_of_Molecules"``), the *template* / *output* paths, and the
    underlying UDFManager error so users can tell which schema field is
    incompatible with their UDFManager version (e.g. OCTA84 vs OCTA85).
    """


def _rewrite_cognac_include(udf_path: str, cognac_version: str) -> None:
    """Rewrite the ``\\include{"cognac<N>.udf"}`` directive in *udf_path*.

    Used so OCTA84 / J-OCTA 9.1 users (which only ship cognac110.udf or
    earlier) can still consume the bundled default_template.udf — which
    requests cognac112.udf — by passing ``--cognac-version 110``.

    The substitution targets the first occurrence of a ``\\include{"cognac
    <digits>.udf"}`` token; non-cognac includes are left untouched.
    """
    import re
    path = Path(udf_path)
    text = path.read_text()
    pattern = re.compile(r'\\include\{\s*"cognac\d+\.udf"\s*\}')
    if pattern.search(text) is None:
        # No cognac include to rewrite — silently leave the template alone
        # so non-cognac templates aren't disturbed.
        return
    new_directive = '\\include{"cognac' + str(cognac_version) + '.udf"}'
    # Pass replacement via lambda to bypass `re.sub`'s backslash interpretation.
    new_text = pattern.sub(lambda _m: new_directive, text, count=1)
    path.write_text(new_text)


@contextmanager
def _section(name: str, template_path: str, out_path: str):
    """Wrap an export section so UDFManager errors carry diagnostic context.

    UDFManager raises bare ``RuntimeError`` / ``IndexError`` with no hint about
    which field path was being written. That makes "OCTA84 schema is missing
    field X" failures very hard to diagnose. This context manager re-raises
    the original exception as :class:`UDFExportError` with the section name,
    template path, and output path attached.
    """
    try:
        yield
    except UDFExportError:
        raise
    except Exception as exc:
        hint = (
            "this often means the template UDF schema is incompatible with your\n"
            "    OCTA's UDFManager. The bundled template requests cognac11.2 (OCTA85).\n"
            "    OCTA8.4 / J-OCTA-9.1-Student only ships cognac10.1, and both the\n"
            "    `\\include` directive AND the data-section structure differ.\n"
            "    Recommended fix:\n"
            "      Open any minimal COGNAC UDF in your OCTA's GOURMET, hit\n"
            "      `File -> Save As` to dump it with your OCTA's native schema,\n"
            "      then pass that file via `--template <path>`. The data section\n"
            "      saved by GOURMET will match the cognac<N> schema available in\n"
            "      your install. See docs/gro2udf.md for full instructions."
        )
        raise UDFExportError(
            f"gro2udf: failed while writing section {name!r}.\n"
            f"  template UDF: {template_path}\n"
            f"  output  UDF: {out_path}\n"
            f"  underlying  : {type(exc).__name__}: {exc}\n"
            f"  hint        : {hint}"
        ) from exc


class TopExporter:
    """
    Orchestrate conversion from GROMACS TOP+GRO to COGNAC UDF.

    Usage::

        exporter = TopExporter()
        exporter.export("system.top", "output.gro",
                        template_path="template.udf",
                        out_path="result.udf")
    """

    def export(
        self,
        top_path: str,
        gro_path: str,
        template_path: str,
        out_path: str,
        mdp_path: Optional[str] = None,
        cognac_version: Optional[str] = None,
    ) -> None:
        """
        Parse *top_path* + *gro_path*, build :class:`TopModel`, write to *out_path*.

        Parameters
        ----------
        top_path      : path to GROMACS .top file (may include .itp via #include)
        gro_path      : path to GROMACS .gro file
        template_path : path to an existing COGNAC UDF file used as schema template
        out_path      : destination UDF path (created / overwritten)
        mdp_path      : optional path to GROMACS .mdp file; when provided,
                        Nose-Hoover Q and Ewald cutoff are computed from its values.
        cognac_version : optional override for the cognac<N>.udf include
                        directive in the template (e.g. ``"110"`` to fall back
                        from the bundled cognac112 default to the cognac110 schema
                        shipped with OCTA84 / J-OCTA 9.1).
        """
        raw = TopParser().parse(top_path)
        model = TopAdapter().build(raw, gro_path, mdp_path=mdp_path)
        self.export_model(model, template_path, out_path,
                          cognac_version=cognac_version)

    def export_model(
        self,
        model: TopModel,
        template_path: str,
        out_path: str,
        frames: Optional[List[GROFrameData]] = None,
        cognac_version: Optional[str] = None,
    ) -> None:
        """
        Write *model* into a new UDF at *out_path* using *template_path* as schema.

        Steps (mirrors convert_gromacs_udf.py ordering):
        1. copy template → out_path
        2. erase existing dynamic records
        3. write Set_of_Molecules (common record)
        4. append one Structure record per frame in *frames* (or ``model.frames``)
        5. set_default_condition  (electrostatic flag)
        6. write Molecular_Attributes (atom types, potentials)
        7. write Interactions (LJ pair potentials)

        Parameters
        ----------
        model : TopModel
            Topology + FF type definitions. ``model.frames`` is used when
            ``frames`` is not provided.
        template_path, out_path : str
            Template COGNAC UDF (schema source) and destination path.
        frames : list of GROFrameData, optional
            Override the structure records. Intended for callers that want
            the topology from a static build but coordinates from a longer
            trajectory (e.g. fcewsmb's ``setupgetcontact_gromacs`` injects
            xtc-sourced frames here via
            :func:`abmptools.amorphous.trajectory_ingest.frames_from_xtc`).
        """
        try:
            from UDFManager import UDFManager
        except ImportError as exc:
            raise UDFExportError(
                "UDFManager module is required but could not be imported. "
                "Install OCTA and set PYTHONPATH to its python3/ directory.\n"
                f"  underlying: {type(exc).__name__}: {exc}"
            ) from exc

        if not os.path.exists(template_path):
            raise UDFExportError(
                f"Template UDF not found: {template_path}.\n"
                f"  hint: use the bundled template at "
                f"`abmptools/gro2udf/default_template.udf`, or supply your own."
            )

        with _section("template-copy", template_path, out_path):
            shutil.copy(template_path, out_path)
            if cognac_version is not None:
                _rewrite_cognac_include(out_path, cognac_version)

        with _section("UDFManager-open", template_path, out_path):
            uobj = UDFManager(out_path)

        with _section("erase-existing-records", template_path, out_path):
            ntotal = uobj.totalRecord()
            if ntotal > 0:
                uobj.eraseRecord(0, ntotal)

        with _section("Set_of_Molecules", template_path, out_path):
            self._write_set_of_molecules(uobj, model)

        frames_to_write = frames if frames is not None else model.frames
        for i, frame in enumerate(frames_to_write):
            with _section(f"Structure[record={i}]", template_path, out_path):
                self._append_structure(uobj, model, frame)

        with _section("default_condition", template_path, out_path):
            self._set_default_condition(uobj, model)
        with _section("Molecular_Attributes", template_path, out_path):
            self._write_molecular_attributes(uobj, model)
        with _section("Interactions", template_path, out_path):
            self._write_interactions(uobj, model)

    # -------------------------------------------------------------------------
    # Writing helpers – all mirror convert_gromacs_udf.py functions
    # -------------------------------------------------------------------------

    @staticmethod
    def _write_set_of_molecules(uobj, model: TopModel) -> None:
        """Port of add_set_of_molecules_byTop."""
        uobj.jump(-1)

        e2Q = 18.224159264
        chargetype = "POINT_CHARGE"
        mol_type_names = model.mol_type_names
        mol_specs = model.mol_specs

        # --- atom / electrostatic / interaction site per molecule instance ---
        ncount = 0
        for imol, mol_name in enumerate(model.mol_instance_list):
            mol_idx = mol_type_names.index(mol_name)
            spec = mol_specs[mol_idx]

            uobj.put(mol_name, "Set_of_Molecules.molecule[].Mol_Name", [imol])

            for iatom, atom in enumerate(spec.atoms):
                uobj.put(ncount,
                         "Set_of_Molecules.molecule[].atom[].Atom_ID",
                         [imol, iatom])
                # Atom_Name と Atom_Type_Name は両方 element symbol を書く:
                # - Atom_Name: ABINIT-MP の read_pdb / obabel の xyz→pdb 変換
                #   で元素として解釈される
                # - Atom_Type_Name: J-OCTA viewer の atom テーブル `type_name`
                #   列に表示され、ここに per-atom unique な ``MOL0_X`` (OpenFF
                #   SMIRNOFF 由来) が入っていると element として認識されず
                #   色付け / フィルタが機能しない (2026-05-26 ユーザー報告)。
                #   よって element symbol を書く。
                # GROMACS 由来の atom 名 (e.g. "c30", "hc1") は GAFF 型名形式
                # で、obabel が atom type を `*` に変換し ABINIT-MP が "Atom
                # type * is not supported" で失敗するため使えない。
                # GAFF / SMIRNOFF の元の per-atom unique type は引き続き
                # ``Set_of_Molecules.molecule[].interaction_Site[].Type_Name``
                # と ``Molecular_Attributes.Atom_Type[]`` 側に保持されるので、
                # LJ Pair_Interaction (per-atom-type の σ/ε) や MD 実行に
                # 必要な per-atom unique typing は壊れない。
                uobj.put(atom.element,
                         "Set_of_Molecules.molecule[].atom[].Atom_Name",
                         [imol, iatom])
                uobj.put(atom.element,
                         "Set_of_Molecules.molecule[].atom[].Atom_Type_Name",
                         [imol, iatom])
                uobj.put(0,
                         "Set_of_Molecules.molecule[].atom[].Chirality",
                         [imol, iatom])
                uobj.put(1,
                         "Set_of_Molecules.molecule[].atom[].Main_Chain",
                         [imol, iatom])

                uobj.put(chargetype,
                         "Set_of_Molecules.molecule[].electrostatic_Site[].Type_Name",
                         [imol, iatom])
                uobj.put(atom.charge * e2Q,
                         "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element",
                         [imol, iatom])
                uobj.put(iatom,
                         "Set_of_Molecules.molecule[].electrostatic_Site[].atom[]",
                         [imol, iatom, 0])

                uobj.put(atom.type_name,
                         "Set_of_Molecules.molecule[].interaction_Site[].Type_Name",
                         [imol, iatom])
                uobj.put(iatom,
                         "Set_of_Molecules.molecule[].interaction_Site[].atom[]",
                         [imol, iatom, 0])
                ncount += 1

        # --- bonds / angles / torsions per molecule instance ---
        for imol, mol_name in enumerate(model.mol_instance_list):
            mol_idx = mol_type_names.index(mol_name)
            spec = mol_specs[mol_idx]

            for idat, bond in enumerate(spec.bonds):
                uobj.put(bond.potential_name,
                         "Set_of_Molecules.molecule[].bond[].Potential_Name",
                         [imol, idat])
                uobj.put(bond.atom1 - 1,
                         "Set_of_Molecules.molecule[].bond[].atom1",
                         [imol, idat])
                uobj.put(bond.atom2 - 1,
                         "Set_of_Molecules.molecule[].bond[].atom2",
                         [imol, idat])
                uobj.put(1.0,
                         "Set_of_Molecules.molecule[].bond[].Order",
                         [imol, idat])

            for idat, angle in enumerate(spec.angles):
                uobj.put(angle.potential_name,
                         "Set_of_Molecules.molecule[].angle[].Potential_Name",
                         [imol, idat])
                uobj.put(angle.atom1 - 1,
                         "Set_of_Molecules.molecule[].angle[].atom1",
                         [imol, idat])
                uobj.put(angle.atom2 - 1,
                         "Set_of_Molecules.molecule[].angle[].atom2",
                         [imol, idat])
                uobj.put(angle.atom3 - 1,
                         "Set_of_Molecules.molecule[].angle[].atom3",
                         [imol, idat])

            for idat, tors in enumerate(spec.torsions):
                uobj.put(tors.potential_name,
                         "Set_of_Molecules.molecule[].torsion[].Potential_Name",
                         [imol, idat])
                uobj.put(tors.atom1 - 1,
                         "Set_of_Molecules.molecule[].torsion[].atom1",
                         [imol, idat])
                uobj.put(tors.atom2 - 1,
                         "Set_of_Molecules.molecule[].torsion[].atom2",
                         [imol, idat])
                uobj.put(tors.atom3 - 1,
                         "Set_of_Molecules.molecule[].torsion[].atom3",
                         [imol, idat])
                uobj.put(tors.atom4 - 1,
                         "Set_of_Molecules.molecule[].torsion[].atom4",
                         [imol, idat])

        uobj.write()

    @staticmethod
    def _append_structure(uobj, model: TopModel, frame: GROFrameData) -> None:
        """Port of append_structure (one GRO frame → one UDF dynamic record)."""
        uobj.newRecord()

        mol_type_names = model.mol_type_names
        mol_specs = model.mol_specs
        ncount = 0

        for imol, mol_name in enumerate(model.mol_instance_list):
            mol_idx = mol_type_names.index(mol_name)
            spec = mol_specs[mol_idx]
            for iatom in range(len(spec.atoms)):
                posx, posy, posz = frame.coord_list[ncount]
                uobj.put(posx, "Structure.Position.mol[].atom[].x",
                         [imol, iatom], "[nm]")
                uobj.put(posy, "Structure.Position.mol[].atom[].y",
                         [imol, iatom], "[nm]")
                uobj.put(posz, "Structure.Position.mol[].atom[].z",
                         [imol, iatom], "[nm]")
                ncount += 1

        cell = frame.cell
        uobj.put(cell[0], "Structure.Unit_Cell.Cell_Size.a", "[nm]")
        uobj.put(cell[1], "Structure.Unit_Cell.Cell_Size.b", "[nm]")
        uobj.put(cell[2], "Structure.Unit_Cell.Cell_Size.c", "[nm]")
        uobj.put(90.0, "Structure.Unit_Cell.Cell_Size.alpha")
        uobj.put(90.0, "Structure.Unit_Cell.Cell_Size.beta")
        uobj.put(90.0, "Structure.Unit_Cell.Cell_Size.gamma")
        uobj.put(frame.step, "Steps")
        uobj.put(frame.time, "Time", "[ps]")
        uobj.put(0.0, "Statistics_Data.Energy.Instantaneous.Hamiltonian")
        uobj.write()

    @staticmethod
    def _set_default_condition(uobj, model: TopModel) -> None:
        """
        Set electrostatic flags, Nose-Hoover Q, and Ewald parameters.

        Nose-Hoover Q
        -------------
        Q [amu·Å²] = g · k_B · T · τ²

        where:
          g   = 3·N_atoms − 3   (degrees of freedom; COM motion removed)
          k_B = 0.83144626      amu·Å² / (ps²·K)  (Boltzmann constant in
                                COGNAC internal units: amu, Å, ps)
          T   = model.ref_t     [K]   (from MDP ref_t, default 300.0)
          τ   = model.tau_t     [ps]  (from MDP tau_t, default 0.1)

        Ewald defaults
        --------------
        Name                = "POINT_CHARGE"
        Algorithm           = "Ewald"
        Scale_1_4_Pair      = 0.83333333333333   (5/6, AMBER convention)
        Ewald.Dielectric_Constant = 0.0
        Ewald.R_cutoff      = model.ewald_r_cutoff [Å]  (Deserno & Holm
                              formula from GRO box; fallback 10.0 Å)
        Ewald.Ewald_Parameters = "Auto"
        """
        uobj.jump(-1)

        # --- electrostatic flag ---
        uobj.put(1, "Simulation_Conditions.Calc_Potential_Flags.Electrostatic")

        # --- Nose-Hoover Q ---
        n = model.n_atoms_total
        if n > 0:
            g = max(1, 3 * n - 3)           # degrees of freedom
            Q = g * KB_AMU_A2_PS2_K * model.ref_t * (model.tau_t ** 2)
            uobj.put(Q,
                     "Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q")

        # --- Ewald electrostatic interaction defaults ---
        loc = "Interactions.Electrostatic_Interaction[0]"
        uobj.put("POINT_CHARGE",      loc + ".Name")
        uobj.put("Ewald",             loc + ".Algorithm")
        uobj.put(0.83333333333333,    loc + ".Scale_1_4_Pair")
        uobj.put(0.0,                 loc + ".Ewald.Dielectric_Constant")
        uobj.put(model.ewald_r_cutoff, loc + ".Ewald.R_cutoff")
        uobj.put("Auto",              loc + ".Ewald.Ewald_Parameters")

        uobj.write()

    @staticmethod
    def _write_molecular_attributes(uobj, model: TopModel) -> None:
        """Port of add_molecular_attributes_byTop."""
        uobj.jump(-1)

        # --- Atom_Type (unique types in order of first appearance) ---
        seen_types: List[str] = []
        for spec in model.mol_specs:
            for atom in spec.atoms:
                if atom.type_name not in seen_types:
                    seen_types.append(atom.type_name)

        for ndata, atype in enumerate(seen_types):
            uobj.put(atype,
                     "Molecular_Attributes.Atom_Type[].Name",
                     [ndata])
            mass = model.mass_dict.get(atype, 0.0)
            uobj.put(mass,
                     "Molecular_Attributes.Atom_Type[].Mass",
                     [ndata])

        # --- Bond_Potential (Harmonic) ---
        for j, bt in enumerate(model.bond_type_specs):
            uobj.put(bt.name,
                     "Molecular_Attributes.Bond_Potential[].Name", [j])
            uobj.put("Harmonic",
                     "Molecular_Attributes.Bond_Potential[].Potential_Type", [j])
            uobj.put(bt.r0,
                     "Molecular_Attributes.Bond_Potential[].R0", [j], "[nm]")
            uobj.put(bt.kb,
                     "Molecular_Attributes.Bond_Potential[].Harmonic.K",
                     [j], "[kJ/mol/nm^2]")

        # --- Angle_Potential (Theta) ---
        for j, at in enumerate(model.angle_type_specs):
            q0 = 180.0 - at.theta0   # COGNAC convention: supplement of GROMACS theta0
            uobj.put(at.name,
                     "Molecular_Attributes.Angle_Potential[].Name", [j])
            uobj.put("Theta",
                     "Molecular_Attributes.Angle_Potential[].Potential_Type", [j])
            uobj.put(q0,
                     "Molecular_Attributes.Angle_Potential[].theta0", [j])
            uobj.put(at.k,
                     "Molecular_Attributes.Angle_Potential[].Theta.K",
                     [j], "[kJ/mol/rad^2]")

        # --- Torsion_Potential ---
        # Each TorsionTypeSpec is already one term (params = [PHASE, PK, PN] for Amber).
        for j, tt in enumerate(model.torsion_type_specs):
            funct = tt.funct
            params = tt.params

            if funct in (1, 9, 4):
                # Amber: params = [PHASE, PK, PN]
                PHASE = params[0]
                PK = params[1]
                PN = params[2]
                uobj.put(tt.name,
                         "Molecular_Attributes.Torsion_Potential[].Name", [j])
                uobj.put("Amber",
                         "Molecular_Attributes.Torsion_Potential[].Potential_Type", [j])
                uobj.put(PK,
                         "Molecular_Attributes.Torsion_Potential[].Amber.PK",
                         [j], "[kJ/mol]")
                uobj.put(1,
                         "Molecular_Attributes.Torsion_Potential[].Amber.IDIVF", [j])
                uobj.put(PN,
                         "Molecular_Attributes.Torsion_Potential[].Amber.PN", [j])
                uobj.put(PHASE,
                         "Molecular_Attributes.Torsion_Potential[].Amber.PHASE", [j])
                uobj.put(1,
                         "Molecular_Attributes.Torsion_Potential[].Amber.trans_is_0",
                         [j])

            elif funct == 3:
                # Ryckaert-Bellemans / Cosine_Polynomial
                nparams = len(params)
                while nparams > 0 and params[nparams - 1] == 0.0:
                    nparams -= 1
                uobj.put(tt.name,
                         "Molecular_Attributes.Torsion_Potential[].Name", [j])
                uobj.put("Cosine_Polynomial",
                         "Molecular_Attributes.Torsion_Potential[].Potential_Type", [j])
                uobj.put(nparams,
                         "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.N",
                         [j])
                uobj.put(1.0,
                         "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.K",
                         [j], "[kJ/mol]")
                for k in range(nparams):
                    uobj.put(params[k],
                             "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.p[]",
                             [j, k])
            else:
                logger.warning("torsion funct=%s is not supported (%s)",
                               funct, tt.name)

        uobj.write()

    @staticmethod
    def _write_interactions(uobj, model: TopModel) -> None:
        """Port of add_interactions."""
        uobj.jump(-1)

        # Collect atom types actually referenced in Set_of_Molecules
        used_types: set = set()
        for spec in model.mol_specs:
            for atom in spec.atoms:
                used_types.add(atom.type_name)

        # Molecular_Attributes.Interaction_Site_Type
        ncount = 0
        for at in model.atom_type_specs:
            if at.name not in used_types:
                continue
            rval = at.sigma * 1.5
            uobj.put(at.name,
                     "Molecular_Attributes.Interaction_Site_Type[].Name",
                     [ncount])
            uobj.put(1,
                     "Molecular_Attributes.Interaction_Site_Type[].Num_of_Atoms",
                     [ncount])
            uobj.put(rval,
                     "Molecular_Attributes.Interaction_Site_Type[].Range",
                     [ncount], "[nm]")
            uobj.put(0,
                     "Molecular_Attributes.Interaction_Site_Type[].Rigid",
                     [ncount])
            ncount += 1

        # Interactions.Pair_Interaction (Lennard-Jones)
        s2c = 2.5
        combrule = model.comb_rule
        scale = model.fudge_lj
        location = "Interactions.Pair_Interaction[]"
        types_used = [at for at in model.atom_type_specs if at.name in used_types]
        ntypes = len(types_used)
        ncount = 0

        for j1 in range(ntypes):
            for j2 in range(j1, ntypes):
                at1 = types_used[j1]
                at2 = types_used[j2]
                s1, s2 = at1.sigma, at2.sigma
                e1, e2 = at1.epsilon, at2.epsilon

                if combrule == 2:       # Lorentz-Berthelot
                    sigma = (s1 + s2) * 0.5
                elif combrule == 3:     # OPLS geometric
                    sigma = (s1 * s2) ** 0.5
                else:
                    sigma = (s1 + s2) * 0.5

                if sigma < 1.0e-10:
                    continue

                epsilon = (e1 * e2) ** 0.5
                cutoff = sigma * s2c
                name = "{}-{}".format(at1.name, at2.name)

                uobj.put(name,          location + ".Name",           [ncount])
                uobj.put("Lennard_Jones", location + ".Potential_Type", [ncount])
                uobj.put(at1.name,      location + ".Site1_Name",     [ncount])
                uobj.put(at2.name,      location + ".Site2_Name",     [ncount])
                uobj.put(cutoff,        location + ".Cutoff",         [ncount], "[nm]")
                uobj.put(scale,         location + ".Scale_1_4_Pair", [ncount])
                uobj.put(sigma,         location + ".Lennard_Jones.sigma",
                         [ncount], "[nm]")
                uobj.put(epsilon,       location + ".Lennard_Jones.epsilon",
                         [ncount], "[kJ/mol]")
                ncount += 1

        uobj.write()
