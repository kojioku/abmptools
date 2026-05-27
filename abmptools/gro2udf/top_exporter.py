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


_OPENFF_TYPE_RE = __import__("re").compile(r"^MOL\d+_(\d+)$")


def _load_trajectory_frames(trajectory_path: str, gro_path: str):
    """Load multi-frame coordinates from a .gro (multi-frame) or .xtc.

    The single-frame .gro at *gro_path* is only used to provide the atom
    count expected by the topology (for sanity-checking the trajectory).
    """
    ext = os.path.splitext(trajectory_path)[1].lower()
    if ext == ".gro":
        from .trajectory_io import frames_from_multi_gro
        return frames_from_multi_gro(trajectory_path)
    if ext == ".xtc":
        # Reuse the existing MDAnalysis-backed loader from abmptools.amorphous.
        from ..amorphous.trajectory_ingest import frames_from_xtc
        # frames_from_xtc requires the topology + xtc + the GROMACS .top.
        # The .top path is reachable via the caller's TopExporter.export
        # invocation; here we pass gro_path (= topology .gro) only.
        return frames_from_xtc(topology_path=gro_path,
                               xtc_path=trajectory_path,
                               top_path=None)
    raise ValueError(
        f"Unsupported --trajectory extension {ext!r}; expected .gro or .xtc"
    )


def _load_energy_series(xvg_path: str):
    """Read a .xvg energy file into ``(times, {legend: values})``."""
    from .trajectory_io import read_xvg
    return read_xvg(xvg_path)


# Mapping from xvg legend names (gmx energy output) to UDF Statistics_Data
# paths. Each value is a (relative-path, unit) tuple; multiple legends can
# fold into the same path (e.g. Proper + Improper -> Energy.Instantaneous
# .Torsion), in which case :func:`_aggregate_statistics_per_frame` sums them.
# Unit is passed to UDFManager via _put_with_unit_fallback so it's converted
# to the schema's native unit (epsilon / P / mass/sigma^3 / T / sigma^3).
_XVG_TO_UDF_STATS = {
    # Energy components (Energy class, native unit [epsilon] = kJ/mol)
    "Bond":          ("Energy.Instantaneous.Bond",          "[kJ/mol]"),
    "Angle":         ("Energy.Instantaneous.Angle",         "[kJ/mol]"),
    "Proper Dih.":   ("Energy.Instantaneous.Torsion",       "[kJ/mol]"),
    "Improper Dih.": ("Energy.Instantaneous.Torsion",       "[kJ/mol]"),
    "LJ-14":         ("Energy.Instantaneous.Nonbonding",    "[kJ/mol]"),
    "LJ (SR)":       ("Energy.Instantaneous.Nonbonding",    "[kJ/mol]"),
    "Disper. corr.": ("Energy.Instantaneous.Nonbonding",    "[kJ/mol]"),
    "Coulomb-14":    ("Energy.Instantaneous.Electrostatic", "[kJ/mol]"),
    "Coulomb (SR)":  ("Energy.Instantaneous.Electrostatic", "[kJ/mol]"),
    "Coul. recip.":  ("Energy.Instantaneous.Electrostatic", "[kJ/mol]"),
    "Potential":     ("Energy.Instantaneous.Potential",     "[kJ/mol]"),
    "Kinetic En.":   ("Energy.Instantaneous.Kinetic",       "[kJ/mol]"),
    "Total Energy":  ("Energy.Instantaneous.Total",         "[kJ/mol]"),
    # Bulk thermodynamic observables (separate top-level Statistics_Data
    # classes with their own native units)
    "Temperature":   ("Temperature.Instantaneous",          "[K]"),
    "Pressure":      ("Pressure.Instantaneous",             "[bar]"),
    "Density":       ("Density.Instantaneous",              "[kg/m^3]"),
    "Volume":        ("Volume.Instantaneous",               "[nm^3]"),
}


def _aggregate_statistics_per_frame(frames, energy_times, energy_series):
    """For each frame, build ``{(relative_udf_path, unit): float}`` from xvg.

    Maps xvg legends through :data:`_XVG_TO_UDF_STATS` to Statistics_Data
    paths (relative to ``Statistics_Data.``). The xvg time grid (typically
    denser than the trajectory grid) is matched to each frame's
    ``frame.time`` via nearest-neighbour lookup. Legends mapped to the same
    UDF path are summed (e.g. Proper Dih. + Improper Dih. -> Torsion;
    LJ-14 + LJ(SR) + Disper.corr. -> Nonbonding).

    Returns ``None`` when no energy data is available.
    """
    if not energy_times or not energy_series:
        return None

    # Pre-compute the xvg row index for each frame time.
    row_for_frame = []
    n_xvg = len(energy_times)
    for frame in frames:
        if frame is None:
            row_for_frame.append(None)
            continue
        target = frame.time
        best = 0
        best_diff = abs(energy_times[0] - target)
        for i in range(1, n_xvg):
            diff = abs(energy_times[i] - target)
            if diff < best_diff:
                best = i
                best_diff = diff
        row_for_frame.append(best)

    per_frame = []
    for row_idx in row_for_frame:
        if row_idx is None:
            per_frame.append({})
            continue
        agg: dict = {}  # (path, unit) -> float
        for legend_name, (udf_path, unit) in _XVG_TO_UDF_STATS.items():
            values = energy_series.get(legend_name)
            if not values or row_idx >= len(values):
                continue
            key = (udf_path, unit)
            agg[key] = agg.get(key, 0.0) + values[row_idx]
        per_frame.append(agg)
    return per_frame


def _display_type_name(type_name: str, element: str) -> str:
    """Convert an OpenFF interchange-style ``MOL0_<N>`` atom-type name to
    an ``<element><N>`` form (e.g. ``MOL0_4`` -> ``C4``), so J-OCTA's
    3D viewer can infer the element from the leading character.

    Non-OpenFF type names (GAFF ``c3``, OPLS ``opls_267``, etc.) are passed
    through unchanged — they already start with the element letter by
    convention.

    Why: OpenFF SMIRNOFF / interchange assigns per-atom unique type names
    starting with the prefix ``MOL0_``. J-OCTA treats the leading character
    as the element symbol; "M" is interpreted as Mg / Mn and the viewer
    falls back to assigning a different random color to every atom type,
    breaking element-based CPK rendering. Prefixing with the actual element
    letter restores correct rendering while preserving the per-atom uniqueness
    needed by LJ Pair_Interaction.
    """
    m = _OPENFF_TYPE_RE.match(type_name)
    if m:
        return f"{element}{m.group(1)}"
    return type_name


def _put_with_unit_fallback(uobj, value, path, indices=None, unit=None):
    """Wrap ``uobj.put`` so callers can request a unit alias (e.g. ``"[nm]"``)
    that may not exist in older cognac schemas.

    cognac11.2 (OCTA85) declares unit aliases ``[nm]`` / ``[ps]``;
    cognac10.1 (OCTA8.4 / J-OCTA-9.1-Student) does not. Passing ``[nm]`` to
    UDFManager on cognac10.1 raises ``RuntimeError: ArgumentError: put data.``.
    We therefore attempt the put with the requested unit first, and if it
    fails fall back to a unit-less put so the value lands as-is in the
    schema's native unit (typically ``[sigma]`` for cognac), keeping older
    OCTA versions working without a separate code path.
    """
    if unit is not None:
        try:
            if indices is None:
                uobj.put(value, path, unit)
            else:
                uobj.put(value, path, indices, unit)
            return
        except Exception:
            pass  # fall through to the no-unit retry
    if indices is None:
        uobj.put(value, path)
    else:
        uobj.put(value, path, indices)


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
        topology_only: bool = False,
        initial_gro_path: Optional[str] = None,
        trajectory_path: Optional[str] = None,
        energy_path: Optional[str] = None,
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
        topology_only : when True, the resulting UDF contains the topology
                        (Set_of_Molecules / Molecular_Attributes / Interactions)
                        but **no** Structure record. Useful when the user
                        loads coordinates from a separate .gro/.xtc and energy
                        from a .xvg directly in J-OCTA Viewer.
        """
        raw = TopParser().parse(top_path)
        model = TopAdapter().build(raw, gro_path, mdp_path=mdp_path)

        # Resolve frames:
        #   --trajectory                : multi-frame trajectory from .gro/.xtc
        #   --topology-only             : single initial frame from initial_gro
        #                                 (or the positional gro_path)
        #   default                     : the 1 frame already loaded into
        #                                 model.frames from gro_path
        if trajectory_path is not None:
            frames: Optional[List[GROFrameData]] = _load_trajectory_frames(
                trajectory_path, gro_path)
            energy_times, energy_series = (
                _load_energy_series(energy_path) if energy_path else (None, None)
            )
        elif topology_only:
            if initial_gro_path is not None and initial_gro_path != gro_path:
                # Read a single frame from the override gro so we never
                # depend on the positional gro_path's coordinates here.
                initial_model = TopAdapter().build(
                    raw, initial_gro_path, mdp_path=mdp_path)
                frames = list(initial_model.frames[:1])
            else:
                frames = list(model.frames[:1])
            energy_series = None
            energy_times = None
        else:
            frames = None
            energy_series = None
            energy_times = None

        self.export_model(model, template_path, out_path,
                          frames=frames,
                          cognac_version=cognac_version,
                          energy_times=energy_times,
                          energy_series=energy_series)

    def export_model(
        self,
        model: TopModel,
        template_path: str,
        out_path: str,
        frames: Optional[List[GROFrameData]] = None,
        cognac_version: Optional[str] = None,
        energy_times: Optional[List[float]] = None,
        energy_series: Optional[dict] = None,
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
        # Pre-aggregate per-frame statistics (energy + temperature + pressure
        # + density + volume) when an .xvg series was given.
        per_frame_stats = _aggregate_statistics_per_frame(
            frames_to_write, energy_times, energy_series,
        )
        for i, frame in enumerate(frames_to_write):
            with _section(f"Structure[record={i}]", template_path, out_path):
                self._append_structure(
                    uobj, model, frame,
                    energy_values=per_frame_stats[i] if per_frame_stats else None,
                )

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
                # Atom_ID is per-molecule local 0-indexed (matches the COGNAC
                # convention used by cognac-shipped sample BDFs; using a
                # global counter here made J-OCTA's atom-table display the
                # last molecule's offset (1617..1649) which looked off).
                uobj.put(iatom,
                         "Set_of_Molecules.molecule[].atom[].Atom_ID",
                         [imol, iatom])
                # Atom_Name と Atom_Type_Name の役割分担 (2026-05-27 最終):
                # - Atom_Name      = element symbol (`C`, `O`, `H`)
                #   ABINIT-MP の read_pdb / obabel の xyz→pdb 変換で
                #   元素として解釈される。J-OCTA viewer の元素ベース描画
                #   (CPK 色, vdW 半径) もここから取得される。
                #   GROMACS 由来の force-field 型名 (例: "c30", "hc1") は
                #   obabel が atom type を `*` に変換し ABINIT-MP が
                #   "Atom type * is not supported" で失敗するため、
                #   Atom_Name には element symbol のみを使う。
                # - Atom_Type_Name = GROMACS top の atomtype 列の値
                #   = 現在使っている力場の型名:
                #     - OpenFF SMIRNOFF (interchange) → "MOL0_0", "MOL0_1", ...
                #       (per-atom unique な仮想 type)
                #     - GAFF / GAFF2 (antechamber / acpype 経由) → "c3", "oh",
                #       "hc", "h1", "os" 等
                #     - OPLS-AA → "opls_267" 等
                #   いずれも各 atom が unique な情報を持つので J-OCTA viewer
                #   の atom テーブルで type 別フィルタ等が可能。
                # - これにより `Atom_Type_Name` / `interaction_Site[].Type_Name`
                #   / `Molecular_Attributes.Atom_Type[].Name` が完全に同一値
                #   となり、UDF 内の atom type 参照整合性が自然に取れる。
                display_type = _display_type_name(atom.type_name, atom.element)
                uobj.put(atom.element,
                         "Set_of_Molecules.molecule[].atom[].Atom_Name",
                         [imol, iatom])
                uobj.put(display_type,
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

                uobj.put(display_type,
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
    def _append_structure(uobj, model: TopModel, frame: GROFrameData,
                          energy_values: Optional[dict] = None) -> None:
        """Port of append_structure (one GRO frame → one UDF dynamic record).

        When *energy_values* is given (mapping UDF field name -> float, e.g.
        ``{"Bond": 1596.5, "Angle": 2220.1, ...}``), those values land in
        ``Statistics_Data.Energy.Instantaneous.<field>`` of this record.
        """
        uobj.newRecord()

        mol_type_names = model.mol_type_names
        mol_specs = model.mol_specs
        ncount = 0

        for imol, mol_name in enumerate(model.mol_instance_list):
            mol_idx = mol_type_names.index(mol_name)
            spec = mol_specs[mol_idx]
            for iatom in range(len(spec.atoms)):
                posx, posy, posz = frame.coord_list[ncount]
                _put_with_unit_fallback(
                    uobj, posx, "Structure.Position.mol[].atom[].x",
                    [imol, iatom], "[nm]")
                _put_with_unit_fallback(
                    uobj, posy, "Structure.Position.mol[].atom[].y",
                    [imol, iatom], "[nm]")
                _put_with_unit_fallback(
                    uobj, posz, "Structure.Position.mol[].atom[].z",
                    [imol, iatom], "[nm]")
                ncount += 1

        cell = frame.cell
        _put_with_unit_fallback(uobj, cell[0],
                                "Structure.Unit_Cell.Cell_Size.a", None, "[nm]")
        _put_with_unit_fallback(uobj, cell[1],
                                "Structure.Unit_Cell.Cell_Size.b", None, "[nm]")
        _put_with_unit_fallback(uobj, cell[2],
                                "Structure.Unit_Cell.Cell_Size.c", None, "[nm]")
        uobj.put(90.0, "Structure.Unit_Cell.Cell_Size.alpha")
        uobj.put(90.0, "Structure.Unit_Cell.Cell_Size.beta")
        uobj.put(90.0, "Structure.Unit_Cell.Cell_Size.gamma")
        uobj.put(frame.step, "Steps")
        _put_with_unit_fallback(uobj, frame.time, "Time", None, "[ps]")

        # Statistics_Data.* (Energy / Temperature / Pressure / Density /
        # Volume). Default Hamiltonian to 0.0 (gmx energy doesn't expose it).
        uobj.put(0.0, "Statistics_Data.Energy.Instantaneous.Hamiltonian")
        if energy_values:
            for (rel_path, unit), value in energy_values.items():
                try:
                    _put_with_unit_fallback(
                        uobj, value,
                        f"Statistics_Data.{rel_path}",
                        None, unit,
                    )
                except Exception:
                    # Schema may lack this field on older cognac releases;
                    # skip silently to keep the rest of the record intact.
                    pass
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

        # --- Dynamics_Conditions.Time (from MDP) ---
        # delta_T native unit is [tau]; we pass [ps] and let UDFManager
        # convert via Unit_Parameter.
        time_loc = "Simulation_Conditions.Dynamics_Conditions.Time"
        if model.dt_ps > 0.0:
            _put_with_unit_fallback(uobj, model.dt_ps,
                                    time_loc + ".delta_T", None, "[ps]")
        if model.nsteps > 0:
            uobj.put(model.nsteps, time_loc + ".Total_Steps")
        # Prefer the trajectory output interval (nstxout-compressed),
        # falling back to nstenergy so frame-aligned energy plots line up.
        out_interval = (model.nstxout_compressed
                        if model.nstxout_compressed > 0
                        else model.nstenergy)
        if out_interval > 0:
            uobj.put(out_interval, time_loc + ".Output_Interval_Steps")

        uobj.write()

    @staticmethod
    def _write_molecular_attributes(uobj, model: TopModel) -> None:
        """Port of add_molecular_attributes_byTop."""
        uobj.jump(-1)

        # --- Atom_Type (unique types in order of first appearance) ---
        # Each entry maps a display name (e.g. "C4", "H19") to the mass of the
        # original OpenFF interchange type (MOL0_X). The display name matches
        # what we wrote in Set_of_Molecules.atom[].Atom_Type_Name and
        # interaction_Site[].Type_Name, preserving UDF cross-references.
        seen_types: List[tuple] = []  # list of (display_name, raw_type_name)
        seen_keys: set = set()
        for spec in model.mol_specs:
            for atom in spec.atoms:
                if atom.type_name not in seen_keys:
                    seen_keys.add(atom.type_name)
                    display = _display_type_name(atom.type_name, atom.element)
                    seen_types.append((display, atom.type_name))

        for ndata, (display, raw_type) in enumerate(seen_types):
            uobj.put(display,
                     "Molecular_Attributes.Atom_Type[].Name",
                     [ndata])
            mass = model.mass_dict.get(raw_type, 0.0)
            uobj.put(mass,
                     "Molecular_Attributes.Atom_Type[].Mass",
                     [ndata])

        # --- Bond_Potential (Harmonic) ---
        for j, bt in enumerate(model.bond_type_specs):
            uobj.put(bt.name,
                     "Molecular_Attributes.Bond_Potential[].Name", [j])
            uobj.put("Harmonic",
                     "Molecular_Attributes.Bond_Potential[].Potential_Type", [j])
            _put_with_unit_fallback(uobj, bt.r0,
                     "Molecular_Attributes.Bond_Potential[].R0", [j], "[nm]")
            _put_with_unit_fallback(uobj, bt.kb,
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
            _put_with_unit_fallback(uobj, at.k,
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
                _put_with_unit_fallback(uobj, PK,
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
                _put_with_unit_fallback(uobj, 1.0,
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
        # and build raw_type -> display_name lookup (e.g. MOL0_4 -> C4).
        used_types: set = set()
        display_map: dict = {}
        for spec in model.mol_specs:
            for atom in spec.atoms:
                used_types.add(atom.type_name)
                display_map[atom.type_name] = _display_type_name(
                    atom.type_name, atom.element)

        # Molecular_Attributes.Interaction_Site_Type
        ncount = 0
        for at in model.atom_type_specs:
            if at.name not in used_types:
                continue
            rval = at.sigma * 1.5
            uobj.put(display_map.get(at.name, at.name),
                     "Molecular_Attributes.Interaction_Site_Type[].Name",
                     [ncount])
            uobj.put(1,
                     "Molecular_Attributes.Interaction_Site_Type[].Num_of_Atoms",
                     [ncount])
            _put_with_unit_fallback(uobj, rval,
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
                d1 = display_map.get(at1.name, at1.name)
                d2 = display_map.get(at2.name, at2.name)
                name = "{}-{}".format(d1, d2)

                uobj.put(name,            location + ".Name",           [ncount])
                uobj.put("Lennard_Jones", location + ".Potential_Type", [ncount])
                uobj.put(d1,              location + ".Site1_Name",     [ncount])
                uobj.put(d2,              location + ".Site2_Name",     [ncount])
                _put_with_unit_fallback(uobj, cutoff,        location + ".Cutoff",         [ncount], "[nm]")
                uobj.put(scale,         location + ".Scale_1_4_Pair", [ncount])
                _put_with_unit_fallback(uobj, sigma,         location + ".Lennard_Jones.sigma",
                         [ncount], "[nm]")
                _put_with_unit_fallback(uobj, epsilon,       location + ".Lennard_Jones.epsilon",
                         [ncount], "[kJ/mol]")
                ncount += 1

        uobj.write()
