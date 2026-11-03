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

import logging
import re
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from ..core.top_atomtypes import fold_atomtypes
from .guard import scan_dihedral_functs, scan_sections

logger = logging.getLogger(__name__)

#: ``[ section ]`` header.
_SECTION_RE = re.compile(r"^\s*\[\s*([A-Za-z_0-9]+)\s*\]")


#: Sections that belong to the ``[ moleculetype ]`` above them. Anything else
#: is file-wide ([ molecules ], [ system ], [ atomtypes ] ...) and is read
#: even while a molecule is being skipped.
_MOLECULE_SCOPED = frozenset({
    "atoms", "bonds", "angles", "dihedrals", "constraints", "exclusions",
    "pairs", "settles", "position_restraints", "distance_restraints",
    "dihedral_restraints", "angle_restraints", "orientation_restraints",
    "virtual_sites1", "virtual_sites2", "virtual_sites3", "virtual_sites4",
    "virtual_sitesn", "dummies1", "dummies2", "dummies3", "dummies4", "cmap",
})


def _molecules_listed(lines) -> "Optional[set]":
    """``[ molecules ]`` に並ぶ分子の名前。

    その節が無ければ ``None`` (単体の ``.itp`` などで「分からない」の意)。
    """
    names = set()
    seen = False
    in_section = False
    for line in lines:
        s = line.strip()
        m = _SECTION_RE.match(s)
        if m:
            in_section = m.group(1).lower() == "molecules"
            seen = seen or in_section
            continue
        if not in_section:
            continue
        s = s.split(";")[0].strip()
        if s and not s.startswith("#"):
            names.add(s.split()[0])
    return names if seen else None


def _is_section(stripped: str, name: str) -> bool:
    """*stripped* が ``[ name ]`` の見出しか。

    **``[x]`` と ``[ x ]`` を同じに扱う。**GROMACS はどちらも読むので、
    片方しか見ないと、その書き方の ``.itp`` の中身が黙って落ちる
    (Martini の脂質ライブラリは空白なしで書かれていて、112 個の
    ``[moleculetype]`` が 1 つも読まれていなかった)。
    """
    m = _SECTION_RE.match(stripped)
    return bool(m) and m.group(1).lower() == name


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
    #: ``[ nonbond_params ]``: ``{(type1, type2): (V, W)}``, the pair written
    #: with the smaller name first. With ``comb-rule`` 2 or 3 the two numbers
    #: are sigma [nm] and epsilon [kJ/mol]; with 1 they are c6 and c12.
    #: A force field that keeps all of its LJ here (coarse-grained ones do)
    #: leaves ``[ atomtypes ]`` at zero, so without this the pairs are lost.
    nonbond_params: Dict = field(default_factory=dict)
    #: ``[ moleculetype ]`` の 2 つめの欄: ``{分子名: nrexcl}``。
    #: **どこまでの隣を非結合から外すか**。読まないと 1-3 の扱いが力場に
    #: 関わらず固定になる (全原子の「1-2 と 1-3 を外す」流儀のまま)。
    nrexcl: Dict = field(default_factory=dict)
    #: ``[ exclusions ]``: ``{molecule: [(i, j), ...]}`` (1-based, i < j).
    #: One line lists several partners, so a line is several pairs.
    exclusions: Dict = field(default_factory=dict)
    #: Phases [deg] of the proper dihedrals that are written as the engine's
    #: Amber form. It only takes 0 or 180, so anything else has to be caught.
    torsion_phases: List = field(default_factory=list)
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

    # section name (lower-case) -> number of data lines, after #include
    # resolution. Used by :mod:`.guard` to spot content that is not read.
    sections: Dict[str, int] = field(default_factory=dict)

    # every ``funct`` seen in [ dihedrals ] / [ dihedraltypes ], read from the
    # raw text because an unsupported funct is lost during parsing.
    dihedral_functs: List[int] = field(default_factory=list)

    # [ constraints ] lines read as harmonic bonds (TopParser(
    # constraints_as_bonds=k)). 0 when the option is off.
    n_constraints_as_bonds: int = 0

    # moleculetype name -> [(ai, aj), ...] (1-based) from [ constraints ].
    # Always read, whether or not they become bonds: they are connectivity,
    # which --make-whole needs to put a constraint-only molecule together.
    constraint_pairs: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)


#: Default force constant for --constraints-as-bonds, kJ/mol/nm^2. A
#: constraint has none, so any value is a stand-in; this one is stiff enough
#: that the bond stays at its constrained length to within ~0.01 nm at
#: room-temperature thermal energy, which is what matters for viewing.
DEFAULT_CONSTRAINT_K = 50000.0


# ---------------------------------------------------------------------------
# Pure helper functions (module level)
# ---------------------------------------------------------------------------

def _is_comment(line: str) -> bool:
    s = line.strip()
    return bool(s) and (s[0] == ";" or s[0] == "#")


def _is_blank(line: str) -> bool:
    return line.strip() == ""


def _fields(line: str) -> List[str]:
    """The whitespace-separated fields of a topology line, with the
    end-of-line comment (from the first ``;``) dropped. Every section is
    split through this: a comment can follow the data on any line (Martini 3
    writes ``1 2 ; sigma-epsilon format of LJ parameters`` under
    ``[ defaults ]``), and a field count or a positional float read from
    ``line.split()`` then picks up a word of the comment."""
    return line.split(";", 1)[0].split()


def get_bond_name(type1: str, type2: str) -> str:
    """2つの原子名からソート済みの結合名を生成する。"""
    if type1 <= type2:
        return type1 + "-" + type2
    return type2 + "-" + type1


def get_angle_name(type1: str, type2: str, type3: str) -> str:
    """3つの原子名からソート済みの角度名を生成する。"""
    if type1 <= type3:
        return type1 + "-" + type2 + "-" + type3
    return type3 + "-" + type2 + "-" + type1


def get_torsion_name(type1: str, type2: str, type3: str, type4: str) -> str:
    """4つの原子名からソート済みの二面角名を生成する。"""
    if type1 < type4:
        return type1 + "-" + type2 + "-" + type3 + "-" + type4
    elif type1 == type4:
        if type2 <= type3:
            return type1 + "-" + type2 + "-" + type3 + "-" + type4
        return type4 + "-" + type3 + "-" + type2 + "-" + type1
    return type4 + "-" + type3 + "-" + type2 + "-" + type1


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
    """Parse a GROMACS .top file (with resolved #include ITP) into TopRawData.

    Parameters
    ----------
    constraints_as_bonds : float, optional
        When given, every ``[ constraints ]`` line inside a ``[ moleculetype ]``
        is read as a harmonic bond (funct 1) of the constrained length with
        this force constant (kJ/mol/nm^2). Without it ``[ constraints ]`` is
        not read, and a molecule held together only by constraints (common in
        Martini: rings, small rigid fragments) ends up with **no bonds** in the
        UDF. A constraint has no force constant, so the value is a stand-in --
        the connectivity and the length are what carry over exactly.
    """

    def __init__(self, constraints_as_bonds: Optional[float] = None) -> None:
        self.constraints_as_bonds = constraints_as_bonds

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def parse(self, top_path: str) -> TopRawData:
        """Parse *top_path* and return a :class:`TopRawData`."""
        lines = self._resolve_includes(top_path)
        lines = self._fold_atomtypes(lines)
        raw = TopRawData()
        raw.sections = scan_sections(lines)
        raw.dihedral_functs = scan_dihedral_functs(lines)

        # Defaults
        raw.comb_rule, raw.fudge_lj, raw.fudge_qq = self._parse_defaults(lines)

        # Global type definitions (atomtypes, bondtypes, …)
        at, bt, agt, tt = self._parse_definitions(lines)
        raw.atomtypes = at
        raw.nonbond_params = self._parse_nonbond_params(lines)
        # The explicit pairs are read now, so they are no longer dropped.
        if raw.nonbond_params:
            raw.sections["nonbond_params"] = 0
        raw.bondtypes = bt
        raw.angletypes = agt
        raw.torsiontypes = tt

        # Per-molecule topology + molecule instances.
        # The names in [ molecules ] are read first so that the molecules the
        # system never contains can be skipped: a force field shipped as one
        # .itp library carries hundreds of them, and reading those only costs
        # time and raises features (functs, sections) the system does not use.
        self._used_molecules = _molecules_listed(lines)
        (raw.mol_types,
         raw.atomlist,
         raw.bondlist,
         raw.anglelist,
         raw.torsionlist,
         raw.mol_instance_list,
         raw.bond_types_from_mol,
         raw.angle_types_from_mol,
         raw.torsion_types_from_mol,
         raw.nrexcl,
         raw.exclusions) = self._parse_topology(lines)
        # [a1, a2, a3, a4, funct, improper, [PHASE, PK, PN]]
        raw.torsion_phases = [
            tt[6][0] for tt in raw.torsion_types_from_mol
            if len(tt) > 6 and tt[4] in (1, 9, 4) and tt[6]]
        raw.n_constraints_as_bonds = self._n_constraints_as_bonds
        raw.constraint_pairs = self._constraint_pairs

        # The constraints have been carried into the bond list, so nothing in
        # that section is lost any more: take it off the guard's list.
        if raw.n_constraints_as_bonds:
            raw.sections["constraints"] = 0

        return raw

    # ------------------------------------------------------------------
    # [ nonbond_params ]
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_nonbond_params(lines: List[str]) -> Dict:
        """``[ nonbond_params ]`` -> ``{(t1, t2): (V, W)}`` (keys sorted).

        Each line is ``type1 type2 funct V W``. Only ``funct = 1`` (LJ) is
        taken; a different funct would mean another functional form, and
        guessing one is worse than leaving the pair to the combination rule.
        """
        out: Dict = {}
        in_section = False
        for line in lines:
            s = line.split(";")[0].strip()
            m = _SECTION_RE.match(line)
            if m:
                in_section = m.group(1).strip().lower() == "nonbond_params"
                continue
            if not in_section or not s or s.startswith("#"):
                continue
            f = s.split()
            if len(f) < 5:
                continue
            try:
                funct = int(f[2])
                v, w = float(f[3]), float(f[4])
            except ValueError:
                continue
            if funct != 1:
                continue
            key = tuple(sorted((f[0], f[1])))
            out[key] = (v, w)
        return out

    # ------------------------------------------------------------------
    # Step 1: resolve #include
    # ------------------------------------------------------------------

    def _resolve_includes(self, top_path: str) -> List[str]:
        """Read *top_path* and inline all ``#include "*.itp"`` files."""
        dirpath = os.path.dirname(os.path.abspath(top_path))
        with open(top_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        return self._inline_includes(lines, dirpath)

    def _fold_atomtypes(self, lines: List[str]) -> List[str]:
        """Interchange が原子ごとに吐く atomtype を、 読み込み時に畳む。

        SMIRNOFF には atom type の概念が無いので、 openff-interchange は
        **原子 1 個につき 1 型** (``MOL0_0``, ``MOL0_1``, ...) を書く。
        パラメータが同じものまで別型になるため、 PVA 10-mer では 75 型が
        出るが中身は 5 種類しかない。 これをそのまま UDF に持ち込むと
        受け取る側が扱えない (原子タイプ欄が意味のない名前で埋まる)。

        2.16 以降の :mod:`abmptools.amorphous` は ``.top`` を書く時点で
        畳むが、 **それ以前に組んだ ``.top`` が既に手元にある**。 ここで
        読み込み時にも畳んで、 作り直さずに UDF 化できるようにする。
        既に畳まれた ``.top`` では何も起きない (no-op)。

        畳み込み自体は :func:`abmptools.core.top_atomtypes.fold_atomtypes`
        が行う。 型名で引くセクション (``bondtypes`` 等) を持つ ``.top``
        は、 引き先が変わって**エラーを出さずに力場が変わる**ため、
        向こう側で素通しされる。
        """
        text, mapping = fold_atomtypes("".join(lines))
        if mapping:
            logger.info(
                "atomtypes folded on read: %d names -> %d types "
                "(openff-interchange writes one type per atom)",
                len(mapping), len(set(mapping.values())),
            )
            return text.splitlines(keepends=True)
        return lines

    def _inline_includes(self, lines: List[str], dirpath: str) -> List[str]:
        result: List[str] = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#include") and ".itp" in stripped:
                parts = _fields(stripped)
                if len(parts) >= 2:
                    itp_name = parts[1].strip('"\'')
                    itp_path = os.path.join(dirpath, itp_name)
                    try:
                        with open(itp_path, "r", encoding="utf-8") as f:
                            itp_lines = f.readlines()
                        result.append("\n")
                        result.extend(itp_lines)
                        logger.info("itp was included from %s", itp_path)
                        continue
                    except OSError as exc:
                        logger.warning("could not include %s: %s", itp_path, exc)
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
            if _is_section(stripped, "defaults"):
                in_section = True
                continue
            if in_section:
                if _is_comment(line) or _is_blank(line):
                    continue
                parts = _fields(line)
                if not parts:
                    continue
                # nbfunc comb_rule [gen-pairs [fudgeLJ [fudgeQQ]]]. The last
                # three are optional in GROMACS and default to no / 1.0 / 1.0
                # (Martini 3 writes only "1 2").
                comb_rule = int(parts[1])
                fudge_lj = float(parts[3]) if len(parts) > 3 else 1.0
                fudge_qq = float(parts[4]) if len(parts) > 4 else 1.0
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

            if _is_section(stripped, "atomtypes"):
                logger.info("reading atomtypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = True, False, False, False
                continue
            elif _is_section(stripped, "bondtypes"):
                logger.info("reading bondtypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = False, True, False, False
                continue
            elif _is_section(stripped, "angletypes"):
                logger.info("reading angletypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = False, False, True, False
                continue
            elif _is_section(stripped, "dihedraltypes"):
                logger.info("reading torsiontypes")
                flg_atom, flg_bond, flg_angle, flg_torsion = False, False, False, True
                continue
            elif stripped.startswith("[") and stripped.endswith("]"):
                # any other section → stop current
                flg_atom = flg_bond = flg_angle = flg_torsion = False
                continue

            if _is_comment(line) or _is_blank(line):
                continue

            parts = _fields(line)
            if not parts:
                continue

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
                    logger.warning("particle type %s detected", ptype)

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

        logger.info("%s atoms, %s bonds, %s angles, %s dihedrals "
                    "have been read (global defs)",
                    len(atomtypes), len(bondtypes),
                    len(angletypes), len(torsiontypes))
        return atomtypes, bondtypes, angletypes, torsiontypes

    # ------------------------------------------------------------------
    # Step 4: parse per-molecule topology
    # ------------------------------------------------------------------

    def _parse_topology(self, lines: List[str]):  # noqa: C901
        """
        Parse [ moleculetype ], [ atoms ], [ bonds ], [ angles ],
        [ dihedrals ], and [ molecules ].

        Returns a 9-tuple:
          (mol_types, atomlist, bondlist, anglelist, torsionlist,
           mol_instance_list,
           bond_types_from_mol, angle_types_from_mol, torsion_types_from_mol)
        """
        nrexcl_map = {}         # moleculetype -> nrexcl
        exclusion_map = {}      # moleculetype -> [(i, j), ...]
        flg_excl = False        # inside [ exclusions ] of a kept molecule
        skip_molecule = False   # inside a molecule this system does not use
        flg1 = False  # moleculetype
        flg2 = False  # atoms
        flg3 = False  # bonds
        flg4 = False  # angles
        flg5 = False  # dihedrals
        flg_con = False  # constraints
        flg_mol = False  # molecules
        self._n_constraints_as_bonds = 0
        self._constraint_pairs: Dict[str, List[Tuple[int, int]]] = {}

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
            if stripped.startswith("[") and stripped.endswith("]"):
                flg_con = False
                if skip_molecule:
                    sec = _SECTION_RE.match(stripped)
                    name = sec.group(1).lower() if sec else ""
                    if name in _MOLECULE_SCOPED:
                        # inside a molecule this system does not contain
                        flg1 = flg2 = flg3 = flg4 = flg5 = False
                        continue
                    if name != "moleculetype":
                        # a file-wide section ends the skipped molecule;
                        # [ molecules ] lives here and must not be skipped
                        skip_molecule = False
            if _is_section(stripped, "exclusions"):
                flg_excl = not skip_molecule
                flg1 = flg2 = flg3 = flg4 = flg5 = False
                continue
            if stripped.startswith("[") and stripped.endswith("]"):
                flg_excl = False
            if _is_section(stripped, "moleculetype"):
                logger.info("reading moleculetype")
                # Flush the molecule that just ended -- but only if it was
                # kept. Flushing a skipped one would push an empty entry into
                # atomlist and put every later molecule at the wrong index.
                if mol_types and not skip_molecule:
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
            elif _is_section(stripped, "atoms"):
                logger.info("reading atoms")
                flg1, flg2, flg3, flg4, flg5 = False, True, False, False, False
                continue
            elif _is_section(stripped, "bonds"):
                logger.info("reading bonds")
                flg1, flg2, flg3, flg4, flg5 = False, False, True, False, False
                continue
            elif _is_section(stripped, "angles"):
                logger.info("reading angles")
                flg1, flg2, flg3, flg4, flg5 = False, False, False, True, False
                continue
            elif _is_section(stripped, "dihedrals"):
                logger.info("reading dihedrals")
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
            elif _is_section(stripped, "molecules"):
                logger.info("reading molecules")
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
            elif _is_section(stripped, "constraints"):
                flg1 = flg2 = flg3 = flg4 = flg5 = False
                # Always read for connectivity (constraint_pairs); turned into
                # bonds only with constraints_as_bonds.
                flg_con = True
                if self.constraints_as_bonds is not None:
                    logger.info("reading constraints as bonds")
                continue
            elif stripped.startswith("[") and stripped.endswith("]"):
                # unknown section -- stop current flags
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
                elif flg_con:
                    flg_con = False
                elif flg_mol:
                    flg_mol = False
                continue

            stmp = _fields(line)
            if not stmp:
                continue

            if flg_excl:
                # A pair listed here must not interact. The engine has no
                # per-pair switch, but it does drop the non-bonded
                # interaction of a **bonded** pair, so each exclusion is
                # written as a bond whose force constant is zero: it removes
                # the pair and adds no energy of its own. This only works
                # while nrexcl is 1 (checked in the guard) -- with a larger
                # nrexcl the extra bond would exclude further neighbours too.
                try:
                    first = int(stmp[0])
                    partners = [int(x) for x in stmp[1:]]
                except ValueError:
                    continue
                here = exclusion_map.setdefault(mol_types[-1] if mol_types else "", [])
                for other in partners:
                    if other == first:
                        continue
                    pair = (min(first, other), max(first, other))
                    if pair in here:
                        continue
                    here.append(pair)
                    tidx = self._put_bond_type(
                        pair[0], pair[1], 1, "0.0", 0.0,
                        atomlist_mol, bondtypes_map, bond_types_from_mol)
                    bondlist_mol.append([pair[0], pair[1], tidx])
                continue

            if flg1:
                name = stmp[0]
                if len(stmp) > 1:
                    try:
                        nrexcl_map[name] = int(stmp[1])
                    except ValueError:
                        pass
                used = getattr(self, "_used_molecules", None)
                if used is not None and name not in used:
                    # Not in this system: skip its sections entirely.
                    skip_molecule = True
                    flg1 = flg2 = flg3 = flg4 = flg5 = False
                    continue
                skip_molecule = False
                mol_types.append(name)

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

            elif flg_con:
                # ai aj funct b0 -- funct 1 and 2 both fix the distance at b0
                # (2 only differs in not generating exclusions). Read as a
                # harmonic bond of that length.
                body = line.split(";", 1)[0].split()
                a1, a2 = int(body[0]), int(body[1])
                if mol_types:
                    self._constraint_pairs.setdefault(
                        mol_types[-1], []).append((a1, a2))
                if self.constraints_as_bonds is not None:
                    tidx = self._put_bond_type(
                        a1, a2, 1, body[3], self.constraints_as_bonds,
                        atomlist_mol, bondtypes_map, bond_types_from_mol)
                    bondlist_mol.append([a1, a2, tidx])
                    self._n_constraints_as_bonds += 1

            elif flg_mol:
                moltype = stmp[0]
                nmols   = int(stmp[1])
                for _ in range(nmols):
                    mol_instance_list.append(moltype)

        return (mol_types, atomlist, bondlist, anglelist, torsionlist,
                mol_instance_list,
                bond_types_from_mol, angle_types_from_mol,
                torsion_types_from_mol, nrexcl_map, exclusion_map)

    # ------------------------------------------------------------------
    # Internal helpers: put*Type -- deduplicate and register type entries
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
