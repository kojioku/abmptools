# -*- coding: utf-8 -*-
"""Which force field a GROMACS topology came from, for ``--ff auto``.

``Unit_Parameter.Comment`` carries ``FF=n`` (see :mod:`.udf_writer`), and a
``.top`` does not say which force field it was built with. This module looks
at what the topology does say and names one -- or says why it cannot.

Evidence, strongest first:

1. ``#include`` names: ``amber14sb.ff/forcefield.itp``, ``oplsaa.ff/...``,
   ``charmm36.ff/...``, an ``.itp`` with ``gaff2`` / ``gaff`` in its name.
2. ``[ atomtypes ]``: lower-case GAFF-style names (``c3``, ``hc``) are told
   apart as GAFF or GAFF2 by their Lennard-Jones values, which differ for
   77 of the 80 types the two share (GAFF 1.81 / GAFF2 2.2.20); upper-case
   AMBER protein types (``CT``, ``HC``, ``N3``); ``opls_123``.
3. ``[ defaults ]``: combination rule and fudge factors -- AMBER-family
   (2, 0.5, 0.8333), OPLS (3, 0.5, 0.5), CHARMM (2, 1.0, 1.0). Used to
   rule out, not to name.

Only the types the molecules use count, and not those of water and ions
(molecules of only O and H with at most four atoms, and single-atom
molecules): every force field brings them.

Nothing is guessed. Two families (a protein in AMBER types and a ligand in
GAFF2), a family without an ``FF=n`` number (CHARMM), types that match no
table, or ``[ defaults ]`` that contradict the types: no name, and the
reasons say which. ``amber20`` is never chosen here; the topology does not
tell it from ``amber``.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set

#: Lennard-Jones (sigma nm, epsilon kJ/mol) of the types whose GAFF and GAFF2
#: values differ by more than 2 % -- gaff.dat 1.81 and gaff2.dat 2.2.20.
_GAFF_LJ = {"br": (0.359923, 1.75728), "c": (0.339967, 0.359824), "c1": (0.339967, 0.87864), "c2": (0.339967, 0.359824), "ca": (0.339967, 0.359824), "cc": (0.339967, 0.359824), "cd": (0.339967, 0.359824), "ce": (0.339967, 0.359824), "cf": (0.339967, 0.359824), "cg": (0.339967, 0.87864), "ch": (0.339967, 0.87864), "cp": (0.339967, 0.359824), "cq": (0.339967, 0.359824), "cu": (0.339967, 0.359824), "cv": (0.339967, 0.359824), "cx": (0.339967, 0.359824), "cy": (0.339967, 0.359824), "cz": (0.339967, 0.359824), "f": (0.311815, 0.255224), "h1": (0.247135, 0.065689), "h2": (0.229317, 0.065689), "h3": (0.211499, 0.065689), "h4": (0.251055, 0.06276), "h5": (0.242146, 0.06276), "ha": (0.259964, 0.06276), "hc": (0.264953, 0.065689), "hn": (0.106908, 0.065689), "hp": (0.106908, 0.065689), "hs": (0.106908, 0.065689), "hx": (0.195998, 0.065689), "n": (0.325, 0.71128), "n1": (0.325, 0.71128), "n2": (0.325, 0.71128), "n3": (0.325, 0.71128), "n4": (0.325, 0.71128), "na": (0.325, 0.71128), "nb": (0.325, 0.71128), "nc": (0.325, 0.71128), "nd": (0.325, 0.71128), "ne": (0.325, 0.71128), "nf": (0.325, 0.71128), "nh": (0.325, 0.71128), "ni": (0.325, 0.71128), "nj": (0.325, 0.71128), "nk": (0.325, 0.71128), "nl": (0.325, 0.71128), "nm": (0.325, 0.71128), "nn": (0.325, 0.71128), "no": (0.325, 0.71128), "np": (0.325, 0.71128), "nq": (0.325, 0.71128), "o": (0.295992, 0.87864), "oh": (0.306647, 0.880314), "op": (0.300001, 0.71128), "oq": (0.300001, 0.71128), "os": (0.300001, 0.71128), "p2": (0.374177, 0.8368), "p3": (0.374177, 0.8368), "p4": (0.374177, 0.8368), "p5": (0.374177, 0.8368), "pb": (0.374177, 0.8368), "pc": (0.374177, 0.8368), "pd": (0.374177, 0.8368), "pe": (0.374177, 0.8368), "pf": (0.374177, 0.8368), "px": (0.374177, 0.8368), "py": (0.374177, 0.8368), "s": (0.356359, 1.046), "s2": (0.356359, 1.046), "s4": (0.356359, 1.046), "s6": (0.356359, 1.046), "sh": (0.356359, 1.046), "sp": (0.356359, 1.046), "sq": (0.356359, 1.046), "ss": (0.356359, 1.046), "sx": (0.356359, 1.046), "sy": (0.356359, 1.046)}
_GAFF2_LJ = {"br": (0.361259, 1.645149), "c": (0.331521, 0.413379), "c1": (0.347896, 0.667766), "c2": (0.331521, 0.413379), "ca": (0.331521, 0.413379), "cc": (0.331521, 0.413379), "cd": (0.331521, 0.413379), "ce": (0.331521, 0.413379), "cf": (0.331521, 0.413379), "cg": (0.347896, 0.667766), "ch": (0.347896, 0.667766), "cp": (0.331521, 0.413379), "cq": (0.331521, 0.413379), "cu": (0.331521, 0.413379), "cv": (0.331521, 0.413379), "cx": (0.339771, 0.451035), "cy": (0.339771, 0.451035), "cz": (0.331521, 0.413379), "f": (0.303422, 0.348109), "h1": (0.2422, 0.087027), "h2": (0.224382, 0.087027), "h3": (0.206564, 0.087027), "h4": (0.253639, 0.067362), "h5": (0.24473, 0.067362), "ha": (0.262548, 0.067362), "hc": (0.260018, 0.087027), "hn": (0.11065, 0.04184), "hp": (0.10746, 0.06025), "hs": (0.108903, 0.051882), "hx": (0.188746, 0.087027), "n": (0.318086, 0.684502), "n1": (0.327352, 0.459403), "n2": (0.338417, 0.393714), "n3": (0.33651, 0.358987), "n4": (0.249951, 16.212163), "na": (0.320581, 0.854373), "nb": (0.338417, 0.393714), "nc": (0.338417, 0.393714), "nd": (0.338417, 0.393714), "ne": (0.338417, 0.393714), "nf": (0.338417, 0.393714), "nh": (0.318995, 0.89956), "ni": (0.318086, 0.684502), "nj": (0.318086, 0.684502), "nk": (0.25886, 10.649535), "nl": (0.25886, 10.649535), "nm": (0.318995, 0.89956), "nn": (0.318995, 0.89956), "no": (0.33651, 0.358987), "np": (0.33651, 0.358987), "nq": (0.33651, 0.358987), "o": (0.304812, 0.612119), "oh": (0.324287, 0.389112), "op": (0.31561, 0.303758), "oq": (0.31561, 0.303758), "os": (0.31561, 0.303758), "p2": (0.369402, 0.960228), "p3": (0.369402, 0.960228), "p4": (0.369402, 0.960228), "p5": (0.369402, 0.960228), "pb": (0.369402, 0.960228), "pc": (0.369402, 0.960228), "pd": (0.369402, 0.960228), "pe": (0.369402, 0.960228), "pf": (0.369402, 0.960228), "px": (0.369402, 0.960228), "py": (0.369402, 0.960228), "s": (0.353241, 1.181562), "s2": (0.353241, 1.181562), "s4": (0.353241, 1.181562), "s6": (0.40584, 0.256898), "sh": (0.353241, 1.181562), "sp": (0.353241, 1.181562), "sq": (0.353241, 1.181562), "ss": (0.353241, 1.181562), "sx": (0.353241, 1.181562), "sy": (0.40584, 0.256898)}

#: Every GAFF / GAFF2 type name. Those not in the tables above (c3, cl, ho,
#: i -- close or zero in both -- and the ones only one of them has) are
#: GAFF-family but do not tell GAFF from GAFF2.
_GAFF_NAMES = frozenset("""
br c c1 c2 c3 c5 c6 ca cc cd ce cf cg ch cl cp cq cs cu cv cx cy cz f h1 h2 h3 h4
h5 ha hc hn ho hp hs hw hx i n n+ n1 n2 n3 n4 n5 n6 n7 n8 n9 na nb nc nd ne nf nh
ni nj nk nl nm nn no np nq ns nt nu nv nx ny nz o oh op oq os ow p2 p3 p4 p5 pb pc
pd pe pf px py s s2 s4 s6 sh sp sq ss sx sy
""".split())

#: Relative tolerance on sigma and epsilon (topologies print 5-6 digits).
_LJ_RTOL = 0.005

#: Upper-case atom types of the AMBER protein / nucleic force fields
#: (parm99 / ff14SB families). A topology made from them prints these names.
_AMBER_TYPES = frozenset("""
C CA CB CC CD CK CM CN CQ CR CT CV CW C* CX CS 2C 3C C4 C5 C8 CO CI
H H1 H2 H3 H4 H5 HA HC HO HP HS HZ HW
N N2 N3 N* NA NB NC NT NY
O O2 OH OS OW P S SH
""".split())

#: The Lennard-Jones (sigma nm, epsilon kJ/mol) pairs of parm10.dat (parm99 +
#: ff99SB + bsc0 + OL3; ff14SB keeps them). An upper-case type counts as
#: AMBER only with one of these: ParmEd renames types (C1, O1, H1 ...) and a
#: name alone would claim them.
_AMBER_LJ = [(0.0, 0.0), (0.106908, 0.065689), (0.141225, 3.743425), (0.195998, 0.065689), (0.211499, 0.065689), (0.229317, 0.065689), (0.242146, 0.06276), (0.247135, 0.065689), (0.251055, 0.06276), (0.259964, 0.06276), (0.264953, 0.065689), (0.295992, 0.87864), (0.300001, 0.71128), (0.30524, 1.923757), (0.306647, 0.880314), (0.311815, 0.255224), (0.325, 0.71128), (0.329633, 0.71128), (0.339967, 0.359824), (0.339967, 0.45773), (0.347094, 1.10876), (0.356359, 1.046), (0.374177, 0.8368), (0.395559, 1.33888), (0.418722, 1.6736)]

#: [ defaults ] (comb-rule, fudgeLJ, fudgeQQ) of each family.
_DEFAULTS = {"amber": (2, 0.5, 0.8333), "opls": (3, 0.5, 0.5), "charmm": (2, 1.0, 1.0)}

#: Family -> the gro2udf force-field name it is written as (None: no FF=n number).
_NAME = {"amber": "amber", "gaff": "gaff", "gaff2": "gaff2", "opls": "oplsaa", "charmm": None,
         "gaff or gaff2": None}


@dataclass
class FFDecision:
    """What ``--ff auto`` decided. ``name`` is a gro2udf force-field name
    (``gaff2``, ``amber`` ...) or ``None`` when no single one fits."""
    name: Optional[str]
    families: Dict[str, List[str]] = field(default_factory=dict)   # family -> evidence
    reasons: List[str] = field(default_factory=list)

    def summary(self) -> str:
        ev = "; ".join("%s: %s" % (f, ", ".join(e[:4]) + (" ..." if len(e) > 4 else ""))
                       for f, e in sorted(self.families.items()))
        return ev or "no evidence"


def include_names(top_path: str, _seen: Optional[Set[str]] = None) -> List[str]:
    """The ``#include`` names of *top_path* and of the files it includes that
    can be read (relative to the including file). A force-field directory
    (``amber14sb.ff/forcefield.itp``) usually lives in GROMACS' own library
    and is not read -- its name is the evidence."""
    seen = _seen if _seen is not None else set()
    path = os.path.abspath(top_path)
    if path in seen:
        return []
    seen.add(path)
    out: List[str] = []
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return out
    for line in lines:
        s = line.strip()
        if not s.startswith("#include"):
            continue
        parts = s.split(None, 1)
        if len(parts) < 2:
            continue
        name = parts[1].strip().strip('"<>\'')
        out.append(name)
        sub = os.path.join(os.path.dirname(path), name)
        if os.path.isfile(sub):
            out.extend(include_names(sub, seen))
    return out


def _family_of_include(name: str) -> Optional[str]:
    n = os.path.basename(os.path.dirname(name)) + "/" + os.path.basename(name)
    n = n.lower()
    if re.search(r"\bamber[\w\-]*\.ff/", n):
        return "amber"
    if "oplsaa.ff/" in n:
        return "opls"
    if re.search(r"charmm[\w\-]*\.ff/", n):
        return "charmm"
    base = os.path.basename(name).lower()
    if "gaff2" in base:
        return "gaff2"
    if "gaff" in base:
        return "gaff"
    return None


def _close(a, b) -> bool:
    return all(abs(x - y) <= _LJ_RTOL * max(abs(x), abs(y), 1e-12) for x, y in zip(a, b))


def _solvent_types(model) -> Set[str]:
    """Types used only by water and single-atom molecules."""
    solvent, other = set(), set()
    for mol in model.mol_specs:
        els = {a.element for a in mol.atoms}
        is_solvent = len(mol.atoms) == 1 or (len(mol.atoms) <= 4 and els <= {"O", "H", ""})
        (solvent if is_solvent else other).update(a.type_name for a in mol.atoms)
    return solvent - other


def detect(model, includes: Sequence[str] = ()) -> FFDecision:
    """Name the force field of *model* (a :class:`TopModel`), or say why not."""
    fam: Dict[str, List[str]] = {}

    def add(f, ev):
        fam.setdefault(f, [])
        if ev not in fam[f]:
            fam[f].append(ev)

    for inc in includes:
        f = _family_of_include(inc)
        if f:
            add(f, "#include " + inc)

    # Only the types the molecules use: a force-field directory's
    # ffnonbonded.itp lists every type it has, used or not.
    skip = _solvent_types(model)
    used = {a.type_name for mol in model.mol_specs for a in mol.atoms}
    unmatched: List[str] = []
    neutral: List[str] = []
    for t in model.atom_type_specs:
        name = t.name
        if name in skip or name not in used:
            continue
        if re.match(r"^opls_\d+$", name):
            add("opls", name)
        elif name in _AMBER_TYPES and any(
                _close((t.sigma, t.epsilon), lj) if lj[1] else t.epsilon == 0.0
                for lj in _AMBER_LJ):
            add("amber", name)
        elif name in _GAFF_LJ:
            lj = (t.sigma, t.epsilon)
            if _close(lj, _GAFF2_LJ[name]):
                add("gaff2", name)
            elif _close(lj, _GAFF_LJ[name]):
                add("gaff", name)
            else:
                unmatched.append(name)
        elif name in _GAFF_NAMES:
            neutral.append(name)
        else:
            unmatched.append(name)
    if neutral and not ({"gaff", "gaff2"} & set(fam)):
        fam["gaff or gaff2"] = neutral

    reasons: List[str] = []
    if unmatched:
        reasons.append("atom types no table names: %s%s"
                       % (", ".join(unmatched[:8]), " ..." if len(unmatched) > 8 else ""))
    d = (model.comb_rule, round(model.fudge_lj, 4), round(model.fudge_qq, 4))
    by_defaults = [f for f, v in _DEFAULTS.items()
                   if v[0] == d[0] and abs(v[1] - d[1]) < 1e-3 and abs(v[2] - d[2]) < 2e-3]

    families = sorted(fam)
    if not families:
        reasons.append("nothing in the topology names a force field")
        return FFDecision(None, fam, reasons)
    if len(families) > 1:
        reasons.insert(0, "more than one force field: %s" % ", ".join(families))
        return FFDecision(None, fam, reasons)
    f = families[0]
    group = "amber" if f in ("amber", "gaff", "gaff2", "gaff or gaff2") else f
    if by_defaults and group not in by_defaults:
        reasons.insert(0, "[ defaults ] %s fit %s, not %s" % (d, by_defaults[0], f))
        return FFDecision(None, fam, reasons)
    if unmatched:
        reasons.insert(0, "%s types, but also types it cannot place" % f)
        return FFDecision(None, fam, reasons)
    if f == "gaff or gaff2":
        reasons.insert(0, "GAFF-family types, but none that tells GAFF from GAFF2 (%s)"
                       % ", ".join(fam[f][:6]))
        return FFDecision(None, fam, reasons)
    if _NAME[f] is None:
        reasons.insert(0, "%s has no FF=n number" % f)
        return FFDecision(None, fam, reasons)
    return FFDecision(_NAME[f], fam, reasons)
