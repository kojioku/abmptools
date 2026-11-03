# -*- coding: utf-8 -*-
"""
abmptools.gro2udf.guard
-----------------------
Detect ``.top`` content that gro2udf cannot carry into COGNAC.

gro2udf reads a fixed set of sections and emits a fixed set of COGNAC
potential types. Everything else used to be dropped **without a message**, so
the conversion appeared to succeed and produced a UDF that was quietly wrong.

The case that motivated this: a Martini 2 topology keeps every non-bonded
parameter in ``[ nonbond_params ]`` and leaves ``[ atomtypes ]`` at
``c6 = c12 = 0``. gro2udf ignores ``[ nonbond_params ]``, so the UDF came out
with **zero** ``Pair_Interaction`` entries, no warning, and an
``Interaction_Site_Type[].Range`` of 0 (it is derived as ``sigma * 1.5``).

Two severities, because the two failures are not equally recoverable:

fatal
    A term is **written with the wrong functional form or the wrong numbers**.
    The UDF looks complete, so nothing downstream can notice.
warning
    A term is **dropped**. The loss is real, but it is visible as an absence
    and this has been gro2udf's behaviour for every all-atom topology so far,
    so it must not start failing conversions that used to run.

This module only detects. Adding real support for these features is separate
work.
"""
from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

__all__ = [
    "UnsupportedTopFeatureError",
    "PARSED_SECTIONS",
    "FATAL_SECTIONS",
    "DROPPED_SECTIONS",
    "SUPPORTED_BOND_FUNCTS",
    "SUPPORTED_ANGLE_FUNCTS",
    "SUPPORTED_TORSION_FUNCTS",
    "PER_MOLECULE_SECTIONS",
    "used_moleculetypes",
    "scan_sections",
    "scan_dihedral_functs",
    "check_top",
    "raise_if_unsupported",
]


class UnsupportedTopFeatureError(ValueError):
    """The ``.top`` uses something gro2udf would mistranslate."""


#: Sections :class:`~.top_parser.TopParser` actually reads.
PARSED_SECTIONS = frozenset({
    "defaults", "atomtypes", "bondtypes", "angletypes", "dihedraltypes",
    "moleculetype", "atoms", "bonds", "angles", "dihedrals", "molecules",
})

#: Dropping these leaves a UDF that cannot describe the system at all.
FATAL_SECTIONS: Dict[str, str] = {
    # [ nonbond_params ] is read since 2026-11: TopParser puts it in
    # raw.nonbond_params and zeroes the section count, so it only reaches
    # this table when the parser could not use it (a funct other than 1).
    "nonbond_params":
        "explicit non-bonded pairs written with a funct this reader does not "
        "know (only funct 1, plain LJ, is read). A field that keeps all of "
        "its LJ here would end up with no Pair_Interaction entries at all",
    "constraints":
        "constrained bonds. They carry no force constant to fall back on, so "
        "the atoms would end up connected by nothing (pass "
        "--constraints-as-bonds to read them as bonds of the same length)",
}

#: Dropping these loses real terms, but the loss is an absence rather than a
#: wrong number, and all-atom topologies have always been converted this way.
DROPPED_SECTIONS: Dict[str, str] = {
    "pairtypes": "per-type 1-4 LJ overrides. The plain [ pairs ] list is\n                  fine -- COGNAC reproduces it from Scale_1_4_Pair, which\n                  gro2udf sets from fudgeLJ -- but a [ pairtypes ] override\n                  cannot be expressed by that single scale factor",
    "settles": "rigid water geometry (note that #ifdef is not evaluated, so a "
               "FLEXIBLE branch in the same .itp is read instead)",
    "virtual_sites1": "virtual sites",
    "virtual_sites2": "virtual sites",
    "virtual_sites3": "virtual sites",
    "virtual_sites4": "virtual sites",
    "virtual_sitesn": "virtual sites",
    "dummies2": "virtual (dummy) sites",
    "dummies3": "virtual (dummy) sites",
    "dummies4": "virtual (dummy) sites",
    "cmaptypes": "the CMAP correction",
    "cmap": "the CMAP correction",
}

# ``#ifdef`` is not evaluated anywhere in gro2udf, so a section inside an
# inactive branch still counts here. Sections that are normally ifdef-gated
# (the ``*_restraints`` family) are therefore left out of the list above:
# warning about them would be wrong more often than right. ``settles`` is
# kept even though it sits in an ``#else`` branch, because that branch is the
# active one -- the ``#ifdef FLEXIBLE`` bonds are read in its place.

#: ``top_exporter`` writes every bond as COGNAC ``Harmonic``, whatever the funct.
SUPPORTED_BOND_FUNCTS = frozenset({1})
#: funct 1 -> COGNAC ``Theta``, funct 2 -> ``Cosine`` (same form, same units).
SUPPORTED_ANGLE_FUNCTS = frozenset({1, 2, 10})
#: funct 1/9/4 -> ``Amber``, 3 -> ``Cosine_Polynomial``; others are skipped.
SUPPORTED_TORSION_FUNCTS = frozenset({1, 3, 4, 9})

_FUNCT_NOTE = {
}
_SECTION_RE = re.compile(r"^\s*\[\s*([A-Za-z_0-9]+)\s*\]")


def _is_data_line(line: str) -> bool:
    s = line.strip()
    return bool(s) and not s.startswith((";", "#", "*"))


#: Sections that live inside one ``[ moleculetype ]``. Lines in a molecule the
#: system never instantiates cost nothing, so they must not fail a conversion:
#: a force field shipped as one ``.itp`` library carries hundreds of molecules
#: the system does not use.
PER_MOLECULE_SECTIONS = frozenset({
    "atoms", "bonds", "angles", "dihedrals", "constraints", "exclusions",
    "pairs", "settles", "position_restraints", "distance_restraints",
    "dihedral_restraints", "angle_restraints", "orientation_restraints",
    "virtual_sites1", "virtual_sites2", "virtual_sites3", "virtual_sites4",
    "virtual_sitesn", "dummies1", "dummies2", "dummies3", "dummies4",
    "cmap",
})


def used_moleculetypes(lines: List[str]) -> Optional[Set[str]]:
    """The molecule names listed in ``[ molecules ]``.

    ``None`` when the file has no ``[ molecules ]`` section (an ``.itp`` on
    its own), which means "cannot tell" -- then everything counts.
    """
    names: Set[str] = set()
    seen = False
    in_section = False
    for line in lines:
        m = _SECTION_RE.match(line)
        if m:
            in_section = m.group(1).lower() == "molecules"
            seen = seen or in_section
            continue
        if in_section and _is_data_line(line):
            parts = line.split(";")[0].split()
            if parts:
                names.add(parts[0])
    return names if seen else None


def scan_sections(lines: List[str], used: Optional[Set[str]] = None) -> Dict[str, int]:
    """Count data (non-comment, non-blank) lines per ``[ section ]``.

    *lines* must already have ``#include`` resolved -- otherwise a section
    living in an ``.itp`` is missed, which is exactly where a coarse-grained
    force field keeps ``[ nonbond_params ]``.

    Sections inside a ``[ moleculetype ]`` are counted **only for the
    molecules the system actually contains** (:data:`PER_MOLECULE_SECTIONS`).
    Pass *used* to override the names read from ``[ molecules ]``. Without
    this, one library ``.itp`` fails every conversion that includes it: the
    Martini lipid library, for instance, carries imidazolium ions whose
    ``[ constraints ]`` have nothing to do with a bilayer of POPC.
    """
    if used is None:
        used = used_moleculetypes(lines)
    counts: Dict[str, int] = {}
    current: Optional[str] = None
    molecule: Optional[str] = None
    expect_name = False
    for line in lines:
        m = _SECTION_RE.match(line)
        if m:
            current = m.group(1).lower()
            if current == "moleculetype":
                expect_name, molecule = True, None
            counts.setdefault(current, 0)
            continue
        if current is None or not _is_data_line(line):
            continue
        if expect_name:
            parts = line.split(";")[0].split()
            if parts:
                molecule = parts[0]
                expect_name = False
        if (used is not None and current in PER_MOLECULE_SECTIONS
                and molecule is not None and molecule not in used):
            continue
        counts[current] += 1
    return counts


def scan_dihedral_functs(lines: List[str],
                         used: Optional[Set[str]] = None) -> List[int]:
    """Collect the ``funct`` column of every dihedral line, from the raw text.

    :class:`~.top_parser.TopParser` cannot be used for this. A dihedral line
    carrying inline parameters for an unsupported funct (e.g. ``2``, a
    harmonic improper) falls through to its "type reference only" branch,
    where column 5 is stored as a **torsion type index**. The funct is gone by
    the time parsing finishes, and what remains points at an unrelated type.

    ``[ dihedrals ]`` is ``ai aj ak al funct ...`` (atoms are integers), while
    ``[ dihedraltypes ]`` names atom types and comes in a 2-name and a 4-name
    flavour, so the funct is taken as the first integer token.

    ``[ dihedrals ]`` inside a ``[ moleculetype ]`` is read **only for the
    molecules the system actually contains**. A force field file carries
    molecules nobody used -- the Martini ion file defines imidazolium ions
    next to the sodium and chloride a bilayer needs -- and counting those
    made every bilayer conversion warn that it was dropping terms that were
    never in the box. ``[ dihedraltypes ]`` is not inside a molecule, so it
    is always read.
    """
    if used is None:
        used = used_moleculetypes(lines)
    functs = set()
    current: Optional[str] = None
    molecule: Optional[str] = None
    expect_name = False
    for line in lines:
        m = _SECTION_RE.match(line)
        if m:
            current = m.group(1).lower()
            if current == "moleculetype":
                expect_name, molecule = True, None
            continue
        if expect_name and _is_data_line(line):
            parts = line.split(";")[0].split()
            if parts:
                molecule, expect_name = parts[0], False
        if current not in ("dihedrals", "dihedraltypes"):
            continue
        if (current == "dihedrals" and molecule is not None
                and molecule not in used):
            continue
        if not _is_data_line(line):
            continue
        tokens = line.split(";")[0].split()
        if current == "dihedrals":
            if len(tokens) >= 5:
                try:
                    functs.add(int(tokens[4]))
                except ValueError:
                    pass
            continue
        for tok in tokens:
            try:
                functs.add(int(tok))
            except ValueError:
                continue
            break
    return sorted(functs)


def _functs(entries, idx) -> List[int]:
    out = set()
    for e in entries:
        if len(e) > idx and isinstance(e[idx], int):
            out.add(e[idx])
    return sorted(out)


#: How far from 0 or 180 degrees a phase may sit and still count as that
#: value. A force field writes the phase through radians, so the number in
#: the file is rounded: ff14SB through parmed gives **180.0000771** degrees
#: (3.14159 rad), which an exact test refuses. The phases this check exists
#: to catch are tens of degrees away, so a hundredth of a degree separates
#: the two cases with room to spare.
PHASE_TOLERANCE_DEG = 0.01


def phase_is_writable(phase) -> bool:
    """True when *phase* is 0 or 180 degrees, which is what the engine takes.

    The value arrives as a number from an inline ``[ dihedrals ]`` line and as
    **text** from a ``[ dihedraltypes ]`` lookup. Taking ``%`` of text is
    string formatting, not a remainder, so the test floats it first. Anything
    that will not float is left to the writer rather than refused here.
    """
    try:
        value = float(phase) % 360.0
    except (TypeError, ValueError):
        return True
    return min(abs(value), abs(value - 180.0),
               abs(value - 360.0)) <= PHASE_TOLERANCE_DEG


def check_top(raw, sections: Dict[str, int],
              dihedral_functs: Optional[List[int]] = None
              ) -> Tuple[List[str], List[str]]:
    """Return ``(fatal, warnings)``, each a list of human-readable messages.

    Two empty lists mean gro2udf can represent everything it found.
    """
    fatal: List[str] = []
    warnings: List[str] = []

    if getattr(raw, "comb_rule", 2) == 1:
        fatal.append(
            "[ defaults ] declares comb-rule 1, so the last two [ atomtypes ] "
            "columns are c6 and c12 -- gro2udf reads them as sigma and epsilon. "
            "Every LJ parameter would be wrong; for Martini they are 0.0, "
            "which also makes Interaction_Site_Type[].Range come out as 0")

    for name, count in sorted(sections.items()):
        if not count:
            continue
        plural = "entry" if count == 1 else "entries"
        if name in FATAL_SECTIONS:
            fatal.append("[ {} ] has {} {} but is not read; this drops {}"
                         .format(name, count, plural, FATAL_SECTIONS[name]))
        elif name in DROPPED_SECTIONS:
            warnings.append("[ {} ] has {} {} but is not read; this drops {}"
                            .format(name, count, plural, DROPPED_SECTIONS[name]))

    # An explicit exclusion is written as a bond with no force constant, and
    # that only removes the one pair while nrexcl is 1. With a larger nrexcl
    # the same bond would also exclude the neighbours it newly puts within
    # reach, which the topology never asked for.
    # Only the pairs nrexcl does NOT already cover are a real request: water
    # out of an Amber force field lists its three pairs although its bonds
    # cover them, and refusing those refused every TIP3P system.
    excl = getattr(raw, "exclusions_needed", None)
    if excl is None:
        excl = getattr(raw, "exclusions", {}) or {}
    nrexcl = getattr(raw, "nrexcl", {}) or {}
    bad = sorted({m for m, pairs in excl.items()
                  if pairs and nrexcl.get(m, 1) >= 2})
    if bad:
        fatal.append(
            "[ exclusions ] in {} with nrexcl >= 2. An exclusion is carried "
            "as a bond with no force constant, which removes exactly the one "
            "pair only while nrexcl is 1; with a larger nrexcl that bond "
            "would exclude further neighbours as well".format(", ".join(bad)))

    # A proper dihedral is written as the engine's Amber form, whose phase
    # must be 0 or 180 degrees. Any other phase carries a sin(n phi) part,
    # which none of the engine's torsion forms can hold -- not the cosine
    # polynomial (even in phi), and not the table either, which is indexed
    # by cos(phi) and so cannot tell phi from -phi.
    odd = [p for p in getattr(raw, "torsion_phases", []) or []
           if not phase_is_writable(p)]
    if odd:
        fatal.append(
            "{} proper dihedral(s) have a phase that is neither 0 nor 180 "
            "degrees ({} distinct, e.g. {:.1f}). The engine's Amber torsion "
            "takes only those two, and its other forms are even in phi, so "
            "such a term cannot be written at all".format(
                len(odd), len({round(p, 3) for p in odd}), odd[0]))

    # A bond or angle is written unconditionally as Harmonic / Theta, so an
    # unsupported funct becomes a wrong number rather than a missing term.
    for kind, supported, target, found in (
            ("bond", SUPPORTED_BOND_FUNCTS, "COGNAC Harmonic",
             _functs(list(raw.bondtypes) + list(raw.bond_types_from_mol), 2)),
            ("angle", SUPPORTED_ANGLE_FUNCTS, "COGNAC Theta",
             _functs(list(raw.angletypes) + list(raw.angle_types_from_mol), 3))):
        for funct in found:
            if funct in supported:
                continue
            note = _FUNCT_NOTE.get((kind, funct))
            fatal.append(
                "{} funct {} is written as {} regardless (only funct {} maps "
                "to it){}".format(
                    kind, funct, target,
                    "/".join(str(f) for f in sorted(supported)),
                    ": " + note if note else ""))

    # A dihedral with an unsupported funct never reaches the exporter, so no
    # wrong number is written -- but the line is misread as a reference to a
    # torsion type that does not exist, so the term is lost.
    found_dihedral = (dihedral_functs if dihedral_functs is not None
                      else _functs(list(raw.torsiontypes)
                                   + list(raw.torsion_types_from_mol), 4))
    for funct in found_dihedral:
        if funct not in SUPPORTED_TORSION_FUNCTS:
            warnings.append(
                "dihedral funct {} is dropped (only funct {} are written); "
                "the line is read as a reference to torsion type {}, which "
                "does not exist".format(
                    funct, "/".join(str(f)
                                    for f in sorted(SUPPORTED_TORSION_FUNCTS)),
                    funct))

    return fatal, warnings


def raise_if_unsupported(raw, sections: Dict[str, int],
                         allow_unsupported: bool = False,
                         dihedral_functs: Optional[List[int]] = None
                         ) -> Tuple[List[str], List[str]]:
    """Run :func:`check_top`; raise on fatal findings unless overridden.

    Warnings are always logged. Returns ``(fatal, warnings)``.
    """
    fatal, warnings = check_top(raw, sections, dihedral_functs)

    if warnings:
        logger.warning("gro2udf will drop the following; the UDF will not "
                       "reproduce them:\n%s",
                       "\n".join("  - " + w for w in warnings))

    if fatal:
        body = "\n".join("  - " + f for f in fatal)
        if allow_unsupported:
            logger.warning(
                "converting anyway (--allow-unsupported); these terms will be "
                "WRONG in the UDF, not merely missing:\n%s", body)
        else:
            raise UnsupportedTopFeatureError(
                "this .top uses features gro2udf writes incorrectly:\n"
                + body
                + "\n\nWithout this check the conversion would have reported "
                  "success and produced a UDF that is quietly wrong. Pass "
                  "--allow-unsupported to convert anyway.")

    return fatal, warnings
