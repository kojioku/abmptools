# -*- coding: utf-8 -*-
"""
abmptools.membrane.lipid_info
-----------------------------
Discovery / reporting helpers for the lipid catalogue.

Two layers:

- :mod:`abmptools.membrane.bilayer` exports :data:`DEFAULT_LIPID_APL`
  (the curated APL table used by :func:`_resolve_apl`),
  :func:`list_known_lipids` (filterable view of the table), and
  :func:`query_packmol_memgen_lipids` (runtime parse of the binary's
  full 259-entry catalogue).

- This module provides a CLI on top:
  ``python -m abmptools.membrane.lipid_info [--known | --available |
                                             --head PC|PE|PG|PS|PA|SM|sterol]``

Examples
--------

List every curated lipid (~60 entries) with APL:

    python -m abmptools.membrane.lipid_info --known

List only sphingomyelins:

    python -m abmptools.membrane.lipid_info --known --head SM

Query packmol-memgen for the full 259-lipid catalogue (incl. those
abmptools has no curated APL for):

    python -m abmptools.membrane.lipid_info --available

Show APL for a single lipid (returns 65.0 fallback if unknown):

    python -m abmptools.membrane.lipid_info --apl POPC
"""
from __future__ import annotations

import argparse
import logging
import sys

from .bilayer import (
    DEFAULT_LIPID_APL,
    _resolve_apl,
    list_known_lipids,
    query_packmol_memgen_lipids,
)
from .models import LipidSpec

logger = logging.getLogger(__name__)


def _format_known_table(rows: list[tuple[str, float, str]]) -> str:
    """Format list_known_lipids output as a fixed-width table."""
    if not rows:
        return "(no lipids match the filter)"
    lines = []
    lines.append(f"  {'resname':<8s} {'APL (Å²)':>10s}  head")
    lines.append("  " + "─" * 30)
    cur_head = None
    for resname, apl, head in rows:
        if cur_head is not None and head != cur_head:
            lines.append("")  # blank line between head groups
        cur_head = head
        lines.append(f"  {resname:<8s} {apl:>10.1f}  {head}")
    return "\n".join(lines)


def _format_available_table(rows: list[tuple[str, int, str]]) -> str:
    """Format query_packmol_memgen_lipids output as a fixed-width table."""
    if not rows:
        return "(packmol-memgen returned no lipids — env / binary issue?)"
    lines = []
    lines.append(f"  {'resname':<8s} {'charge':>6s}  description")
    lines.append("  " + "─" * 80)
    for resname, charge, desc in rows:
        # truncate long descriptions
        if len(desc) > 60:
            desc = desc[:57] + "..."
        lines.append(f"  {resname:<8s} {charge:>+6d}  {desc}")
    return "\n".join(lines)


def _main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Discover and inspect lipid species supported by "
            "abmptools.membrane (curated APL table) and packmol-memgen "
            "(full 259-entry catalogue)."
        ),
    )
    parser.add_argument(
        "--known", action="store_true",
        help="Print the curated DEFAULT_LIPID_APL table "
             "(~60 entries with APL values).",
    )
    parser.add_argument(
        "--available", action="store_true",
        help="Query packmol-memgen for its full 259-lipid catalogue "
             "(incl. lipids that have no curated APL — those default to "
             "65 Å² fallback).",
    )
    parser.add_argument(
        "--head", choices=["PC", "PE", "PG", "PS", "PA", "SM", "sterol"],
        help="Filter --known by head-group classification.",
    )
    parser.add_argument(
        "--apl", metavar="RESNAME",
        help="Print the resolved APL value for a single residue name "
             "(table lookup; falls back to 65.0 for unknown).",
    )
    parser.add_argument(
        "--packmol-memgen-path", default="packmol-memgen",
        help="Override the packmol-memgen binary (default: PATH-resolved).",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")

    if args.apl:
        # Single-lipid lookup
        spec = LipidSpec(resname=args.apl, n_per_leaflet=1)
        apl = _resolve_apl(spec)
        in_table = args.apl in DEFAULT_LIPID_APL
        print(
            f"  resname:        {args.apl}\n"
            f"  APL (Å²):       {apl}\n"
            f"  in curated table: {in_table}"
        )
        return

    if not (args.known or args.available):
        # Default: print --known summary
        args.known = True

    if args.known:
        rows = list_known_lipids(head_group=args.head)
        head_filter = f" ({args.head} only)" if args.head else ""
        print(f"=== Curated DEFAULT_LIPID_APL{head_filter} "
              f"— {len(rows)} entries ===")
        print(_format_known_table(rows))
        if args.available:
            print()  # separator before --available section

    if args.available:
        try:
            rows = query_packmol_memgen_lipids(
                packmol_memgen_path=args.packmol_memgen_path
            )
        except FileNotFoundError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            print(
                f"  (set $AMBERHOME / PATH so packmol-memgen is reachable, "
                f"or pass --packmol-memgen-path /abs/path/to/packmol-memgen)",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"=== packmol-memgen --available_lipids "
              f"— {len(rows)} entries ===")
        print(_format_available_table(rows))


if __name__ == "__main__":
    _main()
