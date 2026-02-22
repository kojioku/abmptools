# -*- coding: utf-8 -*-
"""
cli.py
------
Command-line entry point for gro2udf.

Maintains full backward compatibility with the original::

    python gro2udf.py <input.udf> <input.gro>
    python gro2udf.py <input.udf> <input.gro> <input.xvg>   (xvg ignored)

Can also be called via the installed package::

    python -m abmptools.gro2udf <input.udf> <input.gro>

New mode: convert from GROMACS TOP + GRO to COGNAC UDF::

    python -m abmptools.gro2udf --from-top system.top output.gro \\
        [--template template.udf] [--out result.udf]

    --template  : existing COGNAC UDF to use as schema (default: <topfile_stem>.udf)
    --out       : output UDF path (default: <grofile_stem>_fromtop.udf)
"""
from __future__ import annotations

import os
import sys


def _run_from_top(argv: list) -> None:
    """Handle --from-top mode."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="gro2udf --from-top",
        description="Convert GROMACS TOP+GRO to COGNAC UDF",
        add_help=True,
    )
    parser.add_argument("top_path", help="GROMACS .top file")
    parser.add_argument("gro_path", help="GROMACS .gro file")
    parser.add_argument("--mdp", dest="mdp_path", default=None,
                        help="GROMACS .mdp file (accepted but currently not processed)")
    parser.add_argument("--template", dest="template_path", default=None,
                        help="Existing COGNAC UDF file (schema template)")
    parser.add_argument("--out", dest="out_path", default=None,
                        help="Output UDF file path")

    # Strip the --from-top flag from argv before parsing
    filtered = [a for a in argv[1:] if a != "--from-top"]
    args = parser.parse_args(filtered)

    top_path = args.top_path
    gro_path = args.gro_path

    # Resolve template path
    template_path = args.template_path
    if template_path is None:
        top_stem = os.path.splitext(top_path)[0]
        candidate = top_stem + ".udf"
        if os.path.isfile(candidate):
            template_path = candidate
        else:
            print("ERROR: --template not specified and '{}' not found.".format(candidate))
            sys.exit(1)

    # Resolve output path
    out_path = args.out_path
    if out_path is None:
        gro_stem = os.path.splitext(os.path.basename(gro_path))[0]
        out_path = gro_stem + "_fromtop.udf"

    from .top_exporter import TopExporter
    TopExporter().export(top_path, gro_path, template_path, out_path)
    print("Written: {}".format(out_path))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if "--from-top" in argv:
        _run_from_top(argv)
        return

    if len(argv) < 3:
        print("Usage: {} udffile grofile [xvg]".format(
            os.path.basename(argv[0])
        ))
        print("       {} --from-top topfile grofile [--template t.udf] [--out out.udf]".format(
            os.path.basename(argv[0])
        ))
        raise RuntimeError("Illegal arguments.")

    udf_path = argv[1]
    gro_path = argv[2]
    # argv[3] (xvg) is accepted for CLI compatibility but not processed here

    from .exporter import Exporter
    return Exporter().export(udf_path, gro_path)


if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError:
        print("ERROR: gro2udf failed.")
        sys.exit(1)
