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
        [--template template.udf] [--mdp system.mdp] [--out result.udf]

    --template  : existing COGNAC UDF to use as schema template.
                  Default: <topfile_stem>.udf in the same directory,
                  then the built-in default_template.udf bundled in this package.
    --mdp       : GROMACS .mdp file; when given, Nose-Hoover Q and
                  Ewald R_cutoff are derived from its values.
    --out       : output UDF path (default: <grofile_stem>_fromtop.udf)
"""
from __future__ import annotations

import os
import sys

#: Built-in fallback template bundled with this package.
_BUILTIN_TEMPLATE = os.path.join(os.path.dirname(__file__), "default_template.udf")

#: cognac10.1 (OCTA8.4 / J-OCTA-9.1-Student) compatible bundled template.
#: Used automatically when ``--cognac-version 101`` / ``102`` is given and no
#: ``--template`` override is supplied.
_BUILTIN_TEMPLATE_COGNAC101 = os.path.join(
    os.path.dirname(__file__), "default_template_cognac101.udf"
)


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
                        help="GROMACS .mdp file (ref_t, tau_t, rcoulomb are read)")
    parser.add_argument("--template", dest="template_path", default=None,
                        help="Existing COGNAC UDF file (schema template). "
                             "Defaults to <top_stem>.udf → built-in template.")
    parser.add_argument("--out", dest="out_path", default=None,
                        help="Output UDF file path")
    parser.add_argument("--cognac-version", dest="cognac_version", default=None,
                        help="Override the cognac<N>.udf include in the template "
                             "at runtime. Use e.g. `--cognac-version 110` when "
                             "your OCTA install (e.g. OCTA84 / J-OCTA 9.1) does "
                             "not ship cognac112.udf. The bundled default "
                             "template requests cognac112; lower versions like "
                             "110/111 are compatible because gro2udf only "
                             "writes fields that exist in all of them.")
    parser.add_argument("--topology-only", dest="topology_only", action="store_true",
                        help="Write the topology (Set_of_Molecules / "
                             "Molecular_Attributes / Interactions) but no "
                             "Structure record. Use this when you want to "
                             "load coordinates from a separate trajectory "
                             "(e.g. multi-frame .gro from `gmx trjconv -pbc "
                             "nojump`) and energy from a .xvg directly in "
                             "J-OCTA Viewer, instead of baking a single frame "
                             "into the UDF. The .gro path is still required "
                             "(used only to read the box vector and to confirm "
                             "atom count matches the topology).")

    # Strip the --from-top flag from argv before parsing
    filtered = [a for a in argv[1:] if a != "--from-top"]
    args = parser.parse_args(filtered)

    top_path = args.top_path
    gro_path = args.gro_path

    # --- Resolve template path ---
    template_path = args.template_path
    if template_path is None:
        # 1st priority: <top_stem>.udf in the same directory
        top_stem = os.path.splitext(top_path)[0]
        candidate = top_stem + ".udf"
        if os.path.isfile(candidate):
            template_path = candidate
            print("Template: {} (auto-detected)".format(template_path))
        else:
            # When the user explicitly asked for a cognac10.x schema
            # (OCTA8.4 / J-OCTA-9.1-Student), pick the cognac101-compatible
            # bundled template so its data section parses on that install.
            # NOTE: enumerate cognac10.x explicitly — `str.startswith("10")`
            # would erroneously match `"110"`/`"112"` (those are cognac 11.x,
            # not cognac 10.x).
            cv = args.cognac_version
            cognac10x = {"100", "101", "102"}
            if cv is not None and str(cv) in cognac10x:
                template_path = _BUILTIN_TEMPLATE_COGNAC101
                print("Template: {} (built-in cognac10.x default)".format(
                    template_path))
            else:
                # Default: cognac11.2 (OCTA85)
                template_path = _BUILTIN_TEMPLATE
                print("Template: {} (built-in default)".format(template_path))

    # --- Resolve output path ---
    out_path = args.out_path
    if out_path is None:
        gro_stem = os.path.splitext(os.path.basename(gro_path))[0]
        out_path = gro_stem + "_fromtop.udf"

    from .top_exporter import TopExporter
    TopExporter().export(top_path, gro_path, template_path, out_path,
                         mdp_path=args.mdp_path,
                         cognac_version=args.cognac_version,
                         topology_only=args.topology_only)
    print("Written: {}".format(out_path))
    if args.topology_only:
        print("  (topology-only — no Structure record. Load coordinates "
              "from a .gro/.xtc and energy from a .xvg directly in J-OCTA "
              "Viewer.)")


def main(argv=None):
    """gro2udfのコマンドラインエントリポイント。

    通常モードとTOP変換モード (--from-top) をサポートする。
    """
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
    import traceback
    try:
        main(sys.argv)
    except SystemExit:
        raise
    except Exception as exc:
        # Show the full diagnostic message + traceback so users see the
        # section context attached by top_exporter.UDFExportError instead
        # of a bare "gro2udf failed".
        print(f"ERROR: gro2udf failed: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        print("", file=sys.stderr)
        print("--- traceback ---", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
