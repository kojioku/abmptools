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

#: cognac10.1 (OCTA8.4 / OCTA8.4) compatible bundled template.
#: Used automatically when ``--cognac-version 101`` / ``102`` is given and no
#: ``--template`` override is supplied.
_BUILTIN_TEMPLATE_COGNAC101 = os.path.join(
    os.path.dirname(__file__), "default_template_cognac101.udf"
)


def _from_top_parser():
    """The --from-top option parser.

    Split out so ``--help`` can show the real option list. It used to print
    a two-line usage and tell the reader to run ``--from-top --help`` for
    the rest, which meant ``--trajectory`` -- the reason most people come
    here -- was not in the help they actually read (2026-09-16).
    """
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
    parser.add_argument("--nh-dof", dest="nh_dof", default="3N-3",
                        choices=["3N", "3N-3"],
                        help="Degrees of freedom behind the Nose-Hoover Q, "
                             "and with it whether centre-of-mass motion is "
                             "removed. 3N-3 (default) asks for "
                             "comm-mode = Linear, which is GROMACS' own "
                             "default and stops kinetic energy accumulating "
                             "in the centre of mass. 3N leaves "
                             "comm-mode = None, the usual UDF convention, "
                             "and reproduces its output exactly. The "
                             "difference in Q is 0.03%% at 3050 atoms and "
                             "0.6%% at 80.")
    parser.add_argument("--ff", dest="force_field", default="gaff",
                        help="Force field to declare in Unit_Parameter.Comment "
                             "as FF=n, where downstream converters read it. A .top "
                             "does not say which force field it came from, so "
                             "it has to be named here. Accepts gaff (default), "
                             "gaff2, amber, amber20, dreiding, uff, oplsaa, "
                             "loplsaa, loplsaa2023, pcff, a bare number, or "
                             "'' to leave it unset. An existing value in the "
                             "template is kept.")
    parser.add_argument("--template", dest="template_path", default=None,
                        help="Existing COGNAC UDF file (schema template). "
                             "Defaults to the built-in template; a .udf next to the "
                             ".top is never picked up on its own.")
    parser.add_argument("--out", dest="out_path", default=None,
                        help="Output UDF file path")
    parser.add_argument("--cognac-version", dest="cognac_version", default=None,
                        help="Override the cognac<N>.udf include in the template "
                             "at runtime. Use e.g. `--cognac-version 110` when "
                             "your OCTA install (e.g. OCTA84 / OCTA viewer 9.1) does "
                             "not ship cognac112.udf. The bundled default "
                             "template requests cognac112; lower versions like "
                             "110/111 are compatible because gro2udf only "
                             "writes fields that exist in all of them.")
    parser.add_argument("--topology-only", dest="topology_only", action="store_true",
                        help="Write the topology (Set_of_Molecules / "
                             "Molecular_Attributes / Interactions) only, "
                             "without any Structure record (skeleton UDF). "
                             "Use this when OCTA viewer (GOURMET) will load the "
                             "trajectory and energy separately. Pair with "
                             "`--initial-gro <path>` to also embed a single "
                             "initial frame.")
    parser.add_argument("--initial-gro", dest="initial_gro_path", default=None,
                        help="Path to a .gro file whose first frame is "
                             "written as the single initial Structure record. "
                             "Only meaningful with --topology-only; without "
                             "this option the topology-only UDF has zero "
                             "Structure records.")
    parser.add_argument("--trajectory", dest="trajectory_path", default=None,
                        help="Path to a multi-frame trajectory (a raw .xtc, "
                             "or a multi-frame .gro). Every frame is written "
                             "as a Structure record in the output UDF, so "
                             "the trajectory plays back inside the UDF "
                             "without an external attach step. **It is put "
                             "through `gmx trjconv -pbc nojump` first**, so a "
                             "raw trajectory gives whole molecules without "
                             "being asked -- which is the one thing here that "
                             "needs gmx. See --skip-nojump.")
    parser.add_argument("--energy", dest="energy_path", default=None,
                        help="Path to an .xvg file (e.g. output of `gmx "
                             "energy`). Each frame's row is written to the "
                             "corresponding Structure record's "
                             "`Statistics_Data.Energy.Instantaneous` "
                             "(Bond / Angle / Torsion / Nonbonding / "
                             "Electrostatic terms mapped from xvg legend "
                             "names). Times in the .xvg are matched to "
                             "frame times via nearest-neighbour interpolation.")
    parser.add_argument("--max-frames", dest="max_frames", type=int, default=None,
                        help="Keep at most this many trajectory frames **in "
                             "total**, thinning on read. This is a count, not "
                             "a stride: --max-frames 100 on a 10000-frame xtc "
                             "keeps every 100th. The count is what the UDF "
                             "size tracks, so it is the knob to reach for. "
                             "Does not need gmx. Cannot be combined with "
                             "--frame-step.")
    parser.add_argument("--frame-step", dest="frame_step", type=int, default=None,
                        help="Keep every Nth trajectory frame (a stride). Use "
                             "--max-frames when you care about the total "
                             "instead. Does not need gmx.")
    parser.add_argument("--edr", dest="edr_path", default=None,
                        help="GROMACS .edr. The .xvg is produced on the spot "
                             "(via abmptools.trajectory) and embedded, so you "
                             "do not have to run `gmx energy` yourself. "
                             "**Needs gmx.** Without this flag no energy is "
                             "read at all -- pass --energy instead if you "
                             "already have the .xvg.")
    # --trajectory を渡したら、 既定で `gmx trjconv -pbc nojump` を通す。
    # 生の .xtc をそのまま渡した人が、 何も指定せずに正しい UDF を得られる
    # ようにするため -- 割れた分子は OCTA viewer で開いて初めて分かるので、
    # 既定が OFF だと「気付ける人しか気付けない」形になっていた。
    #
    # **nojump は冪等** (既に nojump 済みの軌跡に掛け直しても座標は変わらない。
    # PVA 30 分子で実測して max |dr| = 0.0000 A / 移動した原子 0 個) なので、
    # gen_for_udf の出力を渡す流れも既定のまま通る。
    #
    # フラグは**この段を飛ばす**とだけ言う。 一度 --already-nojump にして
    # みたが、 あれは「入力はもう nojump 済みだ」と**入力の性質を主張**する。
    # 済んでいなくても意図して飛ばすこと (gmx が無い、 別の後処理で通す、
    # 割れたままの座標を見たい) はあるので、 そのときに --already-nojump と
    # 書かせるのは**嘘を書かせる**ことになる。 動作だけを述べればどちらでも
    # 正しい。
    parser.add_argument("--skip-nojump", dest="skip_nojump",
                        action="store_true",
                        help="Do not run `gmx trjconv -pbc nojump` on "
                             "--trajectory. Without it a trajectory is put "
                             "through nojump first -- molecules that are "
                             "split across the periodic boundary show up "
                             "broken in the OCTA viewer and downstream, and a "
                             "raw .xtc is split. Skip it when the trajectory "
                             "has already been through nojump (running it "
                             "twice changes nothing, so this is about not "
                             "needing gmx), or when there is no gmx here: "
                             "**--trajectory needs gmx unless you pass this.**")
    parser.add_argument("--tpr", dest="tpr_path", default=None,
                        help="Reference for the -pbc nojump step. Optional: "
                             "without it the .gro argument is used, which "
                             "works because -pbc nojump reads no bonded "
                             "information. Given one that this gmx cannot "
                             "read (a newer tpx version), it falls back to "
                             "the .gro and says so.")
    parser.add_argument("--gmx", dest="gmx", default="gmx",
                        help="gmx executable to use for the -pbc nojump step "
                             "and --edr (default: resolved on PATH). Point it "
                             "at your MD environment's gmx when that one is "
                             "not on PATH.")
    parser.add_argument("--allow-unsupported", dest="allow_unsupported",
                        action="store_true",
                        help="Convert even when the .top contains terms "
                             "gro2udf writes incorrectly (comb-rule 1, "
                             "[ nonbond_params ], [ constraints ], non-"
                             "harmonic bond/angle functs). Without this the "
                             "conversion stops instead of reporting success "
                             "and writing a UDF that is quietly wrong. Terms "
                             "that are merely dropped always warn and never "
                             "block.")

    return parser


def _run_from_top(argv: list) -> None:
    """Handle --from-top mode."""
    parser = _from_top_parser()

    # Strip the --from-top flag from argv before parsing
    filtered = [a for a in argv[1:] if a != "--from-top"]
    args = parser.parse_args(filtered)

    top_path = args.top_path
    gro_path = args.gro_path

    # --- Resolve template path ---
    template_path = args.template_path
    if template_path is None:
        # <top_stem>.udf used to be picked up here without being asked for.
        # That is almost always the *pre-MD* input sitting next to the .top,
        # and a template supplies the static structure and box -- so the
        # conversion silently inherited the box the system had before it ran.
        # Say it is there, and let the user ask for it.
        top_stem = os.path.splitext(top_path)[0]
        candidate = top_stem + ".udf"
        if os.path.isfile(candidate):
            print("Note: {} exists but is NOT used. Templates are only used "
                  "when asked for: pass --template {} if that is what you "
                  "want.".format(candidate, candidate))

        # When the user explicitly asked for a cognac10.x schema
        # (OCTA8.4 / OCTA8.4), pick the cognac101-compatible
        # bundled template so its data section parses on that install.
        # NOTE: enumerate cognac10.x explicitly -- `str.startswith("10")`
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

    # --- 間引きの指定は 2 通りあるが、意味が違うので同時には受けない ---
    if args.max_frames is not None and args.frame_step is not None:
        raise RuntimeError(
            "--max-frames and --frame-step do different things and cannot be "
            "combined: --max-frames is a total count, --frame-step is a stride."
        )
    if args.max_frames is not None and args.max_frames < 1:
        raise RuntimeError("--max-frames must be >= 1")
    if args.frame_step is not None and args.frame_step < 1:
        raise RuntimeError("--frame-step must be >= 1")

    # --- gmx が要るのはここだけ。ファイル変換そのものは gmx を呼ばない ---
    trajectory_path = args.trajectory_path
    energy_path = args.energy_path

    # --trajectory があるときだけ走る。 topology だけの変換に gmx を
    # 要求しない。 --skip-nojump は --trajectory に対する指定なので、
    # 軌跡が無いのに書いてあるのは書き間違い -- 黙って通すと
    # **--trajectory を書き忘れた人に topology だけの UDF が出る**。
    if args.skip_nojump and not trajectory_path:
        raise RuntimeError(
            "--skip-nojump applies to --trajectory, which was not given")

    if trajectory_path and not args.skip_nojump:
        from ..trajectory.postprocess import nojump_with_fallback
        # --tpr が無ければ、 位置引数の .gro をそのまま reference にする。
        # -pbc nojump は結合情報を使わないので .gro で成立する。
        reference = args.tpr_path or gro_path
        # 古い gmx は新しい .tpr を読めない (tpx の版違い)。 そのときは
        # .gro へ退避する。 ここでは位置引数の .gro が必ずあるので、
        # 利用者が別途用意する必要はない。
        fallback = gro_path if args.tpr_path else None
        try:
            trajectory_path, used = nojump_with_fallback(
                trajectory=trajectory_path, reference=reference,
                fallback=fallback, gmx=args.gmx,
            )
        except FileNotFoundError as exc:
            # gmx が無い機。 既定で走るようになった段なので、**何を止めれば
            # 変換が通るのか**をここで言う。 言わないと「gro2udf は gmx を
            # 呼ばない変換器」という以前の理解のまま詰まる。
            if "gmx" not in str(exc):
                raise
            raise RuntimeError(
                "{}\n"
                "-pbc nojump runs whenever --trajectory is given, so gmx is "
                "needed here. Point --gmx at it, or pass --skip-nojump to "
                "leave the trajectory alone -- which is what you want when it "
                "has already been through -pbc nojump (for example the .gro "
                "that `abmptools.trajectory gen_for_udf` writes)."
                .format(exc)) from exc
        trajectory_path = str(trajectory_path)
        if used is not None:
            print("Prepared (-pbc nojump): {}\n"
                  "  reference: {} (this gmx could not read {})"
                  .format(trajectory_path, used, args.tpr_path))
        else:
            print("Prepared (-pbc nojump): {}\n  reference: {}"
                  .format(trajectory_path, reference))

    if args.edr_path:
        if energy_path:
            raise RuntimeError(
                "--edr and --energy both give the energy; pass only one "
                "(--edr makes the .xvg here, --energy takes one you have)")
        from ..trajectory.postprocess import gmx_energy
        energy_path = str(gmx_energy(edr=args.edr_path,
                                     output=os.path.splitext(args.edr_path)[0]
                                     + "_energy.xvg",
                                     terms=range(1, 51),
                                     gmx=args.gmx))
        print("Energy dumped: {}".format(energy_path))

    from .top_exporter import TopExporter
    TopExporter().export(top_path, gro_path, template_path, out_path,
                         mdp_path=args.mdp_path,
                         cognac_version=args.cognac_version,
                         topology_only=args.topology_only,
                         initial_gro_path=args.initial_gro_path,
                         trajectory_path=trajectory_path,
                         energy_path=energy_path,
                         max_frames=args.max_frames,
                         frame_step=args.frame_step,
                         allow_unsupported=args.allow_unsupported,
                         force_field=args.force_field,
                         nh_dof=args.nh_dof)
    print("Written: {}".format(out_path))
    if args.topology_only:
        if args.initial_gro_path:
            print(f"  (topology + 1 initial frame from {args.initial_gro_path!s}. "
                  f"Load further trajectory / energy in OCTA viewer (GOURMET).)")
        else:
            print("  (topology-only -- no Structure record. Load trajectory "
                  "/ energy directly in OCTA viewer (GOURMET).)")
    elif args.trajectory_path:
        print("  (embedded {} frames{})".format(
            "trajectory",
            " + energy" if args.energy_path else ""))


def _usage(argv) -> None:
    """Print how to call this. Used by --help and by the argument error."""
    prog = os.path.basename(argv[0])
    print("Usage: {} udffile grofile [xvg]".format(prog))
    print("       {} --from-top topfile grofile "
          "[--template t.udf] [--out out.udf]".format(prog))
    print("       {} --from-top topfile grofile "
          "[--trajectory md.xtc] [--energy e.xvg]".format(prog))
    print()
    print("Legacy mode takes an existing UDF as the schema and replaces its")
    print("coordinates from the .gro. --from-top builds the UDF from a GROMACS")
    print("topology instead, and can embed a whole trajectory:")
    print()
    print("  gro2udf --from-top system.top system.gro \\")
    print("          --trajectory md.xtc --energy energy.xvg --out out.udf")
    print()
    print("--trajectory runs `gmx trjconv -pbc nojump` on it first, so a raw")
    print(".xtc gives whole molecules without being asked. That is the only")
    print("place this needs gmx; --skip-nojump leaves the trajectory alone.")
    print()
    print(_from_top_parser().format_help())


def main(argv=None):
    """gro2udfのコマンドラインエントリポイント。

    通常モードとTOP変換モード (--from-top) をサポートする。
    """
    if argv is None:
        argv = sys.argv

    if "--from-top" in argv:
        _run_from_top(argv)
        return

    # Asking for the usage is not a failure. This printed it correctly and
    # then raised, so `gro2udf --help` came back as 1 and a caller checking
    # the exit code read it as one -- the same thing that was fixed in
    # moldeck's install_minimum.bat and in moldeck.cg.martini_top.
    if any(a in ("-h", "--help") for a in argv[1:]):
        _usage(argv)
        return

    if len(argv) < 3:
        _usage(argv)
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
