"""CLI for ``python -m abmptools.trajectory``.

Subcommands map 1-to-1 to :mod:`abmptools.trajectory.postprocess` functions.
Cross-platform (Windows / Linux / macOS).

Examples
--------
::

    python -m abmptools.trajectory thin_nojump \\
        --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10

    python -m abmptools.trajectory wrap_pbc \\
        --traj 05_npt_final.xtc --tpr 05_npt_final.tpr \\
        --out 05_npt_final_pbc.xtc --ur compact

    python -m abmptools.trajectory nojump \\
        --traj 05_npt_final.xtc --tpr 05_npt_final.tpr \\
        --out 05_npt_final_nojump.gro --group System
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .postprocess import GmxError, nojump, thin, thin_and_nojump, wrap_pbc


def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--traj", required=True, help="入力 trajectory (.xtc 等)")
    p.add_argument("--tpr", required=True, help="reference structure (.tpr 推奨)")
    p.add_argument("--out", default=None,
                   help="出力 path (default: <stem>_<tag>.<ext>)")
    p.add_argument("--group", default="System",
                   help="trjconv の group 名 (default: System)")
    p.add_argument("--ndx", default=None, help="index file (default: なし)")
    p.add_argument("--gmx", default="gmx",
                   help="gmx 実行 path (default: PATH 解決)")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m abmptools.trajectory",
        description="Cross-platform GROMACS trajectory post-processor "
                    "(thin + nojump unwrap + PBC wrap)。 sample/formulation の "
                    "aggregation run + amorphous の wrap_pbc / gen_for_udf 共通基盤。",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # thin_nojump (主用途)
    p_tn = sub.add_parser(
        "thin_nojump",
        help="間引き + nojump unwrap (aggregation 用の基本セット)",
    )
    _add_common_args(p_tn)
    p_tn.add_argument("--skip", type=int, default=10,
                      help="frame 間引き factor (default: 10)")

    # nojump 単発
    p_nj = sub.add_parser(
        "nojump",
        help="-pbc nojump のみ (frame 数そのまま、 OCTA / J-OCTA Viewer 用)",
    )
    _add_common_args(p_nj)

    # thin 単発
    p_th = sub.add_parser(
        "thin",
        help="-skip N のみ (PBC 処理なし)",
    )
    _add_common_args(p_th)
    p_th.add_argument("--skip", type=int, default=10,
                      help="frame 間引き factor (default: 10)")

    # wrap_pbc 単発
    p_wr = sub.add_parser(
        "wrap_pbc",
        help="-pbc mol -ur <ur> (+ -center) で VMD 用 wrap",
    )
    _add_common_args(p_wr)
    p_wr.add_argument("--ur", default="compact",
                      choices=["compact", "rect", "tric"],
                      help="-ur 値 (default: compact)")
    p_wr.add_argument("--center", default=None,
                      help="-center を有効化し、 指定 group を box 中央に置く")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        common = dict(
            trajectory=args.traj,
            tpr=args.tpr,
            output=args.out,
            group=args.group,
            ndx=args.ndx,
            gmx=args.gmx,
        )
        if args.cmd == "thin_nojump":
            out = thin_and_nojump(skip=args.skip, **common)
        elif args.cmd == "nojump":
            out = nojump(**common)
        elif args.cmd == "thin":
            out = thin(skip=args.skip, **common)
        elif args.cmd == "wrap_pbc":
            out = wrap_pbc(ur=args.ur, center=args.center, **common)
        else:  # pragma: no cover - argparse enforces choices
            parser.error(f"unknown command: {args.cmd}")
            return 2
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    except GmxError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return e.returncode if e.returncode != 0 else 1

    print(f"saved: {out}")
    if Path(out).is_file():
        size_mb = Path(out).stat().st_size / 1024 / 1024
        print(f"  size: {size_mb:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
