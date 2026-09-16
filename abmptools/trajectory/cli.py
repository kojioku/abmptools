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

    # OCTA / gro2udf 用の 2 点セット。 stage 名はディレクトリから決まるので、
    # amorphous の md/ でも Tg 計算の出力先でも同じ 1 行で通る。
    python -m abmptools.trajectory gen_for_udf
    python -m abmptools.trajectory gen_for_udf --stage prod
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .postprocess import (GmxError, gen_for_udf, gmx_energy, nojump, thin,
                          thin_and_nojump, wrap_pbc)


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
                    "(thin + nojump unwrap + PBC wrap)。 amorphous の "
                    "wrap_pbc / gen_for_udf の共通基盤。",
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
        help="-pbc nojump のみ (frame 数そのまま、 OCTA viewer (GOURMET) 用)",
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

    # gmx energy
    p_en = sub.add_parser(
        "energy",
        help="gmx energy で .edr → .xvg ダンプ (旧 gen_for_udf.sh の半分)",
    )
    p_en.add_argument("--edr", required=True, help="入力 .edr")
    p_en.add_argument("--out", required=True, help="出力 .xvg")
    p_en.add_argument("--terms-max", type=int, default=50,
                      help="energy term 番号の上限 (1..N、 default: 50)")
    p_en.add_argument("--gmx", default="gmx",
                      help="gmx 実行 path (default: PATH 解決)")

    # gen_for_udf (energy + nojump をまとめて)
    p_gu = sub.add_parser(
        "gen_for_udf",
        help="OCTA / gro2udf 用に <stage>_energy.xvg と <stage>_nojump.gro を "
             "まとめて出力 (stage は自動検出)",
    )
    p_gu.add_argument("--stage", default=None,
                      help="<stage>.edr/.tpr/.xtc の basename "
                           "(default: ディレクトリから自動検出)")
    p_gu.add_argument("--dir", dest="directory", default=".",
                      help="stage ファイルのあるディレクトリ (default: cwd)")
    p_gu.add_argument("--ndx", default=None,
                      help="index file (default: 使わない = group 0 は tpr の "
                           "System)。 系の一部だけを UDF にするとき --group と "
                           "セットで指定")
    p_gu.add_argument("--terms-max", type=int, default=50,
                      help="energy term 番号の上限 (1..N、 default: 50)")
    p_gu.add_argument("--group", default="0",
                      help="trjconv の group (default: 0 = System)")
    p_gu.add_argument("--ref", default=None,
                      help="trjconv -s に渡す構造 (default: <stage>.tpr)。 "
                           "古い gmx が新しい tpr を読めないときは "
                           "<stage>.gro に自動で退避する")
    p_gu.add_argument("--nojump-format", default="gro", choices=["gro", "xtc"],
                      help="nojump trajectory の形式 (default: gro)。 xtc は "
                           "10 倍ほど小さいが、 読むのに MDAnalysis が要る "
                           "(MD 環境付属の Python には無いことがある)")
    p_gu.add_argument("--max-frames", dest="max_frames", type=int, default=None,
                      help="出力 trajectory の **合計 frame 数**の上限 "
                           "(default: 間引かない)。 「何本に 1 本か」ではなく "
                           "「合計何枚か」。 割り切れないので実際の枚数は "
                           "これ以下の別の数になり、 実行時に表示する")
    p_gu.add_argument("--gmx", default="gmx",
                      help="gmx 実行 path (default: PATH 解決)")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        if args.cmd == "gen_for_udf":
            res = gen_for_udf(
                stage=args.stage,
                directory=args.directory,
                ndx=args.ndx,
                n_energy_terms=args.terms_max,
                group=args.group,
                nojump_format=args.nojump_format,
                reference=args.ref,
                gmx=args.gmx,
                max_frames=args.max_frames,
            )
            print(f"stage: {res['stage']}")
            if res["ndx"]:
                print(f"index: {res['ndx']}")
            if res.get("reference_fallback"):
                print(f"reference: {res['reference_fallback']} "
                      "(tpr をこの gmx が読めなかったため)")
            for key, label in (("energy", "gmx energy"),
                               ("trajectory", "trjconv -pbc nojump")):
                if res[key] is None:
                    print(f"  (skipped: {label} -- input missing)")
                    continue
                path = Path(res[key])
                size_mb = path.stat().st_size / 1024 / 1024
                print(f"  {path}  ({label}, {size_mb:.1f} MB)")
            # **実際に何枚になったかを必ず出す。** --max-frames は割り切れない
            # ので、 黙っていると「87 枚だった」に後から気付けない。
            if res.get("n_frames") is not None:
                got, want, skip = res["n_frames"], args.max_frames, res["skip"]
                print(f"  frames: {got} "
                      f"(--max-frames {want}, skip {skip})")
                # skip は整数なので枚数は ceil(total/skip) しか取れない。
                # total が want をわずかに超えているだけだと skip=2 に跳ね、
                # 半分近くまで減る (101 枚に --max-frames 100 で 51 枚)。
                # 黙っていると「減りすぎ」に気付けないので、そこだけ言う。
                if want and got * 3 < want * 2:
                    print(f"  note: 指定 {want} に対して {got} 枚。 skip は整数しか"
                          f"取れないため ({skip} で ceil)、 総数が {want} を"
                          f"わずかに超えるときは大きく減る。 全部残すなら "
                          f"--max-frames を外す")
            return 0
        if args.cmd == "energy":
            out = gmx_energy(
                edr=args.edr,
                output=args.out,
                terms=range(1, args.terms_max + 1),
                gmx=args.gmx,
            )
        else:
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
    except (FileNotFoundError, ValueError) as e:
        # ValueError は gen_for_udf の stage 曖昧エラー。 traceback を出すと
        # 「どの stage か言って」という指示が埋もれる。
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
