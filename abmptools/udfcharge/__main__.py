"""CLI: OCTA/COGNAC UDF の per-atom 電荷操作。

    # 単分子 UDF の電荷をバルク UDF の同名分子へ転写
    python -m abmptools.udfcharge transfer --template mol.udf --bulk bulk.udf \
        [--out bulk_charged.udf] [--mol-name NAME] [--mol-index I] \
        [--no-verify-types] [--non-strict]

    # 中和 (Σq≈0) された 1 分子 UDF の電荷を指定形式電荷に復元
    python -m abmptools.udfcharge restore --udf mol.udf --formal-charge 12 \
        [--out mol_q+12.udf] [--mol-index I] [--mol-name NAME]

後方互換: サブコマンド省略時は ``transfer`` として解釈する。
"""

from __future__ import annotations

import argparse
import logging
import sys

from .core import (
    assign_charges_to_bulk,
    read_molecule_charges,
    restore_formal_charge,
)

_SUBCOMMANDS = {"transfer", "restore"}


def _run_transfer(args) -> int:
    tmpl = read_molecule_charges(
        args.template, mol_name=args.mol_name, mol_index=args.mol_index,
    )
    res = assign_charges_to_bulk(
        args.bulk, tmpl, args.out,
        verify_atom_types=not args.no_verify_types,
        strict=not args.non_strict,
    )
    print(f"template : {tmpl.mol_name} (n_atoms={tmpl.n_atoms}, "
          f"net_charge={tmpl.net_charge:+.4f})")
    print(f"assigned : {res.n_molecules_assigned}/{res.n_molecules_total} molecules")
    if res.skipped_indices:
        print(f"skipped  : {res.skipped_indices}")
    print(f"output   : {res.out_path}")
    return 0


def _run_restore(args) -> int:
    res = restore_formal_charge(
        args.udf, args.formal_charge, args.out,
        mol_index=args.mol_index, mol_name=args.mol_name, mode=args.mode,
    )
    detail = (f"λ={res.lam:.8f}" if res.mode == "proportional"
              else f"shift=S/N={res.shift:.8f}")
    print(f"molecule     : {res.mol_name} (n_atoms={res.n_atoms})")
    print(f"mode         : {res.mode}")
    print(f"input total  : {res.input_total:+.6f}  (中和済み想定)")
    print(f"formal charge: {res.formal_charge:+d}  ({detail})")
    print(f"output total : {res.output_total:+.6f}")
    print(f"output       : {res.out_path}")
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # 後方互換: 旧フラット呼び出し (transfer の --template ...) を許容
    if argv and argv[0] not in _SUBCOMMANDS and not argv[0].startswith("-"):
        pass  # 不明な位置引数 → そのまま argparse にエラーさせる
    elif argv and argv[0].startswith("-") and "--template" in argv:
        argv = ["transfer"] + argv

    p = argparse.ArgumentParser(
        prog="python -m abmptools.udfcharge",
        description="OCTA/COGNAC UDF の per-atom 電荷操作 (転写 / 形式電荷復元)。",
    )
    sub = p.add_subparsers(dest="cmd")

    pt = sub.add_parser("transfer",
                        help="単分子 UDF の電荷をバルク UDF の同名分子へ転写")
    pt.add_argument("--template", required=True, help="電荷を持つ単分子 UDF")
    pt.add_argument("--bulk", required=True, help="割り当て先のバルク UDF")
    pt.add_argument("--out", default=None, help="出力 (省略時 <bulk>_charged.udf)")
    pt.add_argument("--mol-name", default=None)
    pt.add_argument("--mol-index", type=int, default=None)
    pt.add_argument("--no-verify-types", action="store_true")
    pt.add_argument("--non-strict", action="store_true")
    pt.add_argument("-v", "--verbose", action="store_true")
    pt.set_defaults(func=_run_transfer)

    pr = sub.add_parser("restore",
                        help="中和済み 1 分子 UDF の電荷を指定形式電荷に復元")
    pr.add_argument("--udf", required=True, help="中和 (Σq≈0) された 1 分子 UDF")
    pr.add_argument("--formal-charge", type=int, required=True,
                    help="目標の形式電荷 (整数)")
    pr.add_argument("--out", default=None, help="出力 (省略時 <udf>_q<±S>.udf)")
    pr.add_argument("--mode", choices=["proportional", "uniform"],
                    default="proportional",
                    help="中和ルール: proportional (|q| 比例、 既定) / uniform (均等分配)")
    pr.add_argument("--mol-name", default=None)
    pr.add_argument("--mol-index", type=int, default=0)
    pr.add_argument("-v", "--verbose", action="store_true")
    pr.set_defaults(func=_run_restore)

    args = p.parse_args(argv)
    if not getattr(args, "cmd", None):
        p.print_help()
        return 2

    logging.basicConfig(
        level=logging.INFO if getattr(args, "verbose", False) else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
