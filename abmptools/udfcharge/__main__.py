"""CLI: 単分子 UDF の電荷をバルク UDF の同名分子へ割り当てる。

    python -m abmptools.udfcharge --template mol.udf --bulk bulk.udf \
        [--out bulk_charged.udf] [--mol-name NAME] [--mol-index I] \
        [--no-verify-types] [--non-strict]
"""

from __future__ import annotations

import argparse
import logging
import sys

from .core import assign_charges_to_bulk, read_molecule_charges


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m abmptools.udfcharge",
        description="単分子 UDF (電荷あり) の per-atom 電荷を、 バルク UDF の同名分子"
                    "すべてへ割り当てる (OCTA/COGNAC electrostatic_Site)。",
    )
    p.add_argument("--template", required=True,
                   help="電荷を持つ単分子 UDF")
    p.add_argument("--bulk", required=True,
                   help="割り当て先のバルク UDF (同名分子が複数、 電荷なし)")
    p.add_argument("--out", default=None,
                   help="出力 UDF (省略時 <bulk>_charged.udf)")
    p.add_argument("--mol-name", default=None,
                   help="template から抽出する分子の Mol_Name (省略時 先頭分子)")
    p.add_argument("--mol-index", type=int, default=None,
                   help="template から抽出する分子の index (--mol-name より優先)")
    p.add_argument("--no-verify-types", action="store_true",
                   help="Atom_Type_Name 列の一致検証を無効化")
    p.add_argument("--non-strict", action="store_true",
                   help="検証不一致を例外でなく warning + skip で処理")
    p.add_argument("-v", "--verbose", action="store_true", help="INFO ログを出力")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )

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


if __name__ == "__main__":
    sys.exit(main())
