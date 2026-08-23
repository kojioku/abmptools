# -*- coding: utf-8 -*-
"""
cli.py
------
Command-line entry point for udf2gro.

Maintains full backward compatibility with the original::

    python udf2gro.py <input.udf> <output_prefix>

Can also be called via the installed package::

    python -m abmptools.udf2gro <input.udf> <output_prefix>
"""
from __future__ import annotations
import sys
import os


def main(argv=None):
    """udf2groのコマンドラインエントリポイント。

    UDFファイルをGROMACS形式 (gro/top/mdp) に変換する。
    """
    if argv is None:
        argv = sys.argv

    rest = list(argv[1:])
    unit_parameter = None
    positional = []
    i = 0
    while i < len(rest):
        tok = rest[i]
        if tok == "--unit":
            if i + 1 >= len(rest):
                print("ERROR: --unit の値がありません")
                raise RuntimeError("")
            value = rest[i + 1]
            i += 2
        elif tok.startswith("--unit="):
            value = tok.split("=", 1)[1]
            i += 1
        else:
            positional.append(tok)
            i += 1
            continue
        if value == "all_atom":
            unit_parameter = "all_atom"
        else:
            try:
                unit_parameter = tuple(float(x) for x in value.split(","))
            except ValueError:
                print(f"ERROR: --unit {value} を解釈できません "
                      "(all_atom か Mass,Energy,Length)")
                raise RuntimeError("")

    if len(positional) != 2:
        print("Usage: {} in_udf_name output_file_base [--unit all_atom|M,E,L]".format(
            os.path.basename(argv[0])
        ))
        print("  --unit  UDF が Unit_Parameter を持たないときのスケール。")
        print("          all_atom = 長さ Å / エネルギー kcal/mol (J-OCTA / GAFF 系)")
        print("          M,E,L    = Mass[amu], Energy[kJ/mol], Length[nm]")
        raise RuntimeError("")

    in_udf_name, output_file_base = positional

    from .exporter import Exporter
    return Exporter().export(in_udf_name, output_file_base,
                             unit_parameter=unit_parameter)


if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError:
        print("ERROR: Export Gromacs failed.")
        print("may be parameter error.")
