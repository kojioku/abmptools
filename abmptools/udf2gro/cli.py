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
    tau_t = None
    tau_p = None
    barostat = None
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
        elif tok == "--barostat" or tok.startswith("--barostat="):
            _, _, inline = tok.partition("=")
            if inline:
                barostat, i = inline, i + 1
            else:
                if i + 1 >= len(rest):
                    print("ERROR: --barostat の値がありません")
                    raise RuntimeError("")
                barostat, i = rest[i + 1], i + 2
            continue
        elif tok in ("--tau-t", "--tau-p") or tok.startswith(("--tau-t=", "--tau-p=")):
            name, _, inline = tok.partition("=")
            if inline:
                raw, i = inline, i + 1
            else:
                if i + 1 >= len(rest):
                    print("ERROR: {} の値がありません".format(name))
                    raise RuntimeError("")
                raw, i = rest[i + 1], i + 2
            try:
                v = float(raw)
            except ValueError:
                print("ERROR: {} {} を解釈できません (ps)".format(name, raw))
                raise RuntimeError("")
            if name == "--tau-t":
                tau_t = v
            else:
                tau_p = v
            continue
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
        print("Usage: {} in_udf_name output_file_base [--unit all_atom|M,E,L] [--tau-t ps] [--tau-p ps]".format(
            os.path.basename(argv[0])
        ))
        print("  --barostat  圧力浴を指定する (C-rescale / Parrinello-Rahman / Berendsen / no)")
        print("              C-rescale は GROMACS 2021+。安定かつ正しい NPT で拘束とも併用可")
        print("  --tau-t  熱浴の tau_t [ps] を直接指定する (既定は Q から算出)")
        print("  --tau-p  圧力浴の tau_p [ps] を直接指定する (既定 2.0)")
        print("  --unit  UDF が Unit_Parameter を持たないときのスケール。")
        print("          all_atom = 長さ Å / エネルギー kcal/mol (GAFF 系)")
        print("          M,E,L    = Mass[amu], Energy[kJ/mol], Length[nm]")
        raise RuntimeError("")

    in_udf_name, output_file_base = positional

    if barostat is not None:
        # 変換の奥ではなく、 ここで弾く。 綴り違いはそのまま .mdp に流れて
        # grompp で初めて落ちるので、 手元で気付けるようにする。
        from .udf_adapter import canonical_barostat
        try:
            barostat = canonical_barostat(barostat)
        except ValueError as exc:
            print("ERROR: %s" % exc)
            raise RuntimeError("")

    from .exporter import Exporter
    return Exporter().export(in_udf_name, output_file_base,
                             unit_parameter=unit_parameter,
                             tau_t=tau_t, tau_p=tau_p,
                             barostat=barostat)


if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError:
        print("ERROR: Export Gromacs failed.")
        print("may be parameter error.")
