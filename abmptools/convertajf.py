"""AJFファイルのバージョン変換を行うCLIツール。

入力AJFファイルを指定バージョンに変換して出力する。
主な変換: &MCP → &RELPOT (v1dd2026)

Usage:
    python -m abmptools.convertajf -i input.ajf -v v1dd2026
    python -m abmptools.convertajf -i input.ajf -v v1dd2026 -o output.ajf
"""
import argparse
import os
import sys

from .abinit_io import abinit_io


def get_args():
    parser = argparse.ArgumentParser(
        description="Convert AJF file between ABINIT-MP versions")
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Input AJF file")
    parser.add_argument("-o", "--output",
                        default=None,
                        help="Output AJF file (default: overwrite input or stdout with --stdout)")
    parser.add_argument("-v", "--version",
                        default="v1dd2026",
                        help="Target version (default: v1dd2026)")
    parser.add_argument("--stdout",
                        action="store_true",
                        help="Output to stdout instead of file")
    parser.add_argument("--dry-run",
                        action="store_true",
                        help="Show diff without writing")
    return parser.parse_args()


def main():
    args = get_args()
    aio = abinit_io()
    converted = aio.convert_ajf(args.input, args.version)

    if args.dry_run:
        with open(args.input, 'r') as fh:
            original = fh.read()
        if original.strip() == converted.strip():
            print("No changes needed.")
        else:
            import difflib
            diff = difflib.unified_diff(
                original.splitlines(keepends=True),
                converted.splitlines(keepends=True),
                fromfile=args.input,
                tofile=f"{args.input} ({args.version})",
            )
            sys.stdout.writelines(diff)
        return

    if args.stdout:
        print(converted, end='')
    else:
        oname = args.output or args.input
        with open(oname, 'w') as fh:
            fh.write(converted)
        print(f"Converted: {args.input} -> {oname} (version: {args.version})")


if __name__ == "__main__":
    main()
