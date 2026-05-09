"""PDB座標ファイルとsegment_data設定からABINIT-MP入力ファイル（AJF・PDBセット）を生成する。

球体・直方体・周辺カットなど複数のカットモードに対応し、FMO計算用の入力一式を出力する。

Public API
----------
- :func:`run_pdb2fmo` -- in-process driver (param dict 渡し可、abmptools
  内の他サブパッケージから直接呼び出す経路。`abmptools.crystal.builder.
  CrystalOrchestrator` 等で使用)
- :func:`main`        -- 既存 CLI エントリ (`python -m abmptools.pdb2fmo`
  / `abmp-pdb2fmo`)。挙動は完全保持; 内部で run_pdb2fmo を呼ぶ
"""
import argparse
import ast
import os
from typing import Dict, List, Optional

import abmptools as abmp


def get_args():
    """コマンドライン引数を解析する。"""
    parser = argparse.ArgumentParser(
                prog='pdb2fmo.py', # program name
                usage='python pdb2fmo.py -i xxx.pdb', # program usage
                description='generate ABNITMP input (ajf,pdb set) from orig pdb and segment_data file',
                epilog='end',
                add_help=True,
                )

    # add args
    parser.add_argument('-i', '--incoord',
                        help='coordinate file (pdb)',
                        nargs='*',
                        # action='append',
                        required=True)

    parser.add_argument('-p', '--parameter',
                        help='parameter file',
                        default='input_param')

    parser.add_argument('-xyz', '--xyz',
                        help='xyz mode',
                        action='store_true'
                        )

    parser.add_argument('-noreid', '--norefreshresid',
                        help='refreshresid',
                        action='store_false'
                        )

    parser.add_argument('-noreatm', '--norefreshatmtype',
                        help='refreshatmtype',
                        action='store_false'
                        )

    # get args
    args = parser.parse_args()
    return args


def _load_param_file(path: str) -> Dict:
    """``input_param`` ファイル (Python dict literal) を読み込む。"""
    param_read = {}
    with open(path, 'r') as _f:
        _tree = ast.parse(_f.read())
    for _node in _tree.body:
        if isinstance(_node, ast.Assign) and len(_node.targets) == 1 and isinstance(_node.targets[0], ast.Name):
            param_read[_node.targets[0].id] = ast.literal_eval(_node.value)
    return param_read['param']


def _resolve_oname(base: str, aobj) -> str:
    """cutmode に応じて出力 basename に suffix を付与する。

    Extracted from the historical ``main()`` body (line 94-107 of the
    pre-refactor file). Behaviour unchanged.
    """
    if aobj.cutmode == 'sphere':
        return (
            base + '-' + aobj.cutmode + 'p' + str(aobj.tgtpos[0]) + '_'
            + str(aobj.tgtpos[1]) + '_' + str(aobj.tgtpos[2]) + '_ar'
            + str(aobj.criteria)
        )
    if aobj.cutmode == 'around':
        return base + '-' + aobj.cutmode + '_ar' + str(aobj.criteria)
    if aobj.cutmode == 'cube':
        return (
            base + '-' + aobj.cutmode + 'p' + str(aobj.tgtpos[0]) + '_'
            + str(aobj.tgtpos[1]) + '_' + str(aobj.tgtpos[2]) + '_x'
            + str(aobj.criteria[0]) + '_y' + str(aobj.criteria[1]) + '_z'
            + str(aobj.criteria[2])
        )
    if aobj.cutmode == 'neutral':
        return base + '-' + aobj.cutmode + '_ar' + str(aobj.criteria)
    if aobj.cutmode == 'none':
        return base + '-for_abmp'
    return base


def run_pdb2fmo(
    pdb_files: List[str],
    *,
    param: Optional[Dict] = None,
    parameter: Optional[str] = None,
    xyz: bool = False,
    refresh_resid: bool = True,
    refresh_atmtype: bool = True,
    output_dir: str = "for_abmp",
    verbose: bool = True,
) -> List[str]:
    """In-process driver for the pdb2fmo workflow.

    Mirrors the legacy CLI loop body — one fresh :class:`abmp.setfmo`
    instance per input PDB, ``setrfmoparam`` from a dict, optional
    ``is_xyz`` flag, cutmode-aware ``oname`` resolution, then
    ``getcontact_rmapfmopdb``. Suitable for callers inside the same
    Python process (``abmptools.crystal.builder.CrystalOrchestrator``
    is the canonical user); the CLI ``main()`` is a thin wrapper around
    this function.

    Parameters
    ----------
    pdb_files
        PDB filenames to process (one ``setfmo`` invocation each).
    param
        Param dict (mirrors ``input_param``'s ``param =`` literal). One
        of ``param`` / ``parameter`` is required.
    parameter
        Path to an ``input_param``-style file (legacy form). Used when
        ``param`` is not supplied.
    xyz
        If True, ``aobj.is_xyz = True`` -- emits ``Natom={N}`` + ``&XYZ``
        block in the AJF (full-precision, mirrors ``-xyz`` flag).
    refresh_resid
        Maps to legacy CLI: True = no ``-noreid`` flag (default).
    refresh_atmtype
        Maps to legacy CLI: True = no ``-noreatm`` flag (default).
    output_dir
        Directory under cwd to receive ``for_abmp/*.{ajf,pdb}``.
    verbose
        If True, prints the same stdout block the CLI does.

    Returns
    -------
    List[str]
        Per-input output basenames (post-cutmode suffix), without
        extension. Useful for the orchestrator to locate emitted files.
    """
    if param is None and parameter is None:
        raise ValueError(
            "run_pdb2fmo requires either `param` (dict) or `parameter` (path)."
        )
    if param is None:
        param = _load_param_file(parameter)

    onames: List[str] = []
    for fname in pdb_files:
        base, _ext = os.path.splitext(fname)
        aobj = abmp.setfmo()
        aobj.setrfmoparam(param)
        if refresh_resid:
            aobj.refreshresid = True
        if refresh_atmtype:
            aobj.refreshatmtype = True
        if xyz:
            aobj.is_xyz = True

        if verbose:
            print('--- info ---')
            print('setup mode', aobj.cutmode)
            print('piedaflag', aobj.piedaflag)
            print('cmmflag', aobj.cmmflag)
            print('refreshresid', aobj.refreshresid)
            print('refreshatmtype', aobj.refreshatmtype)

        oname = _resolve_oname(base, aobj)

        if verbose and aobj.cutmode == 'around':
            print('molset', aobj.molname)
            print('solute', aobj.solutes)
        if verbose and len(aobj.ionname) != 0:
            print('ion mode: ', aobj.ionmode)
            print('ion name: ', aobj.ionname)

        aobj.getcontact_rmapfmopdb(output_dir, fname, oname)
        onames.append(oname)
    return onames


def main():
    """CLI entry: thin wrapper around :func:`run_pdb2fmo`."""
    args = get_args()

    print('coord(pdb) =', args.incoord)
    print('parameter = ', args.parameter)

    run_pdb2fmo(
        pdb_files=list(args.incoord),
        parameter=args.parameter,
        xyz=args.xyz,
        refresh_resid=args.norefreshresid,
        refresh_atmtype=args.norefreshatmtype,
        verbose=True,
    )


if __name__ == "__main__":
    main()
