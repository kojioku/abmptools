"""Cross-platform Python wrappers around ``gmx trjconv``.

設計方針:
- subprocess を ``shell=False`` で呼び、 group 選択 (``echo 0 | ...`` 相当) は
  ``input=`` で stdin に渡す。 Windows でも動作する。
- 全 path を ``pathlib.Path`` で扱い、 forward slash / backslash の差を吸収。
- 失敗時は :class:`GmxError` (stdout/stderr 添付) で raise。

これらは sample 共通の post-process (aggregation の prod.xtc → 解析用 thin
trajectory、 VMD 向け wrap、 OCTA UDF 向け nojump gro) の基本セット。
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional, Sequence, Union

PathLike = Union[str, Path]


class GmxError(RuntimeError):
    """gmx subprocess の non-zero exit を表す。 stderr / stdout を保持。"""

    def __init__(self, cmd: Sequence[str], returncode: int, stdout: str, stderr: str):
        self.cmd = list(cmd)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        super().__init__(
            f"gmx command failed (exit {returncode}): {' '.join(cmd)}\n"
            f"--- stderr ---\n{stderr}\n--- stdout ---\n{stdout}"
        )


def _resolve_gmx(gmx: str) -> str:
    """gmx 実行 path を返す。 PATH に無ければ :class:`FileNotFoundError`。"""
    if Path(gmx).is_absolute() or gmx.startswith((".", "/", "\\")):
        if not Path(gmx).exists():
            raise FileNotFoundError(f"gmx executable not found: {gmx}")
        return gmx
    found = shutil.which(gmx)
    if found is None:
        raise FileNotFoundError(
            f"'{gmx}' not found in PATH. Activate gmxcudaenv or specify gmx=... explicitly."
        )
    return found


def run_trjconv(
    *,
    trajectory: PathLike,
    reference: PathLike,
    output: PathLike,
    group: Union[str, Sequence[str]] = "System",
    ndx: Optional[PathLike] = None,
    extra_args: Sequence[str] = (),
    gmx: str = "gmx",
) -> Path:
    """``gmx trjconv`` の低レベル wrapper.

    Parameters
    ----------
    trajectory
        入力 trajectory (.xtc / .trr / .gro / .pdb)。
    reference
        reference structure (.tpr 推奨、 .gro でも一応動く)。
    output
        出力 file。 parent dir は自動作成。
    group
        trjconv の対話プロンプトに渡す group 名 (or 番号)。 文字列なら 1 行、
        sequence なら順に 1 行ずつ渡す (``-center`` で center_group + output_group
        2 つ要求される場合は ``["Peptide", "System"]`` のように)。
    ndx
        index ファイル (``-n`` フラグに渡す)。 None なら ``-n`` 不指定。
    extra_args
        ``gmx trjconv`` への追加引数 (``["-pbc", "nojump", "-skip", "10"]`` 等)。
    gmx
        ``gmx`` 実行 path (default は PATH 解決)。

    Returns
    -------
    Path
        生成された ``output`` の絶対 path。

    Raises
    ------
    FileNotFoundError
        ``trajectory`` / ``reference`` / ``gmx`` のいずれかが見つからない。
    GmxError
        ``gmx trjconv`` が非ゼロで終了した。
    """
    traj_p = Path(trajectory)
    ref_p = Path(reference)
    out_p = Path(output)

    if not traj_p.is_file():
        raise FileNotFoundError(f"trajectory not found: {traj_p}")
    if not ref_p.is_file():
        raise FileNotFoundError(f"reference structure not found: {ref_p}")
    out_p.parent.mkdir(parents=True, exist_ok=True)

    gmx_exec = _resolve_gmx(gmx)

    cmd: list[str] = [
        gmx_exec, "trjconv",
        "-f", str(traj_p),
        "-s", str(ref_p),
        "-o", str(out_p),
    ]
    if ndx is not None:
        cmd.extend(["-n", str(Path(ndx))])
    cmd.extend(list(extra_args))

    if isinstance(group, str):
        stdin_text = group + "\n"
    else:
        stdin_text = "".join(f"{g}\n" for g in group)

    proc = subprocess.run(
        cmd,
        input=stdin_text.encode(),
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise GmxError(
            cmd=cmd,
            returncode=proc.returncode,
            stdout=proc.stdout.decode(errors="replace"),
            stderr=proc.stderr.decode(errors="replace"),
        )
    return out_p.resolve()


def _default_output(
    trajectory: PathLike, suffix_tag: str, output: Optional[PathLike]
) -> Path:
    """``<stem><suffix_tag><ext>`` を default 出力名に使う."""
    if output is not None:
        return Path(output)
    traj_p = Path(trajectory)
    return traj_p.with_name(f"{traj_p.stem}{suffix_tag}{traj_p.suffix}")


def thin_and_nojump(
    *,
    trajectory: PathLike,
    tpr: PathLike,
    output: Optional[PathLike] = None,
    skip: int = 10,
    group: str = "System",
    ndx: Optional[PathLike] = None,
    gmx: str = "gmx",
) -> Path:
    """``-pbc nojump -skip <N>`` の組合せ (aggregation 系の基本セット).

    100 ns / 100 ps stride (= 1000 frame) の prod.xtc を skip=10 で間引きすると
    100 frame (= 1 ns stride) の trajectory が得られ、 VMD アニメ + 解析 bandwidth
    に最適、 容量は約 1/10 に縮む。

    ``-pbc nojump`` は first frame を基準に各原子を連続追跡し、 aggregation で
    cluster が box 境界を跨いでも分裂表示にならない。 (注意: 出力 trajectory は
    box 範囲を超える絶対座標を持つので、 wrap が必要なら別 pass で :func:`wrap_pbc`
    を呼ぶ。)

    Default output: ``<stem>_nojump_skip<N>.xtc``。
    """
    out_p = _default_output(trajectory, f"_nojump_skip{skip}", output)
    return run_trjconv(
        trajectory=trajectory, reference=tpr, output=out_p,
        group=group, ndx=ndx,
        extra_args=("-pbc", "nojump", "-skip", str(skip)),
        gmx=gmx,
    )


def nojump(
    *,
    trajectory: PathLike,
    tpr: PathLike,
    output: Optional[PathLike] = None,
    group: str = "System",
    ndx: Optional[PathLike] = None,
    gmx: str = "gmx",
) -> Path:
    """``-pbc nojump`` のみ (frame 数そのまま、 unwrap だけ).

    OCTA / J-OCTA Viewer に出す trajectory はこちら (-pbc mol は分子境界で
    瞬間移動して見える)。 amorphous の旧 ``gen_for_udf.sh`` 相当の処理を 1 file
    分実行。
    """
    out_p = _default_output(trajectory, "_nojump", output)
    return run_trjconv(
        trajectory=trajectory, reference=tpr, output=out_p,
        group=group, ndx=ndx,
        extra_args=("-pbc", "nojump"),
        gmx=gmx,
    )


def thin(
    *,
    trajectory: PathLike,
    tpr: PathLike,
    output: Optional[PathLike] = None,
    skip: int = 10,
    group: str = "System",
    ndx: Optional[PathLike] = None,
    gmx: str = "gmx",
) -> Path:
    """``-skip N`` のみ (PBC 処理なし)."""
    out_p = _default_output(trajectory, f"_skip{skip}", output)
    return run_trjconv(
        trajectory=trajectory, reference=tpr, output=out_p,
        group=group, ndx=ndx,
        extra_args=("-skip", str(skip)),
        gmx=gmx,
    )


def gmx_energy(
    *,
    edr: PathLike,
    output: PathLike,
    terms: Sequence[Union[int, str]] = tuple(range(1, 51)),
    gmx: str = "gmx",
) -> Path:
    """``gmx energy`` で .edr から energy term 一覧を .xvg にダンプする.

    Parameters
    ----------
    edr
        入力 ``.edr`` ファイル (mdrun 出力)。
    output
        出力 ``.xvg`` ファイル。 parent dir は自動作成。
    terms
        ``gmx energy`` の対話プロンプトに渡す term 番号 (or 名前) 列。
        Default は ``1..50``: gmx は存在しない index を silently skip する
        ので、 大きめの上限で全 standard term を一括取得できる
        (旧 ``gen_for_udf.sh`` の ``seq 50`` と同じ慣習)。
    gmx
        ``gmx`` 実行 path (default は PATH 解決)。

    Returns
    -------
    Path
        生成された ``output`` の絶対 path。

    Raises
    ------
    FileNotFoundError
        ``edr`` または ``gmx`` が見つからない。
    GmxError
        ``gmx energy`` が非ゼロで終了した。
    """
    edr_p = Path(edr)
    out_p = Path(output)
    if not edr_p.is_file():
        raise FileNotFoundError(f"edr not found: {edr_p}")
    out_p.parent.mkdir(parents=True, exist_ok=True)

    gmx_exec = _resolve_gmx(gmx)
    cmd = [gmx_exec, "energy", "-f", str(edr_p), "-o", str(out_p)]
    stdin_text = "".join(f"{t}\n" for t in terms)

    proc = subprocess.run(
        cmd,
        input=stdin_text.encode(),
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise GmxError(
            cmd=cmd,
            returncode=proc.returncode,
            stdout=proc.stdout.decode(errors="replace"),
            stderr=proc.stderr.decode(errors="replace"),
        )
    return out_p.resolve()


def wrap_pbc(
    *,
    trajectory: PathLike,
    tpr: PathLike,
    output: Optional[PathLike] = None,
    group: str = "System",
    center: Optional[str] = None,
    ur: str = "compact",
    ndx: Optional[PathLike] = None,
    gmx: str = "gmx",
) -> Path:
    """``-pbc mol -ur <ur>`` (+ optional ``-center``) で VMD 用 wrap.

    aggregation で box 跨ぎを起こした trajectory を box 内に戻す。 VMD で
    compact unit-cell として表示するための定番処理。 amorphous の旧
    ``wrap_pbc.sh`` 相当 (1 file 分)。

    Parameters
    ----------
    center
        指定すると ``-center`` を追加し、 該当 group を box 中央に置く。
        gmx trjconv は 2 つの group を要求するので、 stdin に
        ``center\\noutput_group\\n`` を渡す。

    Default output: ``<stem>_pbc.xtc`` (``_pbc`` は amorphous の慣習)。
    """
    out_p = _default_output(trajectory, "_pbc", output)
    extra: list[str] = ["-pbc", "mol", "-ur", ur]
    if center is not None:
        extra.append("-center")
        group_seq = [center, group]
    else:
        group_seq = group  # type: ignore[assignment]
    return run_trjconv(
        trajectory=trajectory, reference=tpr, output=out_p,
        group=group_seq, ndx=ndx,
        extra_args=tuple(extra),
        gmx=gmx,
    )
