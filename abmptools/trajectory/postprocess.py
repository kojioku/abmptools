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

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Sequence, Union

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


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
            f"{_hint_for(stderr, stdout)}"
        )


def _hint_for(stderr: str, stdout: str) -> str:
    """gmx の失敗に、 直し方が分かる一言を添える。

    ``residuetypes.dat not found`` は「壊れている」のではなく **GMXLIB が
    設定されていない**だけ。 J-OCTA が同梱する gmx は share/ を持っている
    のに環境変数を設定してくれないので、 ``.gro`` を ``-s`` に渡した
    とたんここで止まる。 メッセージ自体は GMXLIB に触れているが、 どこを
    指せばよいかは書いていない。
    """
    text = (stderr or "") + (stdout or "")
    if "residuetypes.dat" in text:
        return (
            "\n--- hint ---\n"
            "gmx could not find its share/top data. Point GMXLIB at it:\n"
            "  cmd : set \"GMXLIB=C:\\J-OCTA-12.0\\additional\\GROMACS\\share\\top\"\n"
            "  sh  : export GMXLIB=/path/to/gromacs/share/top\n"
            "This bites when -s is a .gro (a .tpr carries the data itself)."
        )
    return ""


def _resolve_gmx(gmx: str) -> str:
    """gmx 実行 path を返す。 PATH に無ければ :class:`FileNotFoundError`。"""
    if Path(gmx).is_absolute() or gmx.startswith((".", "/", "\\")):
        if not Path(gmx).exists():
            raise FileNotFoundError(f"gmx executable not found: {gmx}")
        return gmx
    found = shutil.which(gmx)
    if found is None:
        # 「見つからない」だけ言って終わらない。 どう指すかを、 その人が
        # 使っている呼び方 (CLI か API か) で書く。 Windows では J-OCTA に
        # GROMACS が同梱されているので、 その場所も挙げる。
        raise FileNotFoundError(
            f"'{gmx}' not found in PATH.\n"
            "  CLI: --gmx <path to gmx>\n"
            "  API: gmx=\"<path to gmx>\"\n"
            "  or put the directory holding gmx on PATH.\n"
            "  On Windows, J-OCTA ships one:\n"
            "    C:\\J-OCTA-12.0\\additional\\GROMACS\\bin\\gmx.exe"
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

    OCTA viewer (GOURMET) に出す trajectory はこちら (-pbc mol は分子境界で
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


# ---------------------------------------------------------------------------
# gen_for_udf: OCTA / UDF に渡す 2 点セットを 1 呼び出しで作る
# ---------------------------------------------------------------------------

#: stage 自動検出で優先する名前 (amorphous 5-stage protocol の production)。
PREFERRED_STAGES = ("05_npt_final", "prod", "production")


def _is_tpx_version_error(exc: "GmxError") -> bool:
    """``.tpr`` の版が gmx より新しくて読めなかったか。

    gmx は ``reading tpx file (x.tpr) version 138 with version 119 program``
    と言う。 これは入力が壊れているのではなく、 **この gmx が新しすぎる
    tpr を読めない**というだけなので、 reference を差し替えれば先へ進める。
    """
    text = (exc.stderr or "") + (exc.stdout or "")
    return "reading tpx file" in text and "with version" in text


def find_stage(directory: PathLike = ".") -> str:
    """``directory`` から post-process 対象の stage 名を 1 つ決める.

    stage とは ``<stage>.tpr`` / ``<stage>.edr`` / ``<stage>.xtc`` の共通
    basename のこと。 amorphous protocol なら ``05_npt_final`` だが、 Tg 計算
    でも aggregation でも名前は違うだけで構造は同じなので、 ここでは名前を
    決め打ちせずディレクトリの中身から決める。

    候補が複数あって :data:`PREFERRED_STAGES` でも決まらない場合は、
    ``ValueError`` で候補を列挙する (黙って 1 つ選ぶと、 意図しない stage を
    後処理して気付けない)。
    """
    d = Path(directory)
    if not d.is_dir():
        raise FileNotFoundError(f"directory not found: {d}")
    stems = sorted(
        p.stem for p in d.glob("*.tpr")
        if (d / f"{p.stem}.edr").is_file()
        or (d / f"{p.stem}.xtc").is_file()
        or (d / f"{p.stem}.trr").is_file()
    )
    if not stems:
        raise FileNotFoundError(
            f"no MD stage found in {d} "
            "(a stage needs <name>.tpr plus <name>.edr / .xtc / .trr)"
        )
    if len(stems) == 1:
        return stems[0]
    for preferred in PREFERRED_STAGES:
        if preferred in stems:
            return preferred
    raise ValueError(
        f"several stages found in {d}: {', '.join(stems)}. "
        "Pass stage=... (CLI: --stage) to say which one."
    )


def find_ndx(directory: PathLike = ".") -> Optional[Path]:
    """``directory`` の近くにある index file を探す (無ければ ``None``).

    amorphous の md/ からは ``../build/system.ndx``、 単独ディレクトリに置いた
    run なら ``system.ndx``。

    **:func:`gen_for_udf` はこれを自動では呼ばない。** index を渡すと
    ``trjconv`` の group 番号の意味が変わる (group 0 が tpr の System では
    なく、 その index file の最初の group になる) ため、 隣に置かれた無関係な
    ``.ndx`` を拾うと、 原子の一部だけを切り出した ``.gro`` が無警告で
    出来てしまう。 index が要るのは「系の一部だけを UDF にする」場合だけで、
    そのときは呼ぶ側が明示する。
    """
    d = Path(directory)
    for candidate in (d / ".." / "build" / "system.ndx", d / "system.ndx"):
        if candidate.is_file():
            return candidate.resolve()
    return None


def count_frames(trajectory: PathLike, *, gmx: str = "gmx") -> int:
    """``gmx check`` で trajectory の frame 数を数える.

    ``--max-frames`` のように「合計何枚にするか」を指定されたとき、 skip を
    決めるには総数が要る。 MDAnalysis を使えば読めるが、 **gmx だけで完結
    する道を残す** ためにここでは ``gmx check`` を使う
    (MDAnalysis は abmptools の依存ではなく、 J-OCTA 同梱 Python にも
    入っていない -- ``docs/INSTALL.md`` 参照)。

    ``gmx check`` は集計を stderr に出す。 ``Step`` 行の 1 列目が frame 数。
    """
    exe = _resolve_gmx(gmx)
    proc = subprocess.run(
        [exe, "check", "-f", str(trajectory)],
        capture_output=True, text=True,
    )
    # gmx check は読み終えたあと非ゼロで終わることがあるので、 returncode では
    # なく出力で判断する。
    for line in (proc.stderr + proc.stdout).splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "Step":
            try:
                return int(parts[1])
            except ValueError:
                pass
    raise RuntimeError(
        "could not read a frame count out of 'gmx check -f %s'.\n"
        "stderr tail: %s" % (trajectory, "\n".join(
            (proc.stderr or "").splitlines()[-5:]))
    )


def skip_for_max_frames(total: int, max_frames: int) -> int:
    """合計 ``max_frames`` 枚以内に収まる最小の skip を返す.

    ``gmx trjconv -skip S`` は 0, S, 2S, ... を残すので枚数は
    ``ceil(total / S)``。 これを ``max_frames`` 以下にする最小の ``S`` は
    ``ceil(total / max_frames)``。
    """
    if max_frames < 1:
        raise ValueError("max_frames must be >= 1, got %r" % (max_frames,))
    if total < 1:
        return 1
    return max(1, -(-total // max_frames))


def gen_for_udf(
    *,
    stage: Optional[str] = None,
    directory: PathLike = ".",
    ndx: Optional[PathLike] = None,
    n_energy_terms: int = 50,
    group: str = "0",
    nojump_format: str = "gro",
    reference: Optional[PathLike] = None,
    gmx: str = "gmx",
    max_frames: Optional[int] = None,
) -> dict:
    """OCTA viewer / gro2udf 用の ``.xvg`` + nojump ``.gro`` を書き出す.

    amorphous の ``md/gen_for_udf.py`` が行っていた 2 工程をライブラリ側に
    持ってきたもの。 ``stage`` を与えなければ :func:`find_stage` が
    ディレクトリの中身から決めるので、 amorphous の ``05_npt_final`` でも
    Tg 計算後の構造でも同じ呼び方で通る。

    Parameters
    ----------
    stage
        ``<stage>.edr`` / ``<stage>.tpr`` / ``<stage>.xtc`` の basename。
        ``None`` なら :func:`find_stage` で自動検出。
    directory
        stage ファイルが置かれたディレクトリ (default: cwd)。
    ndx
        index file。 default の ``None`` では index を使わず、 group 0 =
        tpr の System (= 全原子) が出力される。 **自動探索はしない** --
        理由は :func:`find_ndx` を参照。 系の一部だけを UDF にしたいときだけ
        ``ndx`` と ``group`` を明示する (その場合、 下流の ``.top`` も同じ
        部分系である必要がある)。
    n_energy_terms
        ``gmx energy`` に渡す term 番号の上限。
    group
        ``trjconv`` の group (default ``"0"`` = System)。
    nojump_format
        nojump trajectory の出力形式。 ``"gro"`` (default) か ``"xtc"``。

        中身は同じで、 入れ物だけが違う。 ``.xtc`` は 10 倍ほど小さい
        (実測 4.27 MB → 0.41 MB) が、 **下流で読むのに MDAnalysis が要る**
        (``.gro`` は不要)。 J-OCTA 同梱の Python には MDAnalysis が入って
        いないので、 Windows では既定の ``"gro"`` が安全。

        どちらを選んでも ``-pbc nojump`` は通す。 あれは入れ物の話ではなく
        **分子を PBC 境界で分断させない**ための処理で、 省くと OCTA viewer
        でも下流の切り出しでも分子が割れる。
    reference
        ``trjconv -s`` に渡す構造。 default では ``<stage>.tpr``。

        **古い gmx は新しい ``.tpr`` を読めない** (J-OCTA 12.0 同梱は
        GROMACS 2020.4 = tpx v119 までで、 GROMACS 2026 が書いた v138 で
        ``reading tpx file ... with version 119 program`` と落ちる)。 その
        ときは ``<stage>.gro`` に自動で切り替えて警告を出す。 明示したい
        ときはここで指定する。
    gmx
        ``gmx`` 実行 path。
    max_frames
        出力する trajectory の **合計 frame 数**の上限。 ``None`` (default) なら
        間引かない。 skip は :func:`skip_for_max_frames` が決める。

        「何本に 1 本か」ではなく「合計何枚か」を指定する点に注意。 下流の
        UDF の大きさが効くのは枚数のほうなので、 こちらを渡せるようにしてある。
        割り切れないので **実際の枚数は指定値以下の別の数になる**。 返り値の
        ``n_frames`` に入れてあり、 CLI は必ず表示する。

    Returns
    -------
    dict
        ``{"stage": str, "energy": Path|None, "trajectory": Path|None,
        "ndx": Path|None, "n_frames": int|None, "skip": int|None}``。
        tpr が読めず ``.gro`` に退避した場合は ``"reference_fallback"`` に
        その path が入る。

    Raises
    ------
    FileNotFoundError
        stage が見つからない、 または energy も trajectory も作れなかった
        (= 後処理として何も成立していない)。
    """
    d = Path(directory)
    if stage is None:
        stage = find_stage(d)

    result = {"stage": stage, "energy": None, "trajectory": None,
              "ndx": Path(ndx) if ndx else None,
              "n_frames": None, "skip": None}

    edr = d / f"{stage}.edr"
    if edr.is_file():
        result["energy"] = gmx_energy(
            edr=edr,
            output=d / f"{stage}_energy.xvg",
            terms=range(1, n_energy_terms + 1),
            gmx=gmx,
        )

    # .trr を優先 (座標 + 速度)、 無ければ .xtc。
    traj = None
    for ext in ("trr", "xtc"):
        candidate = d / f"{stage}.{ext}"
        if candidate.is_file():
            traj = candidate
            break
    if nojump_format not in ("gro", "xtc"):
        raise ValueError(
            f"nojump_format must be 'gro' or 'xtc', got {nojump_format!r}")
    tpr = d / f"{stage}.tpr"
    out_traj = d / f"{stage}_nojump.{nojump_format}"

    def _write(ref):
        """reference を決め打ちで 1 回書き出す (間引きの有無はここで吸収)。"""
        if max_frames is None:
            return nojump(
                trajectory=traj, tpr=ref, output=out_traj,
                group=group, ndx=ndx, gmx=gmx,
            )
        total = count_frames(traj, gmx=gmx)
        skip = skip_for_max_frames(total, max_frames)
        result["skip"] = skip
        result["n_frames"] = -(-total // skip)      # ceil(total / skip)
        return thin_and_nojump(
            trajectory=traj, tpr=ref, output=out_traj,
            skip=skip, group=group, ndx=ndx, gmx=gmx,
        )

    if traj is not None and (reference is not None or tpr.is_file()):
        ref = Path(reference) if reference is not None else tpr
        try:
            result["trajectory"] = _write(ref)
        except GmxError as exc:
            fallback = d / f"{stage}.gro"
            if (reference is not None or not _is_tpx_version_error(exc)
                    or not fallback.is_file()):
                raise
            # gmx が古くて tpr を読めないだけなら、 .gro を reference に
            # 回せば通る。 -pbc nojump は結合情報を使わないので代用できる。
            # ただし**同じ結果にはならない**: nojump は reference から積み
            # 上げるので、 分子まるごとが別の周期イメージに置かれることが
            # ある (分子が割れることはない)。
            logger.warning(
                "%s could not be read by this gmx (tpx version mismatch); "
                "falling back to %s as the trjconv reference. Molecules stay "
                "whole, but a molecule may sit in a different periodic image "
                "than the tpr route would put it in. Pass reference=... "
                "(CLI: --ref) to choose explicitly.",
                tpr.name, fallback.name,
            )
            result["reference_fallback"] = fallback
            result["trajectory"] = _write(fallback)

    if result["energy"] is None and result["trajectory"] is None:
        raise FileNotFoundError(
            f"nothing to export for stage '{stage}' in {d}: "
            f"need {stage}.edr, or {stage}.tpr plus {stage}.xtc/.trr"
        )
    return result
