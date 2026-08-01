# -*- coding: utf-8 -*-
"""
abmptools.core._subprocess
--------------------------
外部コマンド実行の薄いラッパ。

core.acpype が acpype を呼ぶために使う。 元は CG ビルダー側の内部ヘルパ
だったが、 そちらが非公開リポジトリへ移ったため、 公開側で必要な
CommandError / run_command だけをここに持つ。
"""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Optional, Sequence

logger = logging.getLogger(__name__)



class CommandError(Exception):
    """Raised when an external command fails."""

    def __init__(self, cmd: str, returncode: int, stderr: str):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(
            f"Command failed (exit {returncode}): {cmd}\nstderr: {stderr}"
        )


def run_command(
    cmd: Sequence[str],
    cwd: Optional[Path] = None,
    check: bool = True,
    capture: bool = True,
    stdin_text: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """Run an external command with logging.

    Parameters
    ----------
    cmd
        Command arguments.
    cwd
        Working directory for the child process.
    check
        Raise CommandError on non-zero exit (default True).
    capture
        Capture stdout/stderr as text (default True).
    stdin_text
        Optional text fed via stdin (used by ``gmx genion``/``make_ndx``).
    """
    cmd_list = list(cmd)
    cmd_str = " ".join(cmd_list)
    logger.info("Running: %s", cmd_str)

    result = subprocess.run(
        cmd_list,
        cwd=cwd,
        capture_output=capture,
        text=True,
        input=stdin_text,
    )

    if check and result.returncode != 0:
        logger.error("Command failed: %s", cmd_str)
        if result.stderr:
            logger.error("stderr: %s", result.stderr)
        raise CommandError(cmd_str, result.returncode, result.stderr or "")

    return result
