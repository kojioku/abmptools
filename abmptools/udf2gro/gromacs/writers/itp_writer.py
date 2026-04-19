# -*- coding: utf-8 -*-
"""
itp_writer.py
-------------
Stub: writes per-molecule .itp include files (not used in the standard flow,
but provided for extensibility).
"""
from __future__ import annotations
from ....core.system_model import SystemModel
from ._validator import raise_if_cognac_only


class ItpWriter:
    """Optional: writes separate .itp files for each molecule type."""

    def write(self, model: SystemModel, prefix: str) -> None:
        """分子ごとの .itp ファイルを書き出す (未実装)。

        Args:
            model: 中間表現のシステムモデル。
            prefix: 出力ファイル名のプレフィクス。

        Raises:
            ValueError: ``model.ensemble_family == 'cognac_only'`` のとき。
        """
        raise_if_cognac_only(model, kind="itp")
        raise NotImplementedError("ItpWriter is not yet implemented.")
