# -*- coding: utf-8 -*-
"""
itp_writer.py
-------------
Stub: writes per-molecule .itp include files (not used in the standard flow,
but provided for extensibility).
"""
from __future__ import annotations
from ...system_model import SystemModel


class ItpWriter:
    """Optional: writes separate .itp files for each molecule type."""

    def write(self, model: SystemModel, prefix: str) -> None:
        raise NotImplementedError("ItpWriter is not yet implemented.")
