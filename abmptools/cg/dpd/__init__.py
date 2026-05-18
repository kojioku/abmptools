# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd
----------------
DPD 系構築サブパッケージ (`abmptools.cg` namespace 配下、 ``cg_segmenter`` の兄弟)。

Routes
------
- **R0** (existing): ``abmptools.fragmenter.cg_segmenter.export_dpdgen()`` —
  CG segment から ``{name}_monomer`` + ``{name}_calc_sett`` 生成。
- **R1** (this package): ``{name}_monomer`` + ``{name}_calc_sett`` + ``aij.dat``
  → Cognac DPD 入力 UDF (``*_uin.udf``)、 そのまま COGNAC で実行可能。
- **R2** (this package): ``{name}_monomer`` + ``aij.dat`` → OCTA GUI 用 ``*.dpm``
  + ``<dir>/monomer-lib/<seg>/Virtual.mom`` + ``#Message.txt``、 J-OCTA で開いて
  GUI 編集 → save → UDF 出力。

Public API
~~~~~~~~~~
- :class:`~abmptools.cg.dpd.models.AijMatrix`
- :class:`~abmptools.cg.dpd.models.MonomerSpec`
- :class:`~abmptools.cg.dpd.models.CalcSett`
- :class:`~abmptools.cg.dpd.models.DpdSystemSpec`
- :func:`~abmptools.cg.dpd.aij_io.read_aij` / :func:`~abmptools.cg.dpd.aij_io.write_aij`
- :func:`~abmptools.cg.dpd.monomer_io.read_monomer`
- :func:`~abmptools.cg.dpd.calc_sett_io.read_calc_sett`
"""
from __future__ import annotations

from .models import AijMatrix, CalcSett, DpdSystemSpec, MonomerSpec
from .aij_io import read_aij, write_aij, aij_to_dict, create_empty_aij
from .monomer_io import (
    read_monomer, read_monomers_dict, assign_particle_names, build_monomer,
)
from .calc_sett_io import read_calc_sett
from .dpm_writer import (
    patch_dpm, propagate_virtual_mom, write_message_txt,
    DEFAULT_PATCH_FIELDS, MESSAGE_TXT_CONTENT,
)
from .udf_writer import write_dpd_udf
from .orchestrator import CGDpdBuilder


def open_panel(builder):
    """Lazy import of :func:`notebook_ui.open_panel` (ipywidgets が必須)."""
    from .notebook_ui import open_panel as _open_panel
    return _open_panel(builder)


__all__ = [
    "AijMatrix", "CalcSett", "DpdSystemSpec", "MonomerSpec",
    "read_aij", "write_aij", "aij_to_dict", "create_empty_aij",
    "read_monomer", "read_monomers_dict", "assign_particle_names",
    "build_monomer",
    "read_calc_sett",
    "patch_dpm", "propagate_virtual_mom", "write_message_txt",
    "DEFAULT_PATCH_FIELDS", "MESSAGE_TXT_CONTENT",
    "write_dpd_udf",
    "CGDpdBuilder",
    "open_panel",
]
