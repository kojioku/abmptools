# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.calc_sett_io
-----------------------------
cg_segmenter dpdgen_exporter が生成する `{name}_calc_sett` ファイルの読み書き。

`{name}_calc_sett` は Python script で以下を定義する:

    total_num_list = [100000]
    step_list = [100]
    dump_ene = 0
    dump = 0
    density = 3.0
    dt = 0.05
    random = '.true.'
    aij_file = 'aij.dat'
    monomer_file = '{name}_monomer'
    ratio_list = [{'chol': 'solvent'}]
    name_tail_list = ['_dpdgen_from_cg']
    default_file = {...}
    phys_param = {...}
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Union

from .models import CalcSett

logger = logging.getLogger(__name__)


def read_calc_sett(path: Union[str, Path]) -> CalcSett:
    """`{name}_calc_sett` ファイルを読んで :class:`CalcSett` を返す。"""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"calc_sett file not found: {path}")

    ns: Dict[str, Any] = {}
    exec(path.read_text(encoding="utf-8"), ns)  # noqa: S102

    phys_param = ns.get("phys_param", {}) or {}
    box = (
        float(phys_param.get("x", 10.0)),
        float(phys_param.get("y", 10.0)),
        float(phys_param.get("z", 10.0)),
    )

    sett = CalcSett(
        total_num_list=list(ns.get("total_num_list", [100000])),
        step_list=list(ns.get("step_list", [100])),
        monomer_file=str(ns.get("monomer_file", "")),
        aij_file=str(ns.get("aij_file", "aij.dat")),
        box_size=box,
        density=float(ns.get("density", 3.0)),
        dt=float(ns.get("dt", 0.05)),
        ratio_list=list(ns.get("ratio_list", [])),
        name_tail_list=list(ns.get("name_tail_list", [])),
        phys_param=dict(phys_param),
        default_file=dict(ns.get("default_file", {})),
        output_interval=int(ns.get("output_interval_steps", ns.get("dump", 100))),
        aii_val=float(ns.get("aii_val", 25.0)),
    )
    logger.info(
        "read_calc_sett: %s (total_num=%s, step=%s, box=%s, monomer=%s, aij=%s)",
        path, sett.total_num_list, sett.step_list, sett.box_size,
        sett.monomer_file, sett.aij_file,
    )
    return sett
