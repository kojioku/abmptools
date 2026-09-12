# -*- coding: utf-8 -*-
"""
exporter.py
-----------
既にある UDF の**座標とセルだけ**を .gro で差し替える経路。

GROParser -> GROAdapter -> UDFWriter -> UDF。力場・トポロジ・計算条件は
テンプレートの UDF から丸ごと引き継ぐので、**COGNAC で組んだ系を GROMACS で
流して戻す往復**に使う。.top からトポロジごと組み立てたいときは
:class:`~abmptools.gro2udf.top_exporter.TopExporter`。

Usage::

    from abmptools.gro2udf import Exporter
    Exporter().export("test.udf", "output.gro")
    # writes: test_groout.udf  in the current directory

出力名は常に ``{udf_basename}_groout.udf`` で、カレントディレクトリに置く。

Step-filter logic (replicates importStructure / ConvertStructure):
  - read Output_Interval_Steps from the UDF's static section
  - erase all existing records
  - for each GRO frame: write it as a new record only when
      frame.step % output_interval == 0
  - stop after writing more than max_record records (= Total_Steps // interval)
"""
from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)


class Exporter:
    """Thin coordinator: wires together GROParser, GROAdapter, UDFWriter."""

    def export(self, udf_path: str, gro_path: str) -> int:
        """
        Convert *gro_path* frames and write updated UDF to
        ``{udf_basename}_groout.udf`` in the current directory.

        Returns 0 on success, 1 on error.
        """
        try:
            return self._run(udf_path, gro_path)
        except Exception as exc:
            logger.error("%s", exc)
            return 1

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _run(self, udf_path: str, gro_path: str) -> int:
        from UDFManager import UDFManager
        from .gro_parser import GROParser
        from .gro_adapter import GROAdapter
        from .udf_writer import (UDFWriter, set_force_field_comment,
                                 warn_if_template_box_differs,
                                 write_static_cell_abc)

        logger.info("## gro2udf")

        udf = UDFManager(udf_path)

        # Read simulation parameters from static data (before erasing records)
        dt            = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.delta_T", "[ps]"
        )
        total_steps   = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps"
        )
        output_interval = udf.get(
            "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps"
        )
        max_record = int(total_steps // output_interval)

        # Erase all existing records (same as importStructure)
        udf.eraseRecord(0, udf.totalRecord())

        parser  = GROParser()
        adapter = GROAdapter()
        writer  = UDFWriter()

        written = 0
        #: 最初に書いたフレーム。 静的セルの元にする (下記)
        first_frame_cell = None
        for frame in parser.parse_frames(gro_path):
            # Guard: stop if we already wrote more than max_record records
            # (replicates the "while j <= maxRecord" condition)
            if written > max_record:
                break

            # Determine step number (fall back to written count when dt==0)
            steps = frame.step if dt != 0 else written

            # Skip frames not on an output boundary
            if steps % output_interval != 0:
                continue

            udf.newRecord()
            logger.info("steps = %s, record = %s", steps, udf.currentRecord())

            positions, cell = adapter.to_positions_and_cell(frame)
            writer.write_frame(udf, positions, cell, steps, frame.time)
            if first_frame_cell is None:
                first_frame_cell = cell

            written += 1

        logger.info("Total number of records: %s", udf.totalRecord())

        # 静的 Structure.Unit_Cell と Initial_Unit_Cell はテンプレートの値の
        # ままなので、 .gro の箱で揃える。 揃えないと「座標は新しい箱・セルは
        # 古い箱」の UDF になり、 静的側を読む下流 (OCTA の GROMACS
        # コンバータ等) が壊れる。 NVT では箱が変わらないので露見しない。
        if first_frame_cell is not None:
            warn_if_template_box_differs(udf, udf_path, first_frame_cell.a,
                                         first_frame_cell.b, first_frame_cell.c)
            write_static_cell_abc(udf, first_frame_cell.a,
                                  first_frame_cell.b, first_frame_cell.c)

        # 下流は Unit_Parameter.Comment の FF=n で力場を決める。 テンプレート
        # 由来の値があればそれを残す ([[gro2udf]] の既定は GAFF)。
        set_force_field_comment(udf, "gaff")

        # Output file: {udf_basename}_groout.udf in current directory
        output_file = (
            os.path.basename(os.path.splitext(udf_path)[0]) + "_groout.udf"
        )
        logger.info("output file: %s", output_file)
        udf.write(output_file)

        logger.info("Finished!!")
        return 0
