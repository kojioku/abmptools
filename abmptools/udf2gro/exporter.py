# -*- coding: utf-8 -*-
"""
exporter.py
-----------
Orchestrator: UdfAdapter → SystemModel → Writers → output files.

Usage::

    from abmptools.udf2gro import Exporter
    Exporter().export("input.udf", "output_prefix")
"""
from __future__ import annotations
import os


class Exporter:
    """Thin coordinator that wires together adapter and writers."""

    def export(self, udf_path: str, output_prefix: str) -> int:
        """
        Convert *udf_path* and write Gromacs files with *output_prefix*.

        Returns 0 on success, 1 on error (mirrors original return convention).
        """
        from UDFManager import UDFManager
        from .udf_adapter import UdfAdapter
        from .gromacs.writers.gro_writer import GroWriter
        from .gromacs.writers.top_writer import TopWriter
        from .gromacs.writers.mdp_writer import MdpWriter

        udf = UDFManager(udf_path)
        try:
            model = UdfAdapter(udf).build()
        except RuntimeError as exc:
            print(str(exc))
            return 1
        finally:
            udf = None

        file_top = output_prefix + ".top"
        file_mdp = output_prefix + ".mdp"
        file_gro = output_prefix + ".gro"
        file_ndx = output_prefix + ".ndx"

        print(" --- making top ---")
        TopWriter().write(model, file_top)

        print(" --- making gro ---")
        GroWriter().write(model, file_gro)

        if model.ndx_data is not None:
            self._write_ndx(model, file_ndx)

        print(" --- making mdp ---")
        MdpWriter().write(model, file_mdp)

        print("-----------------------------------------")
        print(" Finished!! ")
        print("-----------------------------------------")
        return 0

    # ------------------------------------------------------------------
    @staticmethod
    def _write_ndx(model, filepath: str) -> None:
        from ..core.system_model import NdxData
        ndx: NdxData = model.ndx_data

        with open(filepath, "w") as f:
            f.write("[ System ]\n")
            buf = ""
            for aid in range(1, ndx.atom_id_max + 1):
                buf += "{:5d}".format(aid) if aid < 10000 else "{:7d}".format(aid)
                if len(buf) > 90:
                    f.write(buf + "\n")
                    buf = ""
            if buf:
                f.write(buf + "\n")

            for gro_name, aid_list in ndx.mol_groups.items():
                f.write("[ {} ]\n".format(gro_name))
                buf = ""
                for aid in aid_list:
                    buf += "{:5d}".format(aid) if aid < 10000 else "{:7d}".format(aid)
                    if len(buf) > 90:
                        f.write(buf + "\n")
                        buf = ""
                if buf:
                    f.write(buf + "\n")

            f.write("[ Constraint ]\n")
            buf = ""
            for aid in ndx.constraint_atom_ids:
                buf += "{:5d}".format(aid) if aid < 10000 else "{:7d}".format(aid)
                if len(buf) > 90:
                    f.write(buf + "\n")
                    buf = ""
            if buf:
                f.write(buf + "\n")
