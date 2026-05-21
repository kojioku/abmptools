"""Convert the GROMACS xtc trajectory to a COGNAC BDF for abmptools.hbond.

Run after ``run_all.sh`` + ``wrap_pbc.sh`` finish.
Output: ``05_npt_final.bdf`` next to the xtc.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "..", "build")

from abmptools.gro2udf.top_parser import TopParser
from abmptools.gro2udf.top_adapter import TopAdapter
from abmptools.gro2udf.top_exporter import TopExporter
from abmptools.amorphous.trajectory_ingest import frames_from_xtc
import abmptools.gro2udf as g2u

XTC = os.path.join(HERE, "05_npt_final_pbc.xtc")
GRO = os.path.join(BUILD, "system.gro")
TOP = os.path.join(BUILD, "system.top")
MDP = os.path.join(HERE, "05_npt_final.mdp")
TEMPLATE = os.path.join(os.path.dirname(g2u.__file__), "default_template.udf")
OUT_UDF = os.path.join(HERE, "05_npt_final.udf")

if not os.path.exists(XTC):
    sys.exit(f"missing {XTC}; run wrap_pbc.sh first")

raw = TopParser().parse(TOP)
model = TopAdapter().build(raw, GRO, mdp_path=MDP)
frames = frames_from_xtc(
    topology_path=GRO, xtc_path=XTC, top_path=TOP,
)
print(f"loaded {len(frames)} frames from {XTC}")
TopExporter().export_model(model, TEMPLATE, OUT_UDF, frames=frames)
print(f"wrote {OUT_UDF}")
# UDFManager reads both .udf and .bdf; the hbond CLI accepts either.
