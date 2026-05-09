#!/usr/bin/env bash
# csp7 single-structure smoke pipeline (Phase A/B: legacy CLI invocation).
#
# pdb2fmo runs in `-xyz` mode so the AJF embeds &XYZ coordinates
# directly with full floating-point precision (e.g. 2.891228423347176)
# instead of going through the %8.3f truncation that an external PDB
# ReadGeom path would impose. `pdb_io.exportardxyzfull` reads
# coordinates from a sibling <base>.xyz next to the PDB, so we copy
# the XYZ file emitted by readcif into the PDB directory before
# invoking pdb2fmo.
#
# Phase C will replace this 4-step recipe with a single
#   `python -m abmptools.crystal pipeline --config crystal_config.yaml`
# but until then the legacy CLIs are the supported entry point.
set -euo pipefail

CIF=XXXI-MMFF-R00001.cif
ATOMS_PER_MOL=32
LAYER=5
BASE="XXXI-MMFF-R00001layer${LAYER}Zp1"

# The csp7 cif and UNK.ajf are private and live in abmptools-sample.
# Stage them in this directory if they aren't already present.
SAMPLE_DIR="${ABMPTOOLS_SAMPLE_DIR:-${HOME}/repos/abmptools-sample}"
SRC="${SAMPLE_DIR}/sample/csp7_ciftest/crystal_reference/R00001_layer3_hf_631g"
for f in "${CIF}" UNK.ajf; do
    if [ ! -f "${f}" ]; then
        if [ -f "${SRC}/${f}" ]; then
            cp "${SRC}/${f}" .
            echo "Staged ${f} from ${SRC}/"
        else
            echo "ERROR: ${f} not found here or at ${SRC}/${f}" >&2
            echo "       Set ABMPTOOLS_SAMPLE_DIR to your abmptools-sample checkout." >&2
            exit 1
        fi
    fi
done

# Stage 1: CIF -> supercell PDB + XYZ.
python -m abmptools.readcif -i "${CIF}" -an "${ATOMS_PER_MOL}" -l "${LAYER}"

PDB_DIR="cifout/layer${LAYER}/pdb"
XYZ_DIR="cifout/layer${LAYER}/xyz"

# Stage 2 prep: drivers + matching XYZ next to the layer-5 PDB.
cp input_param segment_data.dat UNK.ajf "${PDB_DIR}/"
cp "${XYZ_DIR}/${BASE}.xyz" "${PDB_DIR}/"

# Stage 2: PDB -> for_abmp/ (-xyz mode = direct &XYZ block in AJF).
cd "${PDB_DIR}"
python -m abmptools.pdb2fmo -i "${BASE}.pdb" -p input_param -xyz

echo
echo "Smoke run complete. FMO inputs are in ${PDB_DIR}/for_abmp/."
