#!/bin/bash

topfile=$1

echo """\
gromber $topfile
outparm ${topfile%.*}.prmtop""" > parmed.in

parmed < parmed.in
