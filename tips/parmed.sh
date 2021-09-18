#!/bin/bash

topfile=$1

echo """\
gromber $topfile
outparm amber.prmtop""" > parmed.in

parmed < parmed.in
