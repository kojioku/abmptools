#!/bin/bash

scrdir=$(cd $(dirname $0); pwd)
abmptdir=$scrdir/../..
ofile=test.`date +%Y%m%d_%H-%M-%S`.log 
# fname=6lu7-covhip-100ns-0514.inpdbs.tar.bz2

# if [ ! -f $fname ]; then
#     bash dl.sh
# fi

# benchmark gly5
bash generate.sh gly5.pdb 2>&1 |tee -a $ofile

# benchmark 1L2Y(PJMY9)
bash generate.sh NMR_1L2Y_1_HET_fixed.pdb 2>&1 |tee -a $ofile

# cmd bash batch.sh > test.`date +%Y%m%d_%H-%M-%S`.log 2>error.log

