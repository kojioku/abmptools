#!/bin/bash

ofile=test.`date +%Y%m%d_%H-%M-%S`.log
fname=TrpCage-1.cpf

if [ ! -f $fname ]; then
    bash dl.sh
fi

python -m abmptools.generate_difie -i TrpCage-xxx.cpf -t 1 5 1 -z 1 -s 0 -f 1-20 -v 23 -np 5 | tee $ofile
