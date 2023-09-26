#!/bin/bash

ofile=test.`date +%Y%m%d_%H-%M-%S`.log 
fname=CS4_ligX_md1_12ns_mod_mp2_631gd.cpf.gz

if [ ! -f $fname ]; then
    bash dl.sh
fi

# python -m abmptools.generate_difie -i CS4_ligX_md1_xxxns_mod_mp2_631gd.cpf.gz -t 12 18 2 -z 1 -s 0 -f 1-299 -v 23 -np 4 # | tee $ofile
python -m abmptools.generate_difie -i CS4_ligX_md1_xxxns_mod_mp2_631gd.cpf.gz -t 12 50 2 -z 1 -s 0 -f 1-299 -v 23 -np 6 | tee $ofile
