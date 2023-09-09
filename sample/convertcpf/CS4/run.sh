
#!/bin/bash
ofile=test.`date +%Y%m%d_%H-%M-%S`.log 
fname=CS4_ligX_md1_12ns_mod_mp2_631gd.cpf.gz

if [ ! -f $fname ]; then
    bash dl.sh
fi

python -m abmptools.convertcpf -i CS4_ligX_md1_12ns_mod_mp2_631gd.cpf.gz -f 1-299 | tee $ofile
