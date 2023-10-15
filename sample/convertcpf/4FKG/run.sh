
#!/bin/bash
ofile=test.`date +%Y%m%d_%H-%M-%S`.log 
fname=4FKQ_MD_96.cpf.gz

if [ ! -f $fname ]; then
    bash dl.sh
fi

python -m abmptools.convertcpf -i 4FKQ_MD_96.cpf.gz -v 23 -f 1-299 | tee $ofile
