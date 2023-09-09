#!/bin/bash
# python < in.txt
ofile=test.`date +%Y%m%d_%H-%M-%S`.log
python -m abmptools.convertcpf -i gly5-10.cpf -o gly5-10to23.cpf | tee -a $ofile
python -m abmptools.convertcpf -i gly5-23.cpf -o gly5-23to23.cpf | tee -a $ofile
python -m abmptools.convertcpf -i gly5-4201.cpf -o gly5-4to23.cpf | tee -a $ofile
