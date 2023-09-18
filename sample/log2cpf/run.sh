#!/bin/bash
# python < in.txt
ofile=test.`date +%Y%m%d_%H-%M-%S`.log
python -m abmptools.log2cpf -i gly5-MP2-STO-3G-nbo.log -o gly5-MP2-STO-3G-nbo.cpf | tee -a $ofile
# python -m abmptools.log2cpf -i gly5-HF-STO-3G-resp.log -o gly5-HF-STO-3G-resp.cpf | tee -a $ofile


