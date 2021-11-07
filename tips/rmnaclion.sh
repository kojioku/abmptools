#!/bin/bash

for f in `ls *.pdb`
do
    cat $f |grep -v "NA"| grep -v "CL" > ${f%.*}-noion.pdb
done
