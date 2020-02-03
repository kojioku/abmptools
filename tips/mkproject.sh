#!/bin/bash
if [ $# -eq 0 ]; then
    echo 'Usage: bash tips/mkproject.sh dir'
    exit
fi
if [ ! -f $1 ]; then
    mkdir $1
fi
cp * $1 2>/dev/null
cp -r tips sample $1
