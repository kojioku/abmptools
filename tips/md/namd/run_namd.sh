#!/bin/bash

chmod +x namdin.sh
bsub -n $1 ./namdin.sh $1
