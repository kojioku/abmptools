#!/bin/bash
module load amber 2> /dev/null
process_mdout.perl *.stdout
