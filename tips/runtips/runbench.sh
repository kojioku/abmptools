#!/bin/bash
if [ $# -ne 1 ]; then
    echo "引数が足りません"
    exit 1
fi

pdborig=$1
pdb=${pdborig##*/}
maindir=`pwd`
# mpi
# i: node
function mpi(){
for i in 1 2 4 6 8 10
do
    mkdir -p "mpi/node$i"
    # j: np
    for j in 1 2 4
    do
        tgtdir="mpi/node$i/np$j"
        mkdir -p $tgtdir
        # generate ajf
        python -m abmptools.generateajf -np $j --method HF --memory 2500 -cmm -i $pdb
        # copy ajf
        mv ${pdb%.*}*.ajf $tgtdir
        # change node and copy sh
        sed -e "s/node=4/node=$i/g" ./ref/runmpi_bindsv1.sh > $tgtdir/runmpi_bindsv1.sh
        # copy pdb
        cp $pdborig $tgtdir
        # run
        cd $tgtdir
        bash runmpi_bindsv1.sh *ajf
        cd $maindir
    done
done
}

function smp(){
for i in 1 2 4 6 8 10
do
    mkdir -p "smp/node$i"
    # j: ompnt 
    for j in 1 2 4
    do
        tgtdir="smp/node$i/ompnt$j"
        mkdir -p $tgtdir
        # generate ajf
        mem=$((2500*$j))
        ompnum_t=$((76/$j))
        python -m abmptools.generateajf -np 1 --method HF --memory $mem -cmm -i $pdb
        # copy ajf
        mv ${pdb%.*}*.ajf $tgtdir
        # change node and copy sh
        sed -e "s/node=4/node=$i/g" ./ref/runsmp_bindsv1.sh | sed -e "s/proc_per_node=19/proc_per_node=$ompnum_t/g" > $tgtdir/runsmp_bindsv1.sh
        # copy pdb
        cp $pdborig $tgtdir
        # run
        cd $tgtdir
        bash runsmp_bindsv1.sh *ajf
        cd $maindir
    done
done
}

# mpi
smp







