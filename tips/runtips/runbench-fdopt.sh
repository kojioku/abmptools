#!/bin/bash
if [ $# -ne 2 ]; then
    echo "引数が足りません"
    exit 1
fi

pdborig=$1
pdb=${pdborig##*/}
ajforig=$2
ajf=${ajforig##*/}
maindir=`pwd`
# mpi
# i: node
function mpi(){
for i in 1 2 4 6 8 10
do
    mkdir -p "optmpi/node$i"
    # j: np
    for j in 1 2 4
    do
        tgtdir="optmpi/node$i/np$j"
        mkdir -p $tgtdir/PDB
        # generate ajf
        sed -e "s/NP=1/NP=$j/g" $ajforig > $tgtdir/$ajf
        # change node and copy sh
        sed -e "s/node=4/node=$i/g" ./ref/runmpi_bindsv1.sh > $tgtdir/runmpi_bindsv1.sh
        # copy pdb
        cp $pdborig $tgtdir
    done
done
}


function smp(){
# for i in 1 2 4 6 8 10
for i in 1
do
    mkdir -p "optsmp/node$i"
    # j: ompnt
    for j in 1 2 4
    do
        tgtdir="optsmp/node$i/ompnt$j"
        mkdir -p $tgtdir/PDB
        # generate ajf
        mem=$((2000*$j))
        ompnum_t=$((76/$j))
        sed -e "s/Memory=2500/Memory=$mem/g" $ajforig > $tgtdir/$ajf
        # change node and copy sh
        sed -e "s/node=4/node=$i/g" ./ref/runsmp_bindsv1.sh | sed -e "s/proc_per_node=19/proc_per_node=$ompnum_t/g" > $tgtdir/runsmp_bindsv1.sh
        # copy pdb
        cp $pdborig $tgtdir
    done
done
}

function sub(){
# ディレクトリ内の runmpi_bindsv1.sh があるすべてのディレクトリを探します。
find $1 -type f -name $2 | while read -r file; do
    # ディレクトリを取得します。
    dir=$(dirname "$file")

    # ディレクトリに移動します。
    cd "$dir" || exit

    # bash runmpi_bindsv1.sh を実行します。
    bash $2 *ajf

    # 元のディレクトリに戻ります。
    cd - || exit
done
}

mpi
# smp
sub optmpi runmpi_bindsv1.sh
# sub optsmp runsmp_bindsv1.sh







