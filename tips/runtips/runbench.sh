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
    done
done
}

function smp(){
# for i in 1 2 4 6 8 10
for i in 2
do
    mkdir -p "smp/node$i"
    # j: ompnt 
    for j in 1 2 4
    do
        tgtdir="smp/node$i/ompnt$j"
        mkdir -p $tgtdir
        # generate ajf
        mem=$((2000*$j))
        ompnum_t=$((76/$j))
        python -m abmptools.generateajf -np 1 --method HF --memory $mem -cmm -i $pdb
        # copy ajf
        mv ${pdb%.*}*.ajf $tgtdir
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

# mpi
# sub mpi runmpi_bindsv1.sh

smp
sub smp runsmp_bindsv1.sh








