#!/bin/bash

##################
sepindexname="!Ion" #!Ion, System
skipinterval=100
################

if [ $# = 0 ]; then
    echo error! please input argment -p and -x
    exit 1
fi
while getopts f:s: OPTION
do
    case $OPTION in
        f) traj=$OPTARG;;
        s) inigro=$OPTARG;;
        *) exit 1
    esac
done

if [ -z "$traj" ];then
    echo "option f was NOT given, exit."
    exit 1
fi

if [ -z "$inigro" ];then
    echo "option s was NOT given, exit."
    exit 1
fi

echo $traj
echo $inigro

dir='gmxpdbs-foropt'
headbuf=${traj%.*}
head=${headbuf##*/}

# 指定した分子種の削除
# echo "!${rminfo_1}
# q" > conv_rm1.in
# gmx make_ndx -f ${traj%.*}.gro < conv_rm1.in

# セパレート
function separate() {
#step0 create .gro for each step
curtime=$(date "+%Y%m%d-%H%M")
if [ -d $dir ]; then
    echo "Backoff! old $dir is renamed to ${dir}.bak.${curtime}"
    mv $dir $dir.bak.${curtime}
    sleep 1
fi

# 作業用ディレクトリの作成
mkdir $dir 2> /dev/null

gmx trjconv -f ${traj} -o ${traj%.*}_.gro -s ${inigro} -skip ${skipinterval} -sep -n index.ndx << EOF
${sepindexname}
EOF

# 各種ファイルの移動とコピー
mv ${traj%.*}_*.gro $dir
cp *.prmtop $dir
cp *.itp $dir
cp *.top $dir
cp *.mdp $dir

echo "$dir/*.gro was generated."
}

separate
