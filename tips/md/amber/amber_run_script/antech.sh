#$ -l f_node=1                 ## * [資源タイプ名]=[使用ノード数個数]
#$ -l h_rt=24:00:00             ## * 経過時間 -> Max 24時間までです。
#$ -p -4                       ## * ジョブの優先度 *デフォルトは-5です。

# カレントディレクトリに移動
#$ -cwd

# Moduleコマンドの初期化
. /etc/profile.d/modules.sh

# CUDA 環境の読み込み
module load amber

# Intel Compiler環境の読み込み
# module load intel/18.0.1.163
# module load intel/19.0.0.117

# Intel MPI 環境の読み込み * コンパイル時の設定環境を指定
# module load intel-mpi/18.1.163
# module load intel-mpi/19.0.117

antechamber -i ${1%.*}.pdb -fi pdb -o ${1%.*}.mol2 -fo mol2 -c bcc -s 2
# antechamber -i ${1%.*}.pdb -fi pdb -o ${1%.*}.mol2 -fo mol2 -c bcc -s 2 -nc 1

