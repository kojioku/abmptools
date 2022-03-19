#!/bin/bash
stime=1
etime=100
interval=1
zeropad=5
script='runab23q2_fgk.sh'
temp='XXXI-MMFF-Rxxxlayer5-around_ar6.0.ajf'
fhead=${temp%xxx*}
ftail=${temp##*xxx}

for i in `seq $stime $interval $etime`
do
    num=`printf "%0${zeropad}d" $i`
    bash $script $fhead$num$ftail
done

#!/bin/bash
# function usage {
#   cat <<EOM
# Usage: $(basename "$0") [OPTION]...
#   -h          Display help
#   -s VALUE    A explanation for arg called a
#   -e VALUE    A explanation for arg called b
#   -inter VALUE    A explanation for arg called c
#   -script VALUE    A explanation for arg called d
#   -temp VALUE    A explanation for arg called d
#
# EOM
#
#   exit 2
# }
#
# while getopts ":a:b:c:d:h" optKey; do
#   case "$optKey" in
#     a)
#       echo "-a = ${OPTARG}"
#       ;;
#     b)
#       echo "-b = ${OPTARG}"
#       ;;
#     c)
#       echo "-c = ${OPTARG}"
#       ;;
#     d)
#       echo "-d = ${OPTARG}"
#       ;;
#     '-h'|'--help'|* )
#       usage
#       ;;
#   esac
# done
