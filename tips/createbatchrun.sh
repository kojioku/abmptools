#!/bin/bash
stime=1
etime=20
interval=1
script='runab24_fgk202107.sh'
temp='1eo8-gmxtrr-wat4ang-xxx-renamed-MP2-6-31Gd-nbo.ajf'
fhead=${temp%xxx*}
ftail=${temp##*xxx}

for i in `seq $stime $interval $etime`
do
    bash $script $fhead$i$ftail
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
