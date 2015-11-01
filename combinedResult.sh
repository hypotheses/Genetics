#!/bin/sh
set -e
USAGE="Combine a series of files with the same prefixname following by an index number. Each file also has a header.\nThis script will use the first file's header and add additional file to create one bigger file.\nsh combinedResult.sh [prefix] [ext] [from] [to]"
if [[ ! $# -eq 4 ]]; then echo -e $USAGE; exit ;fi
prefix=$1
ext=$2
from=$3
to=$4
cat ${prefix}${from} |gzip > ${prefix}.${ext}.gz
for i in `seq $[from+1] $to`; do
    echo -ne "\rCombining ${prefix}${i}.${ext}"
    tail -n +2 ${prefix}${i}.${ext} | gzip >> ${prefix}.${ext}.gz
done

