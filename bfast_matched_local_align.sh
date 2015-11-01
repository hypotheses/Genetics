#!/bin/sh
## Version 06/18/2013
## Match and locally align lot 2 samples
REFERENCE=../TruSeq_Amplicon_Cancer_Panel_Manifest_AFP1_PN15032433.txt.fa
for FILENAME in `ls *.fastq.gz`; do
    PREFIX=$(basename $FILENAME .fastq.gz); 
    if [ ! -e ${PREFIX}.baf ]; then 
        echo "bfast match -t -n 20 -A 0 -l -T /tmp/ -K 100 -M 385 -f $REFERENCE -r ${PREFIX}.fastq.gz -z | bfast localalign -t -n 20 -A 0 -M 10 -f $REFERENCE -M 10 > ${PREFIX}.baf"; 
        bfast match -t -n 20 -A 0 -l -T /tmp/ -K 100 -M 385 -f $REFERENCE -r ${PREFIX}.fastq.gz -z | bfast localalign -t -n 20 -A 0 -M 10 -f $REFERENCE -M 10 > ${PREFIX}.baf
        fi
done
