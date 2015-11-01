#!/bin/bash
set -e
MAPLOCATION=/research/imp_data/users/bhoom/1kg/genotype/map_location.txt
CHRS=`seq 22 -1 1`
## /research/imp_data/users/bhoom/1kg/genotype/map_location.txt
## Containing the start end location for each chromosome
##for i in {1..22}; do START=$(head -n 1 *fw${i}.map | awk '{print $4}'); END=$(tail -n 1 *fw${i}.map | awk '{print $4}'); echo -e "$i,$START,$END"; done > map_location.txt

for CHR in $CHRS; do
    PREPHASEDIR=/research/3rd_imp_data/users/bhoom/1kg/preimputeFile/${RACE}
    PHASEDIR=/research/3rd_imp_data/users/bhoom/1kg/phased/${RACE}
    PEDFILE=${DIR}/${RACE}_b37_fw${CHR}.ped
    MAPFILE=${DIR}/${RACE}_b37_fw${CHR}.map
    CHRSTART=$(grep -w "^$CHR" $MAPLOCATION | awk -F"," '{print $2}')
    CHREND=$(grep -w "^$CHR" $MAPLOCATION | awk -F"," '{print $3}')
    GENMAP=/research/imp_data/users/bhoom/1kg/phaseV3/genetic_map_chr${CHR}_combined_b37.txt
    START=$CHRSTART
    END=$[CHRSTART+5000000]
    pc=1
    while [[ $END -lt $CHREND ]]; do
	echo "$CHR,$pc,$START,$END"
	pc=$[pc+1]
	START=$[END+1]
	END=$[START+5000000]
	if [[ $END -gt $CHREND ]]; then
	    END=$CHREND
	fi
    done
done

