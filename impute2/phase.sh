#!/bin/bash
set -e
CHRS=`seq 22 -1 1`
RUN=$2
GTOOL=/research/imp_data/users/bhoom/bin/gtool
IMPUTE2=/research/imp_data/users/bhoom/bin/impute2
MAPLOCATION=/research/imp_data/users/bhoom/1kg/genotype/map_location.txt
## /research/imp_data/users/bhoom/1kg/genotype/map_location.txt
## Containing the start end location for each chromosome
##for i in {1..22}; do START=$(head -n 1 *fw${i}.map | awk '{print $4}'); END=$(tail -n 1 *fw${i}.map | awk '{print $4}'); echo -e "$i,$START,$END"; done > map_location.txt

for RACE in aa ea; do
    for CHR in $CHRS; do
	DIR=/research/imp_data/users/bhoom/1kg/genotype/${RACE}_chr
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
    echo -e "Generating imputing_${RACE}${CHR}_pc${pc}.sh"
    cat <<EOF > phasing_${RACE}${CHR}_pc${pc}.sh
#!/bin/bash
#BSUB -J ${RACE}${CHR}pc${pc}
#BSUB -o $PHASEDIR/${RACE}${CHR}pc${pc}.%J.out
#BSUB -e $PHASEDIR/${RACE}${CHR}pc${pc}.%J.err
$IMPUTE2 \\
-prephase_g \\
-m $GENMAP \\
-g $PREPHASEDIR/prephased_${RACE}_chr${CHR}.gen \\
-Ne 15000 \\
-int $START $END \\
-o $PHASEDIR/impute2_${RACE}_chr${CHR}_pc${pc}.phased
EOF
if [[ $RUN != N ]]; then bsub < phasing_${RACE}${CHR}_pc${pc}.sh ;fi
    pc=$[pc+1]
    START=$[END+1]
    END=$[START+5000000]
    if [[ $END -gt $CHREND ]]; then
	END=$CHREND
    fi
done
cat <<EOF > phasing_${RACE}${CHR}_pc${pc}.sh
#!/bin/bash
#BSUB -J ${RACE}${CHR}pc${pc}
#BSUB -o $PHASEDIR/${RACE}${CHR}pc${pc}.%J.out
#BSUB -e $PHASEDIR/${RACE}${CHR}pc${pc}.%J.err
$IMPUTE2 \\
-prephase_g \\
-m $GENMAP \\
-g $PREPHASEDIR/prephased_${RACE}_chr${CHR}.gen \\
-Ne 15000 \\
-int $START $END \\
-o $PHASEDIR/impute2_${RACE}_chr${CHR}_pc${pc}.phased
EOF
if [[ $RUN != N ]]; then bsub < phasing_${RACE}${CHR}_pc${pc}.sh ; fi
    done	
done

