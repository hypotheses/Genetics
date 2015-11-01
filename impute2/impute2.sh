#!/bin/sh
## modify from impute.sh to run impute when the phased data file is available
## 10/23/2012
set -e
CHRS=`seq 22 -1 1`
GTOOL=/research/imp_data/users/bhoom/bin/gtool
IMPUTE2=/research/imp_data/users/bhoom/bin/impute2
MAPLOCATION=/research/imp_data/users/bhoom/1kg/genotype/map_location.txt
## /research/imp_data/users/bhoom/1kg/genotype/map_location.txt
## Containing the start end location for each chromosome
##for i in {1..22}; do START=$(head -n 1 *fw${i}.map | awk '{print $4}'); END=$(tail -n 1 *fw${i}.map | awk '{print $4}'); echo -e "$i,$START,$END"; done > map_location.txt

for RACE in aa ea; do
    for CHR in $CHRS ; do
	DIR=/research/imp_data/users/bhoom/1kg/genotype/${RACE}_chr
	PREPHASEDIR=/research/3rd_imp_data/users/bhoom/1kg/preimputeFile/${RACE}
	PHASEDIR=/research/3rd_imp_data/users/bhoom/1kg/phased/${RACE}
	IMPUTEDIR=/research/3rd_imp_data/users/bhoom/1kg/imputed/${RACE}
	PEDFILE=${DIR}/${RACE}_b37_fw${CHR}.ped
	MAPFILE=${DIR}/${RACE}_b37_fw${CHR}.map
	CHRSTART=$(grep -w "^$CHR" $MAPLOCATION | awk -F"," '{print $2}')
	CHREND=$(grep -w "^$CHR" $MAPLOCATION | awk -F"," '{print $3}')
	GENMAP=/research/imp_data/users/bhoom/1kg/phaseV3/genetic_map_chr${CHR}_combined_b37.txt
	HAPFILE=/research/imp_data/users/bhoom/1kg/phaseV3/ALL_1000G_phase1integrated_v3_chr${CHR}_impute.hap.gz
	LEGENDFILE=/research/imp_data/users/bhoom/1kg/phaseV3/ALL_1000G_phase1integrated_v3_chr${CHR}_impute.legend.gz
	START=$CHRSTART
	END=$[CHRSTART+5000000]
	pc=1
	while [[ $END -lt $CHREND ]]; do
	    echo -e "Generating imputing_${RACE}${CHR}_pc${pc}.sh"
	    cat <<EOF > imputing_${RACE}${CHR}_pc${pc}.sh
#!/bin/sh
#BSUB -J ${RACE}${CHR}pc${pc}
#BSUB -o $IMPUTEDIR/${RACE}${CHR}pc${pc}.%J.out
#BSUB -e $IMPUTEDIR/${RACE}${CHR}pc${pc}.%J.err
$IMPUTE2 \\
-m $GENMAP \\
-known_haps_g $PHASEDIR/impute2_${RACE}_chr${CHR}_pc${pc}.phased_haps \\
-h $HAPFILE \\
-l $LEGENDFILE \\
-Ne 15000 \\
-int $START $END \\
-allow_large_regions \\
-o  $IMPUTEDIR/imputed_${RACE}_chr${CHR}_pc${pc}
EOF
	    while [[ ! -e $PHASEDIR/impute2_${RACE}_chr${CHR}_pc${pc}.phased_haps ]]; do
		sleep 30
	    done
	    echo "Imputing _${RACE}${CHR}_pc${pc}"
	    sh imputing_${RACE}${CHR}_pc${pc}.sh > $IMPUTEDIR/${RACE}${CHR}pc${pc}.sh.log
	    pc=$[pc+1]
	    START=$[END+1]
	    END=$[START+5000000]
	    if [[ $END -gt $CHREND ]]; then
		END=$CHREND
	    fi
	done
	cat <<EOF > imputing_${RACE}${CHR}_pc${pc}.sh
#!/bin/sh
#BSUB -J ${RACE}${CHR}pc${pc}
#BSUB -o $IMPUTEDIR/${RACE}${CHR}pc${pc}.%J.out
#BSUB -e $IMPUTEDIR/${RACE}${CHR}pc${pc}.%J.err
$IMPUTE2 \\
-m $GENMAP \\
-known_haps_g $PHASEDIR/impute2_${RACE}_chr${CHR}_pc${pc}.phased_haps \\
-h $HAPFILE \\
-l $LEGENDFILE \\
-Ne 15000 \\
-int $START $END \\
-allow_large_regions \\
-o  $IMPUTEDIR/imputed_${RACE}_chr${CHR}_pc${pc}
EOF
	while [[ ! -e $PHASEDIR/impute2_${RACE}_chr${CHR}_pc${pc}.phased_haps ]]; do
	    sleep 30
	done
	echo "Imputing _${RACE}${CHR}_pc${pc}"
	sh imputing_${RACE}${CHR}_pc${pc}.sh > $IMPUTEDIR/${RACE}${CHR}pc${pc}.sh.log
    done	
done

