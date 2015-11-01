#!/bin/sh
set -e
GTOOL=/research/imp_data/users/bhoom/bin/gtool
IMPUTE2=/research/imp_data/users/bhoom/bin/impute2
MAPLOCATION=/research/imp_data/users/bhoom/1kg/genotype/map_location.txt
## /research/imp_data/users/bhoom/1kg/genotype/map_location.txt
## Containing the start end location for each chromosome
##for i in {1..22}; do START=$(head -n 1 *fw${i}.map | awk '{print $4}'); END=$(tail -n 1 *fw${i}.map | awk '{print $4}'); echo -e "$i,$START,$END"; done > map_location.txt

for RACE in ea aa; do
    for CHR in `seq 1 22` ; do
	DIR=/research/imp_data/users/bhoom/1kg/genotype/${RACE}_chr
	PREPHASEDIR=/research/3rd_imp_data/users/bhoom/1kg/preimputeFile/${RACE}
	PHASEDIR=/research/3rd_imp_data/users/bhoom/1kg/phased/${RACE}
	PEDFILE=${DIR}/${RACE}_b37_fw${CHR}.ped
	MAPFILE=${DIR}/${RACE}_b37_fw${CHR}.map
	CHRSTART=$(grep -w "^$CHR" $MAPLOCATION | awk -F"," '{print $2}')
	CHREND=$(grep -w "^$CHR" $MAPLOCATION | awk -F"," '{print $3}')
cat <<EOF > phasing_${RACE}${CHR}.sh
#!/bin/sh
$GTOOL -P --ped $PEDFILE --map $MAPFILE --og $PREPHASEDIR/prephased_${RACE}_chr${CHR}.gen --os $PREPHASEDIR/prephased_${RACE}_chr${CHR}.sample --discrete_phenotype 0
EOF
sh phasing_${RACE}${CHR}.sh
    done	
done

