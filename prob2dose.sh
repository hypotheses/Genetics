#!/bin/sh
set -e
## Convert genotype probability files from IMPUTE2 to dosage file in csv format with id,snp1,snp2,snp3,...

infile=$1
outfile=$2
idfile=$3
snpfile=$outfile.premlinfo

## create name pipe to store tmp output file

## extract snp info from input file, assume to be in column 2
if [[ `gzip -t $infile` -eq 0 ]]; then
	echo "Extracting from gzip file $infile"
    gzip -dc $infile | cut -d" " -f1-5 > $outfile.premlinfo
else
	echo "Extracting from text file $infile"
    cut -d" " -f1-5 $infile > $outfile.premlinfo
fi

## call prob2dose.pl to calculate dosage from probability data and
## outfile contain the summary statistics of each imputed markers and tmp/final dose file
perl /research/3rd_imp_data/users/bhoom/1kg/dosage/qc_summary.pl $infile $outfile $idfile $snpfile 









