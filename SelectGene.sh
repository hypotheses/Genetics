#!/bin/bash
## Bhoom Suktitipat
## 2019-12-07
## ------------------------------------------
## Usage
## -----------------------------------------
## > bash SelectGene.sh [yourfile.vcf] (INTERVAL) (OUTDIR)
## Ex: >bash SelectGene.sh [yourfile.vcf] chr5:73980969-7417113 HEXB
## -----------------------------------------
## The script currently has been tested when running  GATKv3.8-1 within a docker container. Using -T SelectVariants command and providing -L interval or an interval list file

set -euxo pipefail

## Working for Docker  "broadinstitute/gatk3:3.8-1" 
## docker run -v gatk_bundle:/gatk_bundle -it broadinstitute/gatk3:3.8-1

REF=/gatk_bundle/hg19/ucsc.hg19.fasta
#GATK=/opt/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar

# Use gatk inside docker located at /usr/GenomeAnalysisTK.jar
GATK=/usr/GenomeAnalysisTK.jar
VARIANT_FILE=$1
# HEXB interval chr5:73980969-74017113
# Padded Interval chr5:73975969-74023113
INTERVAL=${2-'chr5:73975969-74023113'}
OUTDIR=${3-HEXB}
PREFIX=$(basename $1 .vcf)

echo Select variants from $1 and saving output to ${PREFIX}-${INTERVAL}.vcf
mkdir -p ${OUTDIR}

function SelectGene {
  java -jar $GATK \
    -T SelectVariants \
    -R $REF \
    -V $1 \
    -L $INTERVAL \
    -o ${OUTDIR}/${PREFIX}-SelectVariants.vcf
}

SelectGene $VARIANT_FILE