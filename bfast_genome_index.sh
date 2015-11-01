#!/bin/sh
CURRENT_DIR=/media/Research/CRC_TruSeq/lot1
cd $CURRENT_DIR

## If you got this error, 
"In function "RGBinaryRead": Fatal Error[OutOfRange]. Variable/Value: original."
## Try converting your input EOL from dos to Unix format
## On ubuntu I used $> fromdos -d filename
## This can be installed using $> sudo apt-get install tofrodos

## Copy reference genome file to ~/Research/CRC_TruSeq folder
## cp -v ./051913-cancer\ panel/051913-Run\ folder/130517_M00468_0007_000000000-A4BTV/Data/Intensities/BaseCalls/Alignment/Reference/TruSeq_Amplicon_Cancer_Panel_Manifest_AFP1_PN15032433.txt.fa .

## convert reference to brg format?
REF_GENOME=../TruSeq_Amplicon_Cancer_Panel_Manifest_AFP1_PN15032433.txt.fa
bfast fasta2brg -f ${REF_GENOME}

## Create indexes with a has width of 14
bfast index -f ${REF_GENOME} -m 1111111111111111111111 -w 14 -i 1 
bfast index -f ${REF_GENOME} -m 1111101110111010100101011011111 -w 14 -i 2
bfast index -f ${REF_GENOME} -m 1011110101101001011000011010001111111 -w 14 -i 3
bfast index -f ${REF_GENOME} -m 10111001101001100100111101010001011111 -w 14 -i 4 
bfast index -f ${REF_GENOME} -m 11111011011101111011111111 -w 14 -i 5 
bfast index -f ${REF_GENOME} -m 111111100101001000101111101110111 -w 14 -i 6 
bfast index -f ${REF_GENOME} -m 11110101110010100010101101010111111 -w 14 -i 7 
bfast index -f ${REF_GENOME} -m 111101101011011001100000101101001011101 -w 14 -i 8 
bfast index -f ${REF_GENOME} -m 1111011010001000110101100101100110100111 -w 14 -i 9
bfast index -f ${REF_GENOME} -m 1111010010110110101110010110111011  -w 14 -i 10

for FILENAME in `ls *.fastq.gz`; do
    PREFIX=$(basename $FILENAME .fastq.gz); 
    if [ ! -e ${PREFIX}.baf ]; then 
        echo "bfast match -t -n 20 -A 0 -l -T /tmp/ -K 100 -M 385 -f $REFERENCE -r ${PREFIX}.fastq.gz -z | bfast localalign -t -n 20 -A 0 -M 10 -f $REFERENCE -M 10 > ${PREFIX}.baf"; 
        #bfast match -t -n 20 -A 0 -l -T /tmp/ -K 100 -M 385 -f $REFERENCE -r ${PREFIX}.fastq.gz -z | bfast localalign -t -n 20 -A 0 -M 10 -f $REFERENCE -M 10 > ${PREFIX}.baf
        bfast postprocess -f $REFERENCE -i ${PREFIX}.baf -a 3 -O 1 | samtools view -bS - | samtools sort - ${PREFIX}
		samtools index ${PREFIX}.bam 
		fi
done

# Filter alignments
bfast postprocess -f ${REF_GENOME} -i ${PREFIX}.baf -a 3 -O 1 | samtools view -bS - ${PREFIX}.sam

## bfast match -t -n 20 -A 0 -l -K 100 -M 385 -f ../TruSeq_Amplicon_Cancer_Panel_Manifest_AFP1_PN15032433.txt.fa -r Blood-DNA.fastq.join | bfast localalign -t -n 2 -A 0 -M 10 -f ../TruSeq_Amplicon_Cancer_Panel_Manifest_AFP1_PN15032433.txt.fa > Blood_join.baf

