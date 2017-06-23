#!/bin/bash
##Written by Lyron Co Ting Keh 6/19/17

## bash fusion.sh hg19.fa Sample_DLBCL021_Normal.singleindex-deduped.sorted.freq.paired.Q30.txt Sample_DLBCL021_Tumor.sorted.bam Sample_DLBCL021-Tumor.singleindex-deduped.sorted.bam 
## bash fusion.sh test.txt blah blah

if [ -z "$3" ]
then
    echo "\nUsage: sh $0 <ref_genome> <normal_freq> <tumor_bam> <tumor_deduped_bam>\n"
    echo "  ref_genome: Reference genome provided will be EDITED\n"
    echo "  normal_freq: Used to substitute normal mutations into the reference genome\n"
    echo "  tumor_bam: Find fusions inside this file, must contain unmapped reads\n"
    exit 1
fi

ref=$1
normal=$2
bam=$3
deduped_bam=$4


# ./editReference $ref $normal
# bwa index $ref


sample_name=$(echo $bam | cut -d '_' -f 2- | cut -d '.' -f 1)
bedtools bamtofastq -i <(samtools view -b -f 4 $bam) -fq $sample_name.unmapped.fq

