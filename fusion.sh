#!/bin/bash
##Written by Lyron Co Ting Keh 6/19/17

## sh fusion.sh hg19.fa Sample_DLBCL021_Normal.singleindex-deduped.sorted.freq.paired.Q30.txt Sample_DLBCL021-Tumor.singleindex-deduped.sorted.bam
## sh fusion.sh test.txt blah blah

if [ -z "$3" ]
then
        echo "\nUsage: sh $0 <ref_genome> <normal_freq> <tumor_bam>\n"
        echo "  ref_genome: Reference genome provided will be EDITED\n"
        echo "  normal_freq: Used to substitute normal mutations into the reference genome\n"
        echo "  tumor_bam: Find fusions inside this file\n"
        exit 1
fi

ref_genome=$1
normal_file_name=$2
bam_file_name=$3

./editReference $ref_genome $normal_file_name