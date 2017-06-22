#!/bin/bash
##Written by Lyron Co Ting Keh 6/19/17

## sh fusion.sh hg19.fa Sample_DLBCL021_Normal.singleindex-deduped.sorted.freq.paired.Q30.txt Sample_DLBCL021-Tumor.singleindex-deduped.sorted.bam

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

clang++ -std=c++11 -stdlib=libc++ editReference.cpp -o editReference #remove this
./editReference $ref_genome $normal_file_name