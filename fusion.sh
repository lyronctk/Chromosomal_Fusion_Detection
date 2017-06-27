#!/bin/bash
##Written by Lyron Co Ting Keh 6/19/17

## bash fusion.sh Sample_DLBCL021_Tumor.sorted.bam Sample_DLBCL021-Tumor.singleindex-deduped.sorted.bam 25 20

if [ -z "$4" ]
then
    echo ""
    echo "Usage: bash $0 <tumor_bam> <tumor_deduped_bam> <end_length> <min_mapping_quality>"
    echo ""
    exit 1
fi


bam=$1
deduped_bam=$2
end_length=$3
min_mapping_quality=$4
sample_name=$(echo $bam | cut -d '_' -f 2- | cut -d '.' -f 1)
ref_genome=/data/indexes/hg19.fa


# ./editReference $ref_genome $normal
# bwa index $ref_genome


echo "----Filtering unmapped reads..."
samtools bam2fq <(samtools view -b -f 4 $bam) > $sample_name.unmapped.fq


echo "----Splitting reads..."
./splitUnmapped $sample_name.unmapped.fq $end_length $sample_name.left.fq $sample_name.right.fq


#get mapped reads with at least 1 mutation on each end--S
echo "----Subset mutated rows..."
samtools view <(samtools fillmd -e -b $deduped_bam $ref_genome) | cut -f1,4,6,10  > $sample_name.refremoved.txt
length_raw=$(head -1 $sample_name.refremoved.txt | cut -f 3)
length=${length_raw%?}
dummy=$(head -c $length < /dev/zero | tr '\0' '\141' | sed -e 's,a,=,g')
awk -v dummy="$dummy" '$4!=dummy' $sample_name.refremoved.txt | sort -k1,1 > $sample_name.mutatedrows.sorted.txt 


echo "----Filtering reads with mutations at each end..."
samtools view $deduped_bam > $sample_name.deduped.sam
./findEndMutations $sample_name.mutatedrows.sorted.txt $end_length $sample_name.deduped.sam $sample_name.left.fq $sample_name.right.fq


echo "----Mapping..."
hg19_index=/home/alizadehlab/lyronctk/shared/lyronctk/hg19_index/bowtie2/hg19
bowtie2 --local --no-unal --no-head $hg19_index $sample_name.left.fq -p 8 > $sample_name.remapped.left.txt
bowtie2 --local --no-unal --no-head $hg19_index  $sample_name.right.fq -p 8 > $sample_name.remapped.right.txt


echo "----Compiling candidates"
./candidates $sample_name.remapped.left.txt $sample_name.remapped.right.txt $min_mapping_quality $sample_name.candidates.txt 


rm $sample_name.deduped.sam $sample_name.unmapped.fq $sample_name.left.fq $sample_name.right.fq $sample_name.refremoved.txt $sample_name.mutatedrows.sorted.txt 
echo "----Done"