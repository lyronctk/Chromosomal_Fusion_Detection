#!/bin/bash
##Written by Lyron Co Ting Keh 6/19/17

## bash fusion.sh Sample_DLBCL021_Tumor.sorted.bam Sample_DLBCL021-Tumor.singleindex-deduped.sorted.bam 
## bash fusion.sh test.txt blah blah

if [ -z "$2" ]
then
    echo ""
    echo "Usage: bash $0 <tumor_bam> <tumor_deduped_bam>"
    echo ""
    exit 1
fi


ref_genome=/home/alizadehlab/lyronctk/shared/lyronctk/hg19_bowtie_index/hg19.fa
# normal=$2
bam=$1
deduped_bam=$2


# ./editReference $ref_genome $normal
# bwa index $ref_genome


echo "----Filtering unmapped reads..."
sample_name=$(echo $bam | cut -d '_' -f 2- | cut -d '.' -f 1)
bedtools bamtofastq -i <(samtools view -b -f 4 $bam) -fq $sample_name.unmapped.fq

echo "----Splitting reads..."
./splitUnmapped $sample_name.unmapped.fq $sample_name.left.fq $sample_name.right.fq


#get mapped reads with at least 1 mutation on each end
echo "----Filtering reads with mutations at each end..."
samtools view <(samtools fillmd -e -b $deduped_bam $ref_genome) | cut -f1,4,6,10  > $sample_name.refremoved.txt
length_raw=$(head -1 $sample_name.refremoved.txt | cut -f 3)
length=${length_raw%?}
dummy=$(head -c $length < /dev/zero | tr '\0' '\141' | sed -e 's,a,=,g')
awk -v dummy="$dummy" '$4!=dummy' $sample_name.refremoved.txt | sort -k1,1 > $sample_name.mutatedrows.sorted.txt 

samtools view $deduped_bam > $sample_name.deduped.sam
./findEndMutations $sample_name.mutatedrows.sorted.txt $sample_name.deduped.sam $sample_name.endmutations.sam