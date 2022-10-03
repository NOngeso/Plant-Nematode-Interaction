#!/usr/bin/bash
set -x # echo each command to standard out before running it

# Generatining Counts of Expressed Genes

for filename in ~/nehe/results/bam/*.sorted.bam
do

echo "Working With File: $filename"

base=$(basename $filename _1.fastq)
echo "Base Name is: $base"

#Path to  annotation_file
annot=~/nehe/data/GFF_file/GCA_900182535.1_Meloidogyne_incognita_V3_genomic.gff
count=~/nehe/results/counts/${base}.txt


#Counting read expression level for each CDS present in the sorted bam file:
htseq-count $filename $annot -f bam -t CDS > $filename.count

done
