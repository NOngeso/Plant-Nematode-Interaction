#!/usr/bin/bash -l

# Load relevant modules
module load hisat2/2.1.0
module load samtools/1.9

# reference genome
genome=~/nehe/data/ref_genome/GCA_900182535.1_Meloidogyne_incognita_V3_genomic.fna

#Build the reference genome  needed for Hisat2 aligner:

hisat2-build $genome

#Path ro indexed reference genome
igenome=~/nehe/data/ref_genome/GCA_900182535.1_Meloidogyne_incognita_V3_genomic

for fq1 in ~/nehe/data/*_1.fastq
do
echo "Working With File $fq1"

base=$(basename $fq1 _1.fastq)
echo "Base Name is $base"

fq1=~/nehe/data/${base}_1.fastq
fq2=~/nehe/data/${base}_2.fastq
sam=~/nehe/results/sam/${base}.aligned.sam
bam=~/nehe/results/bam/${base}.aligned.bam
s_bam=~/nehe/results/bam/${base}.aligned.sorted.bam


# Aligning reads to indexed file
hisat2 -x $igenome -1 $fq1 -2 $fq2 -S ~/nehe/results/sam/${base}.sam

# Converting SAM format to BAM using samtools
samtools view -S -b $sam > $bam
samtools sort -o $s_bam $bam
done


