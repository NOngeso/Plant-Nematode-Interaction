#!/usr/bin/bash

#SBATCH -J Alignment
#SBATCH -c 16
#SBATCH -w compute05
#SBATCH -p batch
#SBATCH --mail-type=All
#SBATCH --mail-user=nehemiahongeso@gmail.com

set -e

cd ~/gitau/data/bovine/align/
mkdir -p sam bam bcf vcf
#relevant modules
#module load hisat2/2.1.0
#module load samtools/1.9
#module load bcftools/1.8

#reference genome

genome=~/gitau/data/bovine/refgenome/GCF_002263795.1_ARS-UCD1.2_genomic.fna

#Build an index of the reference genome  needed for Hisat2 aligner:
hisat2-build $genome ~/gitau/data/bovine/refgenome/GCF_002263795.1_ARS-UCD1.2_genomic

igenome=~/gitau/data/bovine/refgenome/GCF_002263795.1_ARS-UCD1.2_genomic

for filename in ~/gitau/data/bovine/align/*.fastq
do
echo "Working With File $filename"

base=$(basename $filename .fastq)
echo "Base Name is $base"

filename=~/gitau/data/bovine/align/${base}.fastq
sam=~/gitau/data/bovine/align/sam/${base}.aligned.sam
bam=~/gitau/data/bovine/align/bam/${base}.aligned.bam
sorted_bam=~/gitau/data/bovine/align/bam/${base}.aligned.sorted.bam

# Aligning reads to indexed file
hisat2 -x $igenome -U $filename -S ~/gitau/data/bovine/align/sam/${base}.sam
samtools view -S -b $sam > $bam
samtools sort -o $sorted_bam $bam
samtools index $sorted_bam

done

