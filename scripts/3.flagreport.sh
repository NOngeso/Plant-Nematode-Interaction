#!/usr/bin/bash

module load samtools/1.9
# generating a flag report from the sorted file names
# Utilising samtools to learn more about the sorted bam file as well

for filename in ~/nehe/results/bam/*.sorted.bam

do

echo "Flag Report of $filename" | tee -a flag_report.txt

base=$(basename $filename .sorted.bam)

echo "Working on File $base" | tee -a flag_report.txt

samtools flagstat ~/nehe/results/bam/${base}.sorted.bam | tee -a flag_report.txt

done
