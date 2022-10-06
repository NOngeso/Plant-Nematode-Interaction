# Plant Nematode Interaction
**This work is the documentation of RNA-seq data analysis to establish a gene-gene co-expression network of the Root-Knot nematode _Meloidogyne incognita_ which forages on the roots of over 5,500 species of plants causing crop loss. \
The work is licensed under [The MIT License](https://opensource.org/licenses/MIT). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

## Overview
This pipeline consists of the processing of transcriptomic data. 
The pipeline was developed while undertaking my postgraduate project, and it includes the following steps:
- Sourcing of Datasets from NCBI
- Checking of data quality, removal of adapters, and low-quality reads from the RNA-seq data
- Alignment to the reference genome
- Counting the coding sequences/genes present in the expression dataset
- Merging the Count files into one expression dataset
- Generate a weighted gene-gene interaction network, modules, and hub genes.
- Visualize the network, modules, and hub genes in Cytoscape.
- Map identified gene cluester to g: Profiler to determine the biological functions of the modules.

## 1. Data Acquistion 
The datasets that were needed to run the Gene-coexpression network(GCN) analysis: \
Create folders in your home directory to store data and results:
```
mkdir data
mkdir results
```
- **Reference genome and annotation file:** \
 We acquired the reference and annoation genome of _Meloidogyne incognita_ from the [NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=meloidogyne+incognita) genome database and securely downloaded it using `wget`.
 ```
 # Acquiring reference genome
 wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/132/215/GCA_014132215.1_MINJ2/GCA_014132215.1_MINJ2_genomic.fna.gz
 # Acquiring annotation file
 wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/132/215/GCA_014132215.1_MINJ2/GCA_014132215.1_MINJ2_genomic.gbff.gz
 ```
- **RNA-seq data:** \
We acquired openly stored transcriptomic datasets from the [SRA](https://www.ncbi.nlm.nih.gov/sra/?term=meloidogyne+incognita%5Borgn%5D+SRP109232) database stored under Accession `SRP109232`. \
The [sra explorer](https://sra-explorer.info/) was used to extract `file tranfer protocols` (ftp) to download our dataset securely. \
In this illustration we will download one dateset, `SRR5684404`, the ftp links of the forward and reverse reads:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5684404/SRR5684404_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5684404/SRR5684404_2.fastq.gz
```

- **Metadata File:**
We extracted the metadata file to be used for this analysis from [SRA Run selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP109232&o=acc_s%3Aa). The [metadata](https://github.com/NOngeso/Plant-Nematode-Interaction/blob/main/results/metadata.csv) file contained detailed description of each sample file. 

## 2. Checking quality

## 3. Alignment to a reference genome
We performed read alignment or mapping to ascertain where our reads came from in the genome. [HISAT2](http://daehwankimlab.github.io/hisat2/), a fast and sensitive aligner was used for mapping the next-generation reads to the reference genome. \
The alignment process consists of two steps:
- **Indexing the reference genome** \
Indexing of the reference genome was done to speed up the alignment process by enabling HISAT2 to discover potential alignment locations for query sequences.
```  
hisat2-build /data/ref_genome/GCA_900182535.1_Meloidogyne_incognita_V3_genomic.fna
``` 
- **Aligning the reads to the reference genome** \
We will start aligning of reads using one of the samples in our dataset (`SRR5684404`). Later, we will be iterating this whole process on all sample files. \
An example of what a `hisat2` command looks like is below. All index files have this as their base name `GCA_900182535.1_Meloidogyne_incognita_V3_genomic`.
The data had pair end reads, where the forward reads are held in  `SRR5684404_1.fastq.gz` and the reverse reads in `SRR5684404_1.fastq.gz`.

```
hisat2 -x /data/ref_genome/GCA_900182535.1_Meloidogyne_incognita_V3_genomic \
 -1 /data/SRR5684404_1.fastq.gz -2 /data/SRR5684404_2.fastq.gz \
 -S /results/SRR5684404.sam
``` 

You will see output that starts like this:
![screenshot]()

You can have the preview of the [alignment script](https://github.com/NOngeso/Plant-Nematode-Interaction/blob/main/scripts/1.hisat2_align_pe.sh) used to run all the samples.

## 4. Counting the coding sequences/genes present in the expression dataset


The count of sample ```SRR5684404``` has the first 20 ENSEMBL_GeneIDs. count matrix
![SRR5684404 Count](https://github.com/NOngeso/Plant-Nematode-Interaction/blob/main/images/SRR5684404_Count_Matrix.PNG)

Count of each samples are merged into one sample using an [R script](https://github.com/NOngeso/Plant-Nematode-Interaction/blob/main/scripts/4.merge_count.R).
![Merged Counts](https://github.com/NOngeso/Plant-Nematode-Interaction/blob/main/images/Count_Matrix_15_samples.PNG)

You can have the preview of the [count script](https://github.com/NOngeso/Plant-Nematode-Interaction/blob/main/scripts/1.hisat2_align_pe.sh) used to generate gene counts for all samples.
