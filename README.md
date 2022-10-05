# Plant Nematode Interaction
**This work is the documentation of RNA-seq data analysis to establish a gene-gene co-expression network of the Root-Knot nematode. The work is licensed under [The MIT License](https://opensource.org/licenses/MIT). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

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
- **Reference genome:** 
- **RNA-seq data:** 
- Accession-SRP109232
sratoolkit or 
sra explorer
How the ftp protocol for securely downloading the data looks like:

Downloading one dateset we will use for analysis:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5684404/SRR5684404_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR568/004/SRR5684404/SRR5684404_2.fastq.gz
```

Converting the sra compressed file into fastq.gz format readable by analysis tools.
```

```


- **Metadata File:**
Create metadata file containing this: 
( SRR5684403      SRR5684404      SRR5684405      SRR5684406      SRR5684407      SRR5684408      SRR5684409      SRR5684410      SRR5684411      SRR5684412      SRR5684413  SRR5684414       SRR5684415      SRR5684416      SRR5684417)
Metadata file: show how to obtain metdata file using SRA-selector, provide link to the metadata file used


- Annotation file: 


//initialize folder for the user to use 

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
We will start aligning of reads using one of the samples in our dataset (```SRR5684404```). Later, we will be iterating this whole process on all sample files. \
An example of what a ```hisat2``` command looks like is below. All index files have this as their base name ```GCA_900182535.1_Meloidogyne_incognita_V3_genomic```.
The data had pair end reads, where the forward reads are held in  ```SRR5684404_1.fastq.gz``` and the reverse reads in ```SRR5684404_1.fastq.gz```.

```  
hisat2-build /data/ref_genome/GCA_900182535.1_Meloidogyne_incognita_V3_genomic.fna
```

```
hisat2 -x /data/ref_genome/GCA_900182535.1_Meloidogyne_incognita_V3_genomic \
 -1 /data/SRR5684404_1.fastq.gz -2 /data/SRR5684404_2.fastq.gz \
 -S /results/sam/SRR5684404.sam
``` 

You will see output that starts like this:
![screenshot]()

You can have the preview of the [alignment script](https://github.com/NOngeso/Plant-Nematode-Interaction/blob/main/scripts/1.hisat2_align_pe.sh)

## 4. - Counting the coding sequences/genes present in the expression dataset

use cut -f1,6 to show results of individual count matrix
[]()
[]()
show merged dataset 
