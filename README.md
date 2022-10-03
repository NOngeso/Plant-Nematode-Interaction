# Plant Nematode Interaction
** This work is the documentation of RNA-seq data analysis to establish a gene-gene co-expression network of the Root-Knot nematode. The work is licensed under [The MIT License](https://opensource.org/licenses/MIT). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

## Overview
This pipeline consists of the processing of transcriptomic data. 
The pipeline was developed while undertaking my postgraduate project, and it includes the following steps:
- Sourcing of Transcriptomic datasets from NCBI
- Checking of data quality, removal of adapters, and low-quality reads from the RNA-seq data
- Alignment of reads to the reference genome
- Counting the coding sequences/genes present in the expression dataset
- Merging the Count files into one expression dataset
- Generate a weighted gene-gene interaction network, modules, and hub genes.
- Visualize the network, modules, and hub genes in Cytoscape.
- Map identified gene cluester to g: Profiler to determine the biological functions of the modules.

