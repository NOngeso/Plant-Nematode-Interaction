# Plant-Nematode-Interaction
Analysis documentation for RNA-seq Data to establish a gene-gene co-expression network of the Root-Knot nematode


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
