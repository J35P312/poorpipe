# poorpipe
poorly built nanopore pipeline. It performs alignment using Minimap2, SV calling using Sniffles and CNVpythor, expansion detection using TRGT, and SNV calling using Deepvariant.
The SNVs are phased using whatshap, and annotated using VEP.

Command line:

	python poorpipe.py <ubam_folder> <Sample_ID> <config>

ubam_folder is a folder containing bam files. The config is a json file, the file needs to be edited before running the pipeline.

# Dependencies
The pipeline itself requires python 3 and slurmpy. 

# Download
TRGT file

https://harrietdashnow.com/STRchive/data/hg19.STRchive-disease-loci.TRGT.bed

# Todo


