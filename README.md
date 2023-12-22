# poorpipe
poorly built nanopore pipeline. It performs alignment using Minimap2, SV calling using Sniffles and CNVpythor, expansion detection using TRGT, and SNV calling using Deepvariant.
The SNVs are phased using whatshap, and annotated using VEP.

Command line:

	python poorpipe.py <ubam_folder> <Sample_ID>

ubam_folder is a folder containing bam files.

# Dependencies
The pipeline itself requires python 3 and slurmpy. 

# Todo


