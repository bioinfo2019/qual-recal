This repository contains programs, scripts and data for investigating the effects of mapping-quality score calibration on variant detection in low-coverage, whole-genome DNA sequencing data.

The resources in the repository allow for the simulation of diploid plant genomes implanted with variants such as SNP, INDEL and various structural variants. Simulated paried-end reads can be generated from these simulated genomes with the ART read simulation tool. Mapping the simulated reads with Bowtie2 back to their original reference genomes will create BAM files needed for the analysis pipeline described below.

Besides the resources here, a few external tools are needed. The VarSim diploid genome simulation package is available below:

https://github.com/bioinform/varsim/releases/download/v0.7.1/varsim-0.7.1.tar.gz

Setup VarSim as described in the user manual. 

VarSim utilizes the ART Illumina read simulator. Get the latest version here:

https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm




