## Recalibrating Probability Estimates Used to Compute Mapping Quality Scores

This repository contains programs, scripts and data for investigating the effects of mapping-quality score recalibration on variant detection in low-coverage, whole-genome DNA sequencing data.

The resources in the repository allow for the simulation of diploid plant genomes implanted with variants such as SNP, INDEL and various structural variants. Simulated paried-end reads can be generated from these simulated genomes with the ART read simulation tool. Mapping the simulated reads with Bowtie2 back to their original reference genomes will create BAM files needed for the analysis pipeline described below.

**Software**

The C++ folder in this repository contain scripts and programs to perform functions for extracting features from SAM (sequence/Alignment Map) files to use in training machine learning models to detect incorrectly aligned reads. The two files in the *bamlib* subfolder depend on the SeqAn sequence analysis C++ library.

https://www.seqan.de/

Those two files should be compiled as a shared object library using the C++14 dialect. The files in the *recal* subfolder should be compiled as an application that links to the bamlib shared object. The *main* routine can be found in the file nrmap.cc.

The programs in the Python folder train models, implement data sampling schemes and make predictions on unseen samples.

Besides the resources here, a few external tools are needed. The VarSim diploid genome simulation package is available below:

https://github.com/bioinform/varsim/releases/download/v0.7.1/varsim-0.7.1.tar.gz

Setup VarSim as described in the user manual. 

VarSim utilizes the ART Illumina read simulator. Get the latest version here:

https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm


**Data**

Data needed to replicate the results in our paper include the tomato (S. lycopersicon) genome assembly version SL2.50, pepper (C. annuum) version 1.55 and rice (Oryza sativa L. ssp.indica) version ASM465v1:

ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_2.50/S_lycopersicum_chromosomes.2.50.fa.gz

ftp://ftp.solgenomics.net/genomes/Capsicum_annuum/C.annuum_cvCM334/assemblies/Pepper.v.1.55.total.chr.gz

https://www.ncbi.nlm.nih.gov/assembly/GCA_000004655.2/

A VCF file containing the SNPs and INDELS used for simulating the tomato genome can be found here:

ftp://ftp.solgenomics.net/genomes/tomato_150/150_VCFs_2.50/RF_002_SZAXPI009284-57.vcf.gz.snpeff.vcf.gz

VCF files containing the structural variants for tomato can be found in the SV folder, and VCF files containing SNPs and INDELs for rice and pepper are located in the VCF folder. VCF files for structural variants, SNPs and small INDELs were used as input to the VarSim program to created the simulated tomato genome, while only SNPs and INDELs were used in the simulated genomes of pepper and rice.


**Analysis Pipeline**

There is a bash script, *mapping.sh*, that drives the process of 

* Simulating genomes and reads
* Mapping reads
* Extracting features

Another bash script, *recal.sh*, handles

* Generating a recalibrated BAM file
* Calling SNPs with the original and recalibrated BAM files
* Comparing the resulting VCF files

Paths to data and programs can be set in the bash scripts.

The *mapping.sh* script is run first. Then the Python program *ModelTraining.py* can be used to generate models for predicting bad reads and printing out statistics (average precsion, Brier score, F1 score, number of true/false positives, etc.)

The Python program *ModelPredict.py* is used to predict misplaced reads using any of the models produced by *ModelTraining.py*. The program also outputs a very large Tab-delimited text file used as input to the program that produces the recalibrated BAM file.

After running *ModelPredict.py*, the TAB-delimited text file is used with the *recal.sh* script to generate a recalibrated BAM file, call SNPs and compare the VCFs.

The output of the VCF comparison is used with the *eval_snps.py* program to evaluate the performance of the recalibration procedure on SNP calling.
