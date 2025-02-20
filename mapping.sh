#!/bin/bash

# Edit this path to suit your system
BASE_DIR=/your/top/level/path

LD_LIBRARY_PATH=$BASE_DIR/misc/bamlib
export LD_LIBRARY_PATH

# Models codes for the -ml-model paramter to ModelTraining.py and ModelPredict.py are listed below.
# The meaning of each can be found in our paper:

# RULRS
# RULSVMS
# RULR
# RUSVM
# DTBG
# RUBST
# ADB
# XGB

# general paramters

# Comment/uncomment these depending on which genome is being used

## Tomato - low variant density
FASTA_FILE=S_lycopersicum_chromosomes.2.50.fa
SNPS_VCF=RF_002_SZAXPI009284-57.vcf.gz.snpeff.vcf
SAMPLE_ID=AILSACRAIG
BT2_REF=tomato_sl2.50
FEATURES_FILE=tomato_feats.tab

## Tomato - high variant density
#FASTA_FILE=S_lycopersicum_chromosomes.2.50.fa
#SNPS_VCF=RF_037_SZAXPI008747-46.vcf.gz.snpeff.vcf
#SAMPLE_ID=AILSACRAIG
#BT2_REF=tomato_sl2.50
#FEATURES_FILE=tomato_feats.tab

## Pepper
#FASTA_FILE=korean_pepper_genome.fa
#SNPS_VCF=Taeahn_samtools.raw.vcf
#SAMPLE_ID=CM334DS
#BT2_REF=pepper_genome
#FEATURES_FILE=pepper_feats.tab

## Rice
#FASTA_FILE=rice.chromosomes.fa
#SNPS_VCF=rice_sampled_with_geno.vcf
#SAMPLE_ID=RICE_S1
#BT2_REF=rice_genome
#FEATURES_FILE=rice_feats.tab

###################################################

FREEBAYES_PATH=$BASE_DIR/misc_src/freebayes/bin
ART_PATH=$BASE_DIR/varsim_run/ART/art_bin_MountRainier
FASTA_REF=$BASE_DIR/varsim_run/reference/$FASTA_FILE

VARSIM_PATH=$BASE_DIR/varsim_run
TMP_DIR=/tmp  # Any folder on a separate hard disk to minimize IO to a single disk
PYTHON_SCRIPTS_PATH=$BASE_DIR/PythonProgs


# vcf2diploid options
SV_VCF=$VARSIM_PATH/input_vcfs/ril_sv.vcf
SNPS_VCF=$VARSIM_PATH/input_vcfs/$SNPS_VCF


# ART options
READ_LENGTH=100
TOTAL_COVERAGE=3
SEQUENCER="HS20"
MEAN_FRAG_SIZE=400
SD_FRAG_SIZE=100
ART_REF=$VARSIM_PATH/out/$SAMPLE_ID.fa
ART_OPTIONS="'-nf 1 -ss $SEQUENCER -sp'"

# bowtie2 options

BWREF=$BASE_DIR/varsim_run/reference/$BT2_REF
BWOUTFILE=$VARSIM_PATH/alignments/$SAMPLE_ID.bam

# nrmap options

NRMAP_PATH=$BASEDIR/misc/nrmap
NRMAP_INBAM=$BWOUTFILE
NRMAP_OUTBAM=$SAMPLE_ID.sampled.bam
NRMAP_RECALBAM=$SAMPLE_ID.recal.bam

FEATS_FILE=$VARSIM_PATH/$FEATURES_FILE


############# DO NOT CHANGE ANYTHING BELOW THIS LINE!##################


VARSIM="./varsim.py --reference $FASTA_REF --id $SAMPLE_ID --read_length $READ_LENGTH --mean_fragment_size $MEAN_FRAG_SIZE --sd_fragment_size $SD_FRAG_SIZE --nlanes 1 --total_coverage $TOTAL_COVERAGE --simulator_executable $ART_PATH/art_illumina --art_options $ART_OPTIONS --out_dir out --log_dir log --work_dir work --simulator art  --vcfs $SNPS_VCF --disable_rand_vcf --disable_rand_dgv"

eval $VARSIM

bowtie2 -p 52 --rg-id 1 --rg SM:$SAMPLE_ID --rg LB:1 --rg PL:ILLUMINA -x $BWREF -1 $VARSIM_PATH/out/lane0.read1.fq.gz  -2 $VARSIM_PATH/out/lane0.read2.fq.gz | samtools view -bS - > $TMP_DIR/tmp.bam

samtools sort -o $BWOUTFILE -O bam -@ 40 $TMP_DIR/tmp.bam
rm $TMP_DIR/tmp.bam
samtools index $BWOUTFILE

# extract features from bam file
$NRMAP_PATH/nrmap extractfeats --bamfilein $NRMAP_INBAM --readlength $READ_LENGTH --fastafile $FASTA_REF --mapper BT2 --writeclass 1 > $FULL_FEATS_FILE
