#!/bin/bash

LD_LIBRARY_PATH=/home/eliot/Myfiles/PhD/workspace/bamlib/Release
export LD_LIBRARY_PATH

# general paramters

# Comment/uncomment these depending on which genome is being used

## Tomato
FASTA_FILE=S_lycopersicum_chromosomes.2.50.fa
SNPS_VCF=RF_002_SZAXPI009284-57.vcf.gz.snpeff.vcf
SAMPLE_ID=AILSACRAIG
BT2_REF=tomato_sl2.50
FEATURES_FILE=tomato_feats.tab

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

BASE_DIR=/media/eliot/WD4TB/eliot_files/MyFiles/PhD
FREEBAYES_PATH=$BASE_DIR/misc_src/freebayes/bin
ART_PATH=$BASE_DIR/varsim_run/ART/art_bin_MountRainier
FASTA_REF=$BASE_DIR/varsim_run/reference/$FASTA_FILE

VARSIM_PATH=$BASE_DIR/varsim_run
TMP_DIR=/home/eliot  # Any folder on a separate hard disk to minimize IO to a single disk
PYTHON_SCRIPTS_PATH=$BASE_DIR/PythonProgs


# vcf2diploid options
SV_VCF=$VARSIM_PATH/input_vcfs/ril_sv.vcf
SNPS_VCF=$VARSIM_PATH/input_vcfs/$SNPS_VCF


# ART options
READ_LENGTH=100
TOTAL_COVERAGE=3
SEQUENCER="HS20"
MEAN_FRAG_SIZE=400
SD_FRAG_SIZE=30
ART_REF=$VARSIM_PATH/out/$SAMPLE_ID.fa
ART_OPTIONS="'-nf 1 -ss $SEQUENCER -sp'"

# bowtie2 options

BWREF=$BASE_DIR/varsim_run/reference/$BT2_REF
BWOUTFILE=$VARSIM_PATH/alignments/$SAMPLE_ID.bam

# nrmap options

NRMAP_PATH=/home/eliot/cuda-workspace/cuNRmap/Release
NRMAP_INBAM=$BWOUTFILE
NRMAP_SAMPLEDBAM=$SAMPLE_ID.sampled.bam
NRMAP_RECALBAM=$SAMPLE_ID.recal.bam
NRMAP_SCORESFILE=$BASE_DIR/models/10_scores.tab # Not used for anything now, but needs to be left in

SAMPLED_FEATS_FILE=$BASE_DIR/tomato_sampled_feats.tab
FULL_FEATS_FILE=$VARSIM_PATH/$FEATURES_FILE

############# DO NOT CHANGE ANYTHING BELOW THIS LINE!##################

# recalibrate full bam file
$NRMAP_PATH/cuNRmap -c recal -b $NRMAP_INBAM -o $TMP_DIR/$NRMAP_RECALBAM -f $FASTA_REF -s /media/eliot/WD4TB/eliot_files/MyFiles/PhD/models/tom_recal_scores.tab -d $NRMAP_INBAM -l $READ_LENGTH

samtools sort -o $VARSIM_PATH/alignments/$NRMAP_RECALBAM -O bam -@ 12 $TMP_DIR/$NRMAP_RECALBAM
rm $TMP_DIR/$NRMAP_RECALBAM
samtools index $VARSIM_PATH/alignments/$NRMAP_RECALBAM


#################################### Call and compare SNPs #################################################


$FREEBAYES_PATH/freebayes -f $FASTA_REF -b $VARSIM_PATH/alignments/$NRMAP_RECALBAM -v $VARSIM_PATH/variant_calls/recal.vcf --use-mapping-quality
$FREEBAYES_PATH/freebayes -f $FASTA_REF -b $NRMAP_INBAM -v $VARSIM_PATH/variant_calls/norecal.vcf --use-mapping-quality

java -jar $VARSIM_PATH/VarSim.jar vcfcompare -true_vcf $VARSIM_PATH/out/$SAMPLE_ID.truth.vcf -prefix $SAMPLE_ID.recal $VARSIM_PATH/variant_calls/pepper_recal.vcf
java -jar $VARSIM_PATH/VarSim.jar vcfcompare -true_vcf $VARSIM_PATH/out/$SAMPLE_ID.truth.vcf -prefix $SAMPLE_ID.norecal $VARSIM_PATH/variant_calls/pepper_norecal.vcf 
