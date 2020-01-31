#!/bin/bash

LD_LIBRARY_PATH=/home/eliot/Myfiles/PhD/workspace/bamlib/Release
export LD_LIBRARY_PATH

# general paramters

BASE_DIR=/media/eliot/WD4TB/eliot_files/MyFiles/PhD
FREEBAYES_PATH=$BASE_DIR/misc_src/freebayes/bin
ART_PATH=$BASE_DIR/varsim_run/ART/art_bin_MountRainier
#FASTA_REF=$BASE_DIR/varsim_run/reference/korean_pepper_genome.fa
FASTA_REF=$BASE_DIR/varsim_run/reference/S_lycopersicum_chromosomes.2.50.fa
#FASTA_REF=$BASE_DIR/varsim_run/rice/rice.chromosomes.fa
SAMPLE_ID=naive25xsnps
#SAMPLE_ID=RICE_S1
#SAMPLE_ID=CM334DS
VARSIM_PATH=$BASE_DIR/varsim_run
TMP_DIR=/home/eliot
PYTHON_SCRIPTS_PATH=$BASE_DIR/PythonProgs


# vcf2diploid options

SV_VCF=$VARSIM_PATH/input_vcfs/ril_sv.vcf
SNPS_VCF=$VARSIM_PATH/input_vcfs/RF_002_SZAXPI009284-57.vcf.gz.snpeff.vcf
#SNPS_VCF=$VARSIM_PATH/input_vcfs/Taeahn_samtools.raw.vcf
#SNPS_VCF=$VARSIM_PATH/rice/rice_sampled_with_geno.vcf


# ART options
READ_LENGTH=100
TOTAL_COVERAGE=3
SEQUENCER="HS20"
MEAN_FRAG_SIZE=300
SD_FRAG_SIZE=90
ART_REF=$VARSIM_PATH/out/$SAMPLE_ID.fa
ART_OPTIONS="'-nf 1 -ss $SEQUENCER -sp'"

# bowtie2 options

BWREF=$BASE_DIR/varsim_run/reference/tomato_sl2.50
#BWREF=$BASE_DIR/varsim_run/reference/korean_pepper_genome
BWOUTFILE=$VARSIM_PATH/alignments/$SAMPLE_ID.bam

# nrmap options

NRMAP_PATH=/home/eliot/cuda-workspace/cuNRmap/Release
NRMAP_INBAM=$BWOUTFILE
NRMAP_SAMPLEDBAM=$SAMPLE_ID.sampled.bam
NRMAP_RECALBAM=$SAMPLE_ID.recal.bam
NRMAP_SCORESFILE=$BASE_DIR/models/10_scores.tab #$BASE_DIR/naive_scores.tab
#NRMAP_SCORESFILE=/home/eliot/cuda-workspace/cuNRmap/Release/rice_scores.tab
SAMPLED_FEATS_FILE=$BASE_DIR/tomato_sampled_feats.tab
FULL_FEATS_FILE=$BASE_DIR/models/tom_all_feats_4base.tab

#SAMPLED_FEATS_FILE=$BASE_DIR/varsim_run/pepper_sampled_feats.tab
#FULL_FEATS_FILE=$BASE_DIR/varsim_run/pepper_full_feats.tab

############# DO NOT CHANGE ANYTHING BELOW THIS LINE!##################


VARSIM="./varsim.py --reference $FASTA_REF --id $SAMPLE_ID --read_length $READ_LENGTH --mean_fragment_size 240 --sd_fragment_size 120 --nlanes 1 --total_coverage $TOTAL_COVERAGE --simulator_executable $ART_PATH/art_illumina --art_options $ART_OPTIONS --out_dir out --log_dir log --work_dir work --simulator art  --vcfs $SNPS_VCF --disable_rand_vcf --disable_rand_dgv"

eval $VARSIM

#$ART_PATH/art_illumina -i $ART_REF -p -l $READ_LENGTH -f $TOTAL_COVERAGE -m $MEAN_FRAG_SIZE -s $SD_FRAG_SIZE -nf 1 -ss $SEQUENCER -sp -o $VARSIM_PATH/out/simu.read

bowtie2 -p 52 --rg-id 1 --rg SM:$SAMPLE_ID --rg LB:1 --rg PL:ILLUMINA -x $BWREF -1 $VARSIM_PATH/out/lane0.read1.fq.gz  -2 $VARSIM_PATH/out/lane0.read2.fq.gz | samtools view -bS - > $TMP_DIR/tmp.bam

samtools sort -o $BWOUTFILE -O bam -@ 40 $TMP_DIR/tmp.bam
rm $TMP_DIR/tmp.bam
samtools index $BWOUTFILE


# sample bam file
$NRMAP_PATH/cuNRmap -c samplebam -b $NRMAP_INBAM -o $TMP_DIR/$NRMAP_SAMPLEDBAM -f $FASTA_REF -s $NRMAP_SCORESFILE -d $NRMAP_INBAM -l $READ_LENGTH

samtools sort -o $VARSIM_PATH/alignments/$NRMAP_SAMPLEDBAM -O bam -@ 10 $TMP_DIR/$NRMAP_SAMPLEDBAM
rm $TMP_DIR/$NRMAP_SAMPLEDBAM
samtools index $VARSIM_PATH/alignments/$NRMAP_SAMPLEDBAM

# extract features from sampled bam file
$NRMAP_PATH/cuNRmap -c extractfeats -b $VARSIM_PATH/alignments/$NRMAP_SAMPLEDBAM -o $NRMAP_SAMPLEDBAM -f $FASTA_REF -s $NRMAP_SCORESFILE -d $NRMAP_INBAM -l $READ_LENGTH > $SAMPLED_FEATS_FILE

# extract features from full bam file
$NRMAP_PATH/cuNRmap -c extractfeats -b $NRMAP_INBAM -o $NRMAP_SAMPLEDBAM -f $FASTA_REF -s $NRMAP_SCORESFILE -d $NRMAP_INBAM -l $READ_LENGTH > $FULL_FEATS_FILE




#################################### PYTHON CALLS HERE #########################################
#
/hdd1/PhD/PythonProgs/train_lr.py -s $SAMPLED_FEATS_FILE -f $FULL_FEATS_FILE -c 100
#
################################################################################################


# recalibrate full bam file
$NRMAP_PATH/cuNRmap -c recal -b $NRMAP_INBAM -o $TMP_DIR/$NRMAP_RECALBAM -f $FASTA_REF -s /media/eliot/WD4TB/eliot_files/MyFiles/PhD/models/tom_recal_scores.tab -d $NRMAP_INBAM -l $READ_LENGTH

samtools sort -o $VARSIM_PATH/alignments/$NRMAP_RECALBAM -O bam -@ 12 $TMP_DIR/$NRMAP_RECALBAM
rm $TMP_DIR/$NRMAP_RECALBAM
samtools index $VARSIM_PATH/alignments/$NRMAP_RECALBAM


#################################### Call and compare SNPs #################################################


$FREEBAYES_PATH/freebayes -f $FASTA_REF -b $VARSIM_PATH/alignments/$NRMAP_RECALBAM -v $VARSIM_PATH/variant_calls/recal.vcf --use-mapping-quality
$FREEBAYES_PATH/freebayes -f $FASTA_REF -b $NRMAP_INBAM -v $VARSIM_PATH/variant_calls/norecal.vcf --use-mapping-quality

java -jar $VARSIM_PATH/VarSim.jar vcfcompare -true_vcf $VARSIM_PATH/out/$SAMPLE_ID.truth.vcf -prefix pepperrecal $VARSIM_PATH/variant_calls/pepper_recal.vcf
java -jar $VARSIM_PATH/VarSim.jar vcfcompare -true_vcf $VARSIM_PATH/out/$SAMPLE_ID.truth.vcf -prefix peppernorecal $VARSIM_PATH/variant_calls/pepper_norecal.vcf 



#################################### PYTHON CALLS HERE #########################################
#
/media/eliot/WD4TB/eliot_files/MyFiles/PhD/PythonProgs/vcfqual_t-test.py 
#
################################################################################################
