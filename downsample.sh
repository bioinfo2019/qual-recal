#!/bin/bash

BASE_DIR=/your/top/level/path

export LD_LIBRARY_PATH=$BASE_DIR/misc/bamlib
export PATH=$PATH:$BASE_DIR/misc/freebayes/bin:$BASE_DIR/misc/freebayes/vcflib/bin
ulimit -Sn 1000000
READLENGTH=100
PREFXI=CUCLW
CUCFASTA=Cucumber_201809_Chr.fa

FASTA_FILE=S_lycopersicum_chromosomes.2.50.fa
SNPS_VCF=RF_037_SZAXPI008747-46.vcf.gz.snpeff.vcf #tom_046.subsample.vcf #RF_002_SZAXPI009284-57.vcf.gz.snpeff.vcf
TOM_SAMPLE_ID=TOM046A #AILSA
BT2_REF=tomato_sl2.50 #cucumber_genome_bt2 #tomato_sl2.50
FEATURES_FILE=cuclw_feats.tab #ailsa_feats.tab

FREEBAYES_PATH=$BASE_DIR/misc_src/freebayes/bin
ART_PATH=$BASE_DIR/varsim_run/ART/art_bin_MountRainier
FASTA_REF=$BASE_DIR/varsim_run/reference/$FASTA_FILE

VARSIM_PATH=$BASE_DIR/varsim_run
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


#---------------------------------------------------#

# Simulate high variant density tomato genome and 100bp reads. Once aligned, features will be extracted and used to train a AdaBoost model to 
# detect misaligned reads. This model will be used to subsequently predict misaligned reads in the real, downsampled cucumber data.

BWOUTFILE=$VARSIM_PATH/alignments/$TOM_SAMPLE_ID.bam
cd $BASE_DIR/varsim_run
echo "simulating tomato genome with high variant density...."
VARSIM="./varsim.py --reference $FASTA_REF --id $SAMPLE_ID --read_length $READ_LENGTH --mean_fragment_size $MEAN_FRAG_SIZE --sd_fragment_size $SD_FRAG_SIZE \
       --nlanes 1 --total_coverage $TOTAL_COVERAGE --simulator_executable $ART_PATH/art_illumina --art_options $ART_OPTIONS --out_dir out --log_dir log \
       --work_dir work --simulator art  --vcfs $SNPS_VCF --disable_rand_vcf --disable_rand_dgv"

eval $VARSIM

echo "aligning simulated tomato reads...."
bowtie2 -p 12 --rg-id 1 --rg SM:$TOM_SAMPLE_ID --rg LB:1 --rg PL:ILLUMINA -x $BWREF -1 $VARSIM_PATH/out/lane0.read1.fq.gz -2 $VARSIM_PATH/out/lane0.read2.fq.gz \
| samtools view -ub - | samtools sort -o $BWOUTFILE -O bam -@ 12 -
samtools index $BWOUTFILE

echo "extracting features...."
cd $BASE_DIR/misc/nrmap
./nrmap extractfeats --bamfilein $BWOUTFILE --readlength 100 --fastafile $FASTA_REF --mapper BT2 --writeclass 1 > $BASEDIR/models/$TOM_SAMPLE_ID.feats.tab
sed -i "s/\t$//" $BASEDIR/models/$TOM_SAMPLE_ID.feats.tab

echo "training model (this will take a LONG time) ...."


cd $BASE_DIR/PythonProgs
./ModelTraining.py --base-dir $BASE_DIR --ml-model ADB --features MAPQ,PAIR_ALIGNMENT_TYPE,SECONDARY_ALIGNMENT_SCORE,N_LOW_QUAL_BASES,ALIGNMENT_SCORES_DIFF \
--feats-file $BASEDIR/models/$TOM_SAMPLE_ID.feats.tab


# Now align the cucumber reads. The reads should be in $BASE_DIR/varsim_run/fastq

echo "aligning 41x cucumber reads...."
BWREF=$BASE_DIR/varsim_run/reference/cucumber_genome_bt2
BWOUTFILE=$VARSIM_PATH/alignments/$PREFIX.BT2BWOUTFILE=$VARSIM_PATH/alignments/$TOM_SAMPLE_ID.bam.bam
bowtie2 -p 12 --rg-id 1 --rg SM:$PREFIX --rg LB:1 --rg PL:ILLUMINA -x $BWREF -1 $VARSIM_PATH/fastq/cuc_1.fq -2 $VARSIM_PATH/fastq/cuc_2.fq \
| samtools view -ub - | samtools sort -o $BWOUTFILE -O bam -@ 12 -
samtools index $BWOUTFILE


declare -a arr=("15" "20" "25" "30" "35" "40" "45" "50" "55" "60" "65" "70" "85" "95")

for f in "${arr[@]}"
do

	# DOWNSAMPLING

	echo "downsampling ...."
	samtools view -s $f.078 -b -o $BASE_DIR/varsim_run/alignments/$PREFIX.$f.bam $BASE_DIR/varsim_run/alignments/$BWOUTFILE

	# FEATURE EXTRACTION

	echo "extracting features...."
	cd $BASE_DIR/misc/nrmap
	./nrmap extractfeats --bamfilein $BWOUTFILE --readlength 100 --fastafile $BASE_DIR/varsim_run/reference/$CUCFASTA --mapper BT2 --writeclass 0 > $BASEDIR/models/$PREFIX.feats.tab

	sed -i "s/\t$//" $BASEDIR/models/$PREFIX.feats.tab

	# PREDICTION

	echo "predicting bad reads...."
	cd $BASE_DIR/PythonProgs
	./ModelTraining.py --base-dir $BASE_DIR --ml-model ADB --features MAPQ,PAIR_ALIGNMENT_TYPE,SECONDARY_ALIGNMENT_SCORE,N_LOW_QUAL_BASES,ALIGNMENT_SCORES_DIFF \
	--feats-file $BASEDIR/models/$PREFIX.feats.tab

	sed "s/'//g;s/\[//;s/\]//" -i $BASE_DIR/models/$PREFIX.scores.tab
	cd $BASE_DIR

	# RECALIBRATION

	cd $BASE_DIR/misc/nrmap
	echo "recalibrating BAM file...."
	./nrmap recalbam --bamfilein $BASE_DIR/varsim_run/alignments/$PREFIX.$f.bam --bamfileout $BASE_DIR/varsim_run/alignments/tmp.bam --qualsfile $BASE_DIR/models/$PREFIX.scores.tab
	samtools sort -o $BASE_DIR/varsim_run/alignments/$PREFIX.$f.recal.bam -O bam -@ 12 $BASE_DIR/varsim_run/alignments/tmp.bam
	samtools index $BASE_DIR/varsim_run/alignments/$PREFIX.$f.recal.bam
	rm /home/eliot/tmp.bam


	# SNP CALLING WITH BCFTOOLS

	echo "calling SNPs...."


	cd $BASE_DIR/varsim_run/variant_calls

	./call_snps.sh $BASE_DIR/varsim_run/alignments/$prefix.$f.bam $BASE_DIR/varsim_run/reference/$CUCFASTA nocal.vcf
	./call_snps.sh $BASE_DIR/varsim_run/alignments/$prefix.$f.recal.bam $BASE_DIR/varsim_run/reference/$CUCFASTA recal.vcf

	mv nocal.vcf $BASE_DIR/varsim_run/variant_calls
	mv recal.vcf $BASE_DIR/varsim_run/variant_calls


	# SNP CALLING WITH FREEBAYES

	#cd $BASE_DIR/misc_src/freebayes/scripts

	# ./freebayes-parallel <(./fasta_generate_regions.py $BASE_DIR/varsim_run/reference/$fasta.fai 1000000) 16 --use-mapping-quality \
    # -f $BASE_DIR/varsim_run/reference/$fasta $BASE_DIR/varsim_run/alignments/$prefix.$f.bam > $BASE_DIR/varsim_run/variant_calls/cuc_subsampled/CUC_nocal.vcf
	# ./freebayes-parallel <(./fasta_generate_regions.py $BASE_DIR/varsim_run/reference/$fasta.fai 1000000) 16 --use-mapping-quality \
    # -f $BASE_DIR/varsim_run/reference/$fasta $BASE_DIR/varsim_run/alignments/$prefix.$f.recal.bam > $BASE_DIR/varsim_run/variant_calls/cuc_subsampled/CUC_recal.vcf

	# VCF COMPARE

	echo "comparing VCF files...."
	cd $BASE_DIR/varsim_run/variant_calls
	java -jar ../VarSim.jar vcfcompare -true_vcf $BASE_DIR/varsim_run/input_vcfs/cuc_implant.vcf -prefix recal recal.vcf
	java -jar ../VarSim.jar vcfcompare -true_vcf $BASE_DIR/varsim_run/input_vcfs/cuc_implant.vcf -prefix nocal nocal.vcf

	sed -i '/\t\t/d' recal_FP.vcf
	sed -i '/\t\t/d' nocal_FP.vcf

	sed "s/^1/chr1/" -i recal_TP.vcf
	sed "s/^2/chr2/" -i recal_TP.vcf
	sed "s/^3/chr3/" -i recal_TP.vcf
	sed "s/^4/chr4/" -i recal_TP.vcf
	sed "s/^5/chr5/" -i recal_TP.vcf
	sed "s/^6/chr6/" -i recal_TP.vcf
	sed "s/^7/chr7/" -i recal_TP.vcf
	sed "s/^1/chr1/" -i recal_FP.vcf
	sed "s/^2/chr2/" -i recal_FP.vcf
	sed "s/^3/chr3/" -i recal_FP.vcf
	sed "s/^4/chr4/" -i recal_FP.vcf
	sed "s/^5/chr5/" -i recal_FP.vcf
	sed "s/^6/chr6/" -i recal_FP.vcf
	sed "s/^7/chr7/" -i recal_FP.vcf

	sed "s/^1/chr1/" -i nocal_TP.vcf
	sed "s/^2/chr2/" -i nocal_TP.vcf
	sed "s/^3/chr3/" -i nocal_TP.vcf
	sed "s/^4/chr4/" -i nocal_TP.vcf
	sed "s/^5/chr5/" -i nocal_TP.vcf
	sed "s/^6/chr6/" -i nocal_TP.vcf
	sed "s/^7/chr7/" -i nocal_TP.vcf
	sed "s/^1/chr1/" -i nocal_FP.vcf
	sed "s/^2/chr2/" -i nocal_FP.vcf
	sed "s/^3/chr3/" -i nocal_FP.vcf
	sed "s/^4/chr4/" -i nocal_FP.vcf
	sed "s/^5/chr5/" -i nocal_FP.vcf
	sed "s/^6/chr6/" -i nocal_FP.vcf
	sed "s/^7/chr7/" -i nocal_FP.vcf

	cd $BASE_DIR/PythonProgs

	./snp_stats.py --base-dir $BASE_DIR

done

