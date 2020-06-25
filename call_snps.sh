#!/bin/bash
bamfile=$1
ref=$2
vcfname=$3
ncpus=8


cut -f 1 $ref.fai | parallel --gnu -j $ncpus ./call_snps_by_contig.sh {} "$bamfile" "$ref"

bcftools concat -O v -o raw.bcf *.bcf

bcftools view raw.bcf > $vcfname

rm -rf *.bcf


