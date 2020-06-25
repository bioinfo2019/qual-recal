#!/bin/bash
region=$1
bamfile=$2
ref=$3

bcftools mpileup -A -C 0 -Q 5 --ignore-overlaps --no-BAQ -r $region -f $ref $bamfile | bcftools call --variants-only --multiallelic-caller -f GQ -O b - > $region.bcf
