#!/usr/local/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 10:54:33 2017

@author: eliot
"""

import sys
import numpy as np
from scipy import stats
import vcf
import matplotlib.pyplot as plt

#def main(argv):

BASE_DIR = '/hdd1/PhD'
recal_tp_quals = []
recal_fp_quals = []

norecal_tp_quals = []
norecal_fp_quals = []

norecal_quals = {}
recal_quals = {}
r_fp_gt_20 = 0
r_tp_gt_20 = 0
nr_fp_gt_20 = 0
nr_tp_gt_20 = 0
r_total_fp = 0
r_total_tp = 0
nr_total_fp = 0
nr_total_tp = 0
total_snps = 0



vcf_reader = vcf.Reader(open(BASE_DIR + '/varsim_run/out/naive25xsnps.truth.vcf', 'r'))
for record in vcf_reader:
    if record.is_snp:
        total_snps += 1
        
        
vcf_reader = vcf.Reader(open(BASE_DIR + '/overlapper/calling_results/unmerged.vcf', 'r'))
for record in vcf_reader:
    if record.is_snp:
        key = record.CHROM + str(record.POS)
        norecal_quals[key] = record.QUAL


vcf_reader = vcf.Reader(open(BASE_DIR + '/overlapper/calling_results/unmerged_FP.vcf', 'r'))
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in norecal_quals:
        norecal_fp_quals.append(norecal_quals[key])
        nr_total_fp += 1
        if norecal_quals[key] > 20:
            nr_fp_gt_20 += 1
 
vcf_reader = vcf.Reader(open(BASE_DIR + '/overlapper/calling_results/unmerged_TP.vcf', 'r'))
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in norecal_quals:
        norecal_tp_quals.append(norecal_quals[key])
        nr_total_tp += 1
        if norecal_quals[key] > 20:
            nr_tp_gt_20 += 1

vcf_reader = vcf.Reader(open(BASE_DIR + '/overlapper/calling_results/merged.vcf', 'r'))
for record in vcf_reader:
    if record.is_snp:
        key = record.CHROM + str(record.POS)
        recal_quals[key] = record.QUAL
        
vcf_reader = vcf.Reader(open(BASE_DIR + '/overlapper/calling_results/merged_FP.vcf', 'r'))
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in recal_quals:
        recal_fp_quals.append(recal_quals[key])
        r_total_fp += 1
        if recal_quals[key] > 20:
            r_fp_gt_20 += 1
            
vcf_reader = vcf.Reader(open(BASE_DIR + '/overlapper/calling_results/merged_TP.vcf', 'r'))
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in recal_quals:
        recal_tp_quals.append(recal_quals[key])
        r_total_tp += 1
        if recal_quals[key] > 20:
            r_tp_gt_20 += 1

print "Total SNPS: " + str(total_snps)
print "No Recal:"
print "TP > 20: " + str(nr_tp_gt_20)
print "FP > 20: " + str(nr_fp_gt_20)
print "Total TP: "+ str(nr_total_tp)
print "Total FP: "+ str(nr_total_fp)
print ""
print "Recal:"
print "TP > 20: " + str(r_tp_gt_20)
print "FP > 20: " + str(r_fp_gt_20)
print "Total TP: "+ str(r_total_tp)
print "Total FP: "+ str(r_total_fp)

boxplot_data = [norecal_fp_quals, recal_fp_quals]


# Create a figure instance
fig = plt.figure(1, figsize=(9, 6))

# Create an axes instance
ax = fig.add_subplot(111)



# Create the boxplot
bp = ax.boxplot(boxplot_data)
## Custom x-axis labels
ax.set_xticklabels(['No recal FP', 'Recal FP'])

# Save the figure
fig.savefig(BASE_DIR + '/varsim_run/FP.png', bbox_inches='tight')

plt.show()
#if __name__ == "__main__":
#   main(sys.argv[1:])





























#print "Truepositive: ", stats.ttest_ind(recal_tp_quals, norecal_tp_quals)

#print np.mean(recal_quals)
#print np.mean(norecal_quals)