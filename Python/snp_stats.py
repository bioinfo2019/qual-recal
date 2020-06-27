#!/path/to/anaconda3/bin/python

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 10:54:33 2017

@author: eliot
"""

import argparse
import vcf

################## Parse command line args #############################
    
parser = argparse.ArgumentParser(description='Train model and write it to disk')

parser.add_argument('-b','--base-dir', help='Top level folder containing scripts and data', required=True)
args = parser.parse_args()

top_level = args.base_dir

########################################################################

BASE_DIR = top_level + '/varsim_run/variant_calls/'
before_fn = 'nocal_FN.vcf'
after_fn = 'recal_FN.vcf'
before_vcf = 'nocal.vcf'
after_vcf = 'recal.vcf'
before_fp = 'nocal_FP.vcf'
before_tp = 'nocal_TP.vcf'
after_fp = 'recal_FP.vcf'
after_tp = 'recal_TP.vcf'

quality_threshold = 20

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

before_fn_calls = 0
after_fn_calls = 0

vcf_reader = vcf.Reader(open(BASE_DIR + before_fn, 'r'))
for record in vcf_reader:
    if record.is_snp:
        before_fn_calls += 1

vcf_reader = vcf.Reader(open(BASE_DIR + after_fn, 'r'))
for record in vcf_reader:
    if record.is_snp:
        after_fn_calls += 1
        
vcf_reader = vcf.Reader(open(BASE_DIR + before_vcf, 'r'))
for record in vcf_reader:
    if record.is_snp:
        key = record.CHROM + str(record.POS)
        norecal_quals[key] = record.QUAL

vcf_reader = vcf.Reader(open(BASE_DIR + before_fp, 'r'))
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in norecal_quals:
        norecal_fp_quals.append(norecal_quals[key])
        nr_total_fp += 1
        if norecal_quals[key] >= quality_threshold:
            nr_fp_gt_20 += 1
 
vcf_reader = vcf.Reader(open(BASE_DIR + before_tp, 'r'))
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in norecal_quals:
        norecal_tp_quals.append(norecal_quals[key])
        nr_total_tp += 1
        if norecal_quals[key] >= quality_threshold:
            nr_tp_gt_20 += 1

foo = 0
bar = 0
vcf_reader = vcf.Reader(open(BASE_DIR + after_vcf, 'r'))
for record in vcf_reader:
    if record.is_snp:
        key = record.CHROM + str(record.POS)
        recal_quals[key] = record.QUAL
        if record.QUAL >= 20:
            foo += 1
        bar += 1  

vcf_reader = vcf.Reader(open(BASE_DIR + after_fp, 'r'))
#vcf_writer = vcf.Writer(open(BASE_DIR + '/varsim_run/corrected_FP_gt_20.vcf', 'w'), vcf_reader)
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in recal_quals:
        recal_fp_quals.append(recal_quals[key])
        r_total_fp += 1
        if recal_quals[key] >= quality_threshold:
            #if key not in norecal_quals or (key in norecal_quals and norecal_quals[key] <= 20):
                #vcf_writer.write_record(record)
            r_fp_gt_20 += 1
            
vcf_reader = vcf.Reader(open(BASE_DIR + after_tp, 'r'))
for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    if key in recal_quals:
        recal_tp_quals.append(recal_quals[key])
        r_total_tp += 1
        if recal_quals[key] >= quality_threshold:
            r_tp_gt_20 += 1

before_precision = float(nr_tp_gt_20)/(float(nr_tp_gt_20) + float(nr_fp_gt_20))
after_precision = float(r_tp_gt_20)/(float(r_tp_gt_20) + float(r_fp_gt_20))

# Recall is not actually correct here. VarSim seems to count False Negatives differently for the before and after.
# For the paper, I just recalculated recall using the total number of SNPs implanted into the simulated genomes

before_recall = float(nr_tp_gt_20)/(float(nr_tp_gt_20) + float(before_fn_calls))
after_recall = float(r_tp_gt_20)/(float(r_tp_gt_20) + float(after_fn_calls))

print("No Recal:")
print("TP > 20: " + str(nr_tp_gt_20))
print("FP > 20: " + str(nr_fp_gt_20))
print("Total TP: "+ str(nr_total_tp))
print("Total FP: "+ str(nr_total_fp))
print("Precision: " + str(before_precision))
print("Recall: " + str(before_recall))
print("")
print("Recal:")
print("TP > 20: " + str(r_tp_gt_20))
print("FP > 20: " + str(r_fp_gt_20))
print("Total TP: "+ str(r_total_tp))
print("Total FP: "+ str(r_total_fp))
print("Precision: " + str(after_precision))
print("Recall: " + str(after_recall))

prec_pct_change = (float(after_precision) / float(before_precision) - 1) * 100.0
recall_pct_change = (float(after_recall) / float(before_recall) - 1) * 100.0

#this just appends each line to a file called stats.tab. primitive, I know ....
with open(top_level + "/models/stats.tab", "a") as f:
    
    l = [str(nr_tp_gt_20), str(nr_fp_gt_20), str(nr_total_tp), str(nr_total_fp), str(before_fn_calls), str(r_tp_gt_20), str(r_fp_gt_20), str(r_total_tp), str(r_total_fp), str(after_fn_calls)]

    l.append(str(before_precision))
    l.append(str(before_recall))
    l.append(str(after_precision))
    l.append(str(after_recall))
    l.append(str(prec_pct_change))
    l.append(str(recall_pct_change))

    f.write('\t'.join(l) + '\n')
