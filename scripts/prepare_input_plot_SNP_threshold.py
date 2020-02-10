#!/usr/bin/env python3

import pandas as pd

tbl = pd.read_csv('snp_comparison/snp_comparisons.tsv', sep = '\t')

print("Method", "SNP_threshold", "Comparison", "Number_isolate_pairs", sep = '\t')

for snp_threshold in range(1,21):
	for comparison in [ 'between_traveler', 'within_traveler_within_timepoint' ]:
		isolate_pairs = tbl.query('comparison == @comparison & SNPs_corrected < @snp_threshold').shape[0]
		print("MPA", snp_threshold, comparison, isolate_pairs, sep = '\t')

for snp_threshold in range(1,21):
	for comparison in [ 'between_traveler', 'within_traveler_within_timepoint' ]:
		isolate_pairs = tbl.query('comparison == @comparison & SNPs_not_corrected < @snp_threshold').shape[0]
		print("MPA_nocorr", snp_threshold, comparison, isolate_pairs, sep = '\t')

for snp_threshold in range(1,21):
	for comparison in [ 'between_traveler', 'within_traveler_within_timepoint' ]:
		isolate_pairs = tbl.query('comparison == @comparison & SNPs_no_gaps < @snp_threshold').shape[0]
		print("Core_genome_nogaps", snp_threshold, comparison, isolate_pairs, sep = '\t')
