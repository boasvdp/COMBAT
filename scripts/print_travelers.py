#!/usr/bin/env python3

import sys
import pandas as pd

T1 = float(sys.argv[1])
T2 = float(sys.argv[2])

snp_comparisons = pd.read_csv("snp_comparison/snp_comparisons.tsv", sep = '\t')

snp_comparisons = snp_comparisons.query('comparison == "within_traveler_between_timepoint"')

print("Traveler", "SNPs_corrected", "Type", sep = '\t')

for traveler in list(snp_comparisons.traveler1.drop_duplicates()):
	min_SNPs = min(snp_comparisons.query('traveler1 == @traveler')['SNPs_corrected'])
	if min_SNPs < T1:
		type = "Very_likely_clonal"
	elif min_SNPs >= T1 and min_SNPs < T2:
		type = "Likely_clonal"
	else:
		type = "Not_clonal"
	print(traveler, "%.3f" % min_SNPs, type, sep = '\t')
