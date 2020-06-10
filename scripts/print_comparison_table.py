#!/usr/bin/env python3

import pandas as pd
import sys

path_tbl_ATCC = str(sys.argv[1])
path_tbl_K12 = str(sys.argv[2])
path_tbl_DSM = str(sys.argv[3])

tbl_ATCC = pd.read_csv(path_tbl_ATCC, sep = '\t')
tbl_K12 = pd.read_csv(path_tbl_K12, sep = '\t')
tbl_DSM = pd.read_csv(path_tbl_DSM, sep = '\t')

print("SNPs/Mbp threshold used", "Persistent travellers (ref: ATCC25922)", "Persistent travellers (ref: K12)", "Persistent travellers (ref: DSM30083)", "Travellers sharing strains (ref: ATCC25922)", "Travellers sharing strains (ref: K12)", "Travellers sharing strains (ref: DSM30083)", "F1 score (ref: ATCC25922)", "F1 score (ref: K12)", "F1 score (ref: DSM30083)", sep = '\t')

for snp_threshold in range(1,1001):
	h = tbl_ATCC.query('comparison == "between_traveler"  & SNPs_corrected < @snp_threshold')[['traveler1', 'traveler2']].drop_duplicates()
	df = []
	for ind in h.index:
		list = [h['traveler1'][ind], h['traveler2'][ind]]
		list.sort()
		df.append(list)
	ATCC_FP = pd.DataFrame(df).drop_duplicates().shape[0]
	ATCC_TP = tbl_ATCC.query('comparison == "within_traveler_between_timepoint"  & SNPs_corrected < @snp_threshold').traveler1.drop_duplicates().shape[0]
	h = tbl_K12.query('comparison == "between_traveler"  & SNPs_corrected < @snp_threshold')[['traveler1', 'traveler2']].drop_duplicates()
	df = []
	for ind in h.index:
		list = [h['traveler1'][ind], h['traveler2'][ind]]
		list.sort()
		df.append(list)
	K12_FP = pd.DataFrame(df).drop_duplicates().shape[0]
	K12_TP = tbl_K12.query('comparison == "within_traveler_between_timepoint"  & SNPs_corrected < @snp_threshold').traveler1.drop_duplicates().shape[0]
	h = tbl_DSM.query('comparison == "between_traveler"  & SNPs_corrected < @snp_threshold')[['traveler1', 'traveler2']].drop_duplicates()
	df = []
	for ind in h.index:
		list = [h['traveler1'][ind], h['traveler2'][ind]]
		list.sort()
		df.append(list)
	DSM_FP = pd.DataFrame(df).drop_duplicates().shape[0]
	DSM_TP = tbl_DSM.query('comparison == "within_traveler_between_timepoint"  & SNPs_corrected < @snp_threshold').traveler1.drop_duplicates().shape[0]
	if ATCC_TP == 0:
		ATCC_FSCORE = 0
	else:
		ATCC_FSCORE = 1 / (((1/(ATCC_TP / (ATCC_TP + ATCC_FP))) + (1/(ATCC_TP / 19))) / 2)
	if K12_TP == 0:
		K12_FSCORE = 0
	else:
		K12_FSCORE = 1 / (((1/(K12_TP / (K12_TP + K12_FP))) + (1/(K12_TP / 19))) / 2)
	if DSM_TP == 0:
		DSM_FSCORE = 0
	else:
		DSM_FSCORE = 1 / (((1/(DSM_TP / (DSM_TP + DSM_FP))) + (1/(DSM_TP / 19))) / 2)

	print(snp_threshold, ATCC_TP, K12_TP, DSM_TP, ATCC_FP, K12_FP, DSM_FP, "%.3f" % ATCC_FSCORE, "%.3f" % K12_FSCORE, "%.3f" % DSM_FSCORE, sep = '\t')
