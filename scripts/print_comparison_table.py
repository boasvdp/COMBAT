#!/usr/bin/env python3

import pandas as pd
import sys

path_tbl_ATCC = str(sys.argv[1])
path_tbl_K12 = str(sys.argv[2])
path_tbl_DSM = str(sys.argv[3])

tbl_ATCC = pd.read_csv(path_tbl_ATCC, sep = '\t')
tbl_K12 = pd.read_csv(path_tbl_K12, sep = '\t')
tbl_DSM = pd.read_csv(path_tbl_DSM, sep = '\t')

print("SNPs/Mbp threshold used", "Persistent travellers (ref: ATCC25922)", "Persistent travellers (ref: K12)", "Persistent travellers (ref: DSM30083)", "Travellers sharing strains (ref: ATCC25922)", "Travellers sharing strains (ref: K12)", "Travellers sharing strains (ref: DSM30083)", "F1 score (ref: ATCC25922)", "F1 score (ref: K12)", "F1 score (ref: DSM30083)", "Precision (ref: ATCC25922)", "Precision (ref: K12)", "Precision (ref: DSM30083)", "Recall (ref: ATCC25922)", "Recall (ref: K12)", "Recall (ref: DSM30083)", sep = '\t')

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
	
	ATCC_PRECISION = ATCC_TP / (ATCC_TP + ATCC_FP)
	ATCC_RECALL = ATCC_TP / 19
	ATCC_FSCORE_2 = 2 * ((ATCC_PRECISION * ATCC_RECALL) / (ATCC_PRECISION + ATCC_RECALL))
	K12_PRECISION = K12_TP / (K12_TP + K12_FP)
	K12_RECALL = K12_TP / 19
	K12_FSCORE_2 = 2 * ((K12_PRECISION * K12_RECALL) / (K12_PRECISION + K12_RECALL))
	DSM_PRECISION = DSM_TP / (DSM_TP + DSM_FP)
	DSM_RECALL = DSM_TP / 19
	DSM_FSCORE_2 = 2 * ((DSM_PRECISION * DSM_RECALL) / (DSM_PRECISION + DSM_RECALL))
	
	print(snp_threshold, ATCC_TP, K12_TP, DSM_TP, ATCC_FP, K12_FP, DSM_FP, "%.3f" % ATCC_FSCORE, "%.3f" % K12_FSCORE, "%.3f" % DSM_FSCORE, "%.3f" % ATCC_PRECISION, "%.3f" % K12_PRECISION, "%.3f" % DSM_PRECISION, "%.3f" % ATCC_RECALL, "%.3f" % K12_RECALL, "%.3f" % DSM_RECALL, sep = '\t')
