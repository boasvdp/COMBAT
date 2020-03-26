#!/usr/bin/env python3

import pandas as pd
import numpy as np
import itertools
import sys

path_list_strains = str(sys.argv[1])
path_snpmat_standard = str(sys.argv[2])
path_snpmat_nogaps = str(sys.argv[3])
path_alnlengths = str(sys.argv[4])
path_output = str(sys.argv[5])

# open list_strains.txt as list_strains
list_strains = open(path_list_strains, "r").read().splitlines()

# use itertools to make all unique combinations of list_strains with itself, convert to dataframe and give column names
combos = pd.DataFrame(list(itertools.combinations(list_strains,2)), columns = ['strain1', 'strain2'])

# read metadata file and give NO column names
metadata = pd.read_csv("COMBAT_metadata.tsv", sep = '\t', header = None)

# match metadata values to the combinations in the first column of combos (strain1)
metadata.columns = ['traveler', 'strain1', 'timepoint']
df = pd.merge(combos, metadata, on = 'strain1')
df = df.rename(columns={'traveler': 'traveler1', 'timepoint': 'timepoint1'})

# Change metadata name for matching, and repeat above but for strain2
metadata.columns = ['traveler', 'strain2', 'timepoint']
df = pd.merge(df, metadata, on = 'strain2')
df = df.rename(columns={'traveler': 'traveler2', 'timepoint': 'timepoint2'})

# Define three conditions
conditions = [
	(df['traveler1'] == df['traveler2']) & (df['timepoint1'] == df['timepoint2']),
	(df['traveler1'] == df['traveler2']) & (df['timepoint1'] != df['timepoint2']),
	df['traveler1'] != df['traveler2']]

# Define values per condition
choices = [
	'within_traveler_within_timepoint',
	'within_traveler_between_timepoint',
	'between_traveler']

# Assign comparison types using numpy and previously defined conditions and choices
df['comparison'] = np.select(conditions, choices, default=np.nan)

# Read the files with molten SNPs (standard and no gaps) and alignment lengths
snpmat_standard = pd.read_csv(path_snpmat_standard, sep = '\t', names = ['strain1','strain2','SNPs_not_corrected'])
snpmat_nogaps = pd.read_csv(path_snpmat_nogaps, sep = '\t', names = ['strain1','strain2','SNPs_no_gaps'])
alnlengths = pd.read_csv(path_alnlengths, sep = '\t', names = ['strain1','strain2','alignment_length'])

# Merge the SNP and aln lengths data with the traveler/strain metadata
snp_comparisons = pd.merge(df, alnlengths, how='left', left_on=['strain1','strain2'], right_on = ['strain1','strain2'])
snp_comparisons = pd.merge(snp_comparisons, snpmat_nogaps, how='left', left_on=['strain1','strain2'], right_on = ['strain1','strain2'])
snp_comparisons = pd.merge(snp_comparisons, snpmat_standard, how='left', left_on=['strain1','strain2'], right_on = ['strain1','strain2'])

# Calculate SNPs per mbp
snp_comparisons['SNPs_corrected'] = (snp_comparisons['SNPs_not_corrected'] / snp_comparisons['alignment_length']) * 1000000

# Write to output tsv file
snp_comparisons.to_csv(path_output, sep ='\t', index=False)
