#!/usr/bin/env python3

import pandas as pd
import numpy as np
import itertools

# open list_strains.txt as list_strains
list_strains = open("list_strains.txt", "r").read().splitlines()

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

snpmat_standard = pd.read_csv('snp_comparison/snp_distances_standard.tsv', sep = '\t', index_col = 0)
snpmat_standard = snpmat_standard.stack().reset_index()
snpmat_standard.columns = ['strain1','strain2','SNPs_not_corrected']

snpmat_nogaps = pd.read_csv('snp_comparison/snp_distances_no_gaps.tsv', sep = '\t', index_col = 0)
snpmat_nogaps = snpmat_nogaps.stack().reset_index()
snpmat_nogaps.columns = ['strain1','strain2','SNPs_nogaps']

snp_comparisons = pd.merge(df, snpmat_nogaps, how='left', left_on=['strain1','strain2'], right_on = ['strain1','strain2'])
snp_comparisons = pd.merge(snp_comparisons, snpmat_standard, how='left', left_on=['strain1','strain2'], right_on = ['strain1','strain2'])
snp_comparisons.to_csv('snp_comparisons_intermediate.tsv', sep ='\t', index=False)
