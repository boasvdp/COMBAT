#!/usr/bin/env python3

import pandas as pd
import argparse
import itertools
import numpy as np

# Parse argument
parser = argparse.ArgumentParser(description='Calculated raw and scaled SNPs between isolates.')

parser.add_argument('indir', help="Input data dir", type=str)
parser.add_argument('ST', help="Sequence type to process", type=str)
parser.add_argument('table', help="Mapping table indicating which isolates should be compared", type=str)
parser.add_argument('output', help="Output file", type=str)

args = parser.parse_args()

# Make paths based on input directory
path_df_snps_masked = args.indir + '/snps_ST' + args.ST + '_masked.tsv'
path_df_snps_normal = args.indir + '/snps_ST' + args.ST + '_normal.tsv'
path_df_lengths_masked = args.indir + '/lengths_ST' + args.ST + '_masked.tsv'
path_df_lengths_normal = args.indir + '/lengths_ST' + args.ST + '_normal.tsv'

# Read dataframes with SNP and alignment length numbers
df_snps_masked = pd.read_csv(path_df_snps_masked, sep = '\t', names = ['isolate_T1', 'isolate_T5', 'SNPs_masked'])
df_snps_normal = pd.read_csv(path_df_snps_normal, sep = '\t', names = ['isolate_T1', 'isolate_T5', 'SNPs_normal'])
df_lengths_masked = pd.read_csv(path_df_lengths_masked, sep = '\t', names = ['isolate_T1', 'isolate_T5', 'lengths_masked'])
df_lengths_normal = pd.read_csv(path_df_lengths_normal, sep = '\t', names = ['isolate_T1', 'isolate_T5', 'lengths_normal'])

# Read in metadata from mapping table
df_metadata = pd.read_csv(args.table, sep = '\t')

# Create list of strains
list_strains_T1 = list(df_metadata.query('timepoint == "T1"')['isolate'])
list_strains_T5 = list(df_metadata.query('timepoint == "T5"')['isolate'])

# use itertools to make all unique combinations of list_strains with itself, convert to dataframe and give column names
combos = pd.DataFrame(list(itertools.product(list_strains_T1,list_strains_T5)), columns = ['isolate_T1', 'isolate_T5'])

# match metadata values to the combinations in the first column of combos (strain1)
df_metadata.columns = ['isolate_T1', 'ST', 'traveler', 'timepoint']
df_out = pd.merge(combos, df_metadata, on = 'isolate_T1')
df_out = df_out.rename(columns={'traveler': 'traveler_T1', 'timepoint': 'timepoint_T1'})

# Change metadata name for matching, and repeat above but for strain2
df_metadata = df_metadata.drop(columns = 'ST')
df_metadata.columns = ['isolate_T5', 'traveler', 'timepoint']
df_out = pd.merge(df_out, df_metadata, on = 'isolate_T5')
df_out = df_out.rename(columns={'traveler': 'traveler_T5', 'timepoint': 'timepoint_T5'})

# Define three conditions
conditions = [
	(df_out['traveler_T1'] == df_out['traveler_T5']) & (df_out['timepoint_T1'] == df_out['timepoint_T5']),
	(df_out['traveler_T1'] == df_out['traveler_T5']) & (df_out['timepoint_T1'] != df_out['timepoint_T5']),
	df_out['traveler_T1'] != df_out['traveler_T5']]

# Define values per condition
choices = [
	'within_traveler_within_timepoint',
	'within_traveler_between_timepoint',
	'between_traveler']

# Assign comparison types using numpy and previously defined conditions and choices
df_out['comparison'] = np.select(conditions, choices, default=np.nan)

#df_all_snps = pd.merge(df_snps_masked, df_lengths_masked)
#df_all_snps = pd.merge(df_all_snps, df_snps_normal)
#df_all_snps = pd.merge(df_all_snps, df_lengths_normal)

#df_all_snps_ordered = pd.DataFrame(columns = ['isolate_T1', 'isolate_T5', 'SNPs_masked', 'lengths_masked', 'SNPs_normal', 'lengths_normal'])

#for index, row in df_out.iterrows():
#  ISOL_T1 = row['isolate_T1']
#  ISOL_T5 = row['isolate_T5']
#  df_tmp = df_all_snps.query('(strain1 == @ISOL_T1 & strain2 == @ISOL_T5) | (strain1 == @ISOL_T5 & strain2 == @ISOL_T1)')
#  print('printing dftmp')
#  print(df_tmp)
#  if df_tmp['strain1'] == ISOL_T1:
#    df_all_snps_ordered = df_all_snps_ordered.append({'isolate_T1': strain1, 'isolate_T5': strain2, 'SNPs_masked': df_tmp['SNPs_masked'], 'lengths_masked': df_tmp['lengths_masked'], 'SNPs_normal': df_tmp['SNPs_normal'], 'lengths_normal': df_tmp['lengths_normal']}, ignore_index=True)
#  else:
#    df_all_snps_ordered = df_all_snps_ordered.append({'isolate_T1': strain2, 'isolate_T5': strain1, 'SNPs_masked': df_tmp['SNPs_masked'], 'lengths_masked': df_tmp['lengths_masked'], 'SNPs_normal': df_tmp['SNPs_normal'], 'lengths_normal': df_tmp['lengths_normal']}, ignore_index=True)

# Merge the SNP and aln lengths data with the traveler/strain metadata
df_out = pd.merge(df_out, df_snps_masked, how='left', left_on=['isolate_T1', 'isolate_T5'], right_on = ['isolate_T1', 'isolate_T5'])
df_out = pd.merge(df_out, df_lengths_masked, how='left', left_on=['isolate_T1', 'isolate_T5'], right_on = ['isolate_T1', 'isolate_T5'])
df_out = pd.merge(df_out, df_snps_normal, how='left', left_on=['isolate_T1', 'isolate_T5'], right_on = ['isolate_T1', 'isolate_T5'])
df_out = pd.merge(df_out, df_lengths_normal, how='left', left_on=['isolate_T1', 'isolate_T5'], right_on = ['isolate_T1', 'isolate_T5'])

#df_out = pd.merge(df_out, df_all_snps_ordered, how='left', left_on=['isolate_T1', 'isolate_T5'], right_on = ['isolate_T1', 'isolate_T5'])

# Calculate scaled SNPs, masked and unmasked
df_out['SNPs_masked_scaled'] = (df_out['SNPs_masked'] / df_out['lengths_masked']) * 5000000
df_out['SNPs_normal_scaled'] = (df_out['SNPs_normal'] / df_out['lengths_normal']) * 5000000

# Output final directory
df_out.to_csv(args.output, sep = '\t', index=False)
