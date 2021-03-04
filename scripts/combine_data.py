#!/usr/bin/env python3

import pandas as pd
import argparse
import itertools

#Parse argument
parser = argparse.ArgumentParser(description='Combine MLST and SNP typing data to print traveler outcome.')

parser.add_argument('-m', '--metadata', dest='metadata', help='Metadata table', required=True, type=str)
parser.add_argument('-l', '--mlst', dest='mlst', help='MLST table', required=True, type=str)
parser.add_argument('-s', '--snps', dest='snps', help='SNPs all table', required=True, type=str)
parser.add_argument('-t', '--threshold', dest='threshold', help='SNP threshold to use', required=True, type=int)
parser.add_argument('-o', '--main-output', dest='out', help='Main output table', required=True, type=str)
parser.add_argument('-c', '--conclusion-output', dest='dfconclout', help='Conclusions output table', required=True, type=str)

args = parser.parse_args()

# 
path_metadata = args.metadata
path_mlst = args.mlst
path_snps = args.snps
SNP_threshold = args.threshold

df_metadata = pd.read_csv(path_metadata, sep = '\t', names=['traveler', 'isolate', 'timepoint'])
df_mlst = pd.read_csv(path_mlst, sep = '\t', index_col='strain')

mlst_dict = df_mlst['st'].to_dict()

df_combos = pd.DataFrame(columns = ['traveler', 'isolate_T1', 'isolate_T5'])

for traveler in df_metadata['traveler'].unique():
  list_strains_T1 = list(df_metadata.query('traveler == @traveler & timepoint == "T1"')['isolate'])
  list_strains_T5 = list(df_metadata.query('traveler == @traveler & timepoint == "T5"')['isolate'])
  combos_tmp = pd.DataFrame(list(itertools.product(list_strains_T1,list_strains_T5)), columns = ['isolate_T1', 'isolate_T5'])
  combos_tmp.insert(0, 'traveler', traveler)
  df_combos = df_combos.append(combos_tmp, ignore_index=True)

df_combos['ST_T1'] = df_combos['isolate_T1'].map(mlst_dict)
df_combos['ST_T5'] = df_combos['isolate_T5'].map(mlst_dict)

df_snps = pd.read_csv(path_snps, sep = '\t', usecols = ['isolate_T1', 'isolate_T5', 'SNPs_masked', 'lengths_masked', 'SNPs_normal', 'lengths_normal', 'SNPs_masked_scaled', 'SNPs_normal_scaled'])

df_out = df_combos.merge(df_snps, left_on = ['isolate_T1', 'isolate_T5'], right_on = ['isolate_T1', 'isolate_T5'], how = 'left').fillna('NA')

df_out.to_csv(args.out, sep = '\t', index=False)

