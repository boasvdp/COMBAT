#!/usr/bin/env python3

import pandas as pd
import argparse

#Parse argument
parser = argparse.ArgumentParser(description='Combine MLST and SNP typing data to print traveler outcome.')

parser.add_argument('-i', '--isolate-comparisons', dest='isol', help='Isolate comparisons file', required=True, type=str)
parser.add_argument('-t', '--threshold', dest='threshold', help='SNP threshold to use', required=True, type=int)
parser.add_argument('-e', '--exclude-st', dest='exclude', nargs='+', help='ExPEC sequence types to ignore', required=True, type=str)

args = parser.parse_args()

df = pd.read_csv(args.isol, sep = '\t')
threshold = args.threshold

print('traveler', 'MLST_verdict', 'MLST', 'SNP_verdict', 'least_SNPs_masked', 'isolate_T1', 'isolate_T5', sep = '\t')

for traveler in df['traveler'].unique():
  df_trav = df.query('traveler == @traveler')

  common_ST = set(df_trav['ST_T1']) & set(df_trav['ST_T5'])
  if len(common_ST) > 0:
    MLST_verdict = 'Yes'
    MLST = str(common_ST).strip('{}\'')
    df_trav_mlst = df_trav.query('ST_T1 == @MLST & ST_T5 == @MLST')
    least_SNPs_masked = min(df_trav_mlst['SNPs_masked_scaled'])
    if (least_SNPs_masked <= threshold) & (MLST not in args.exclude):
      isolate_T1 = df_trav_mlst.query('SNPs_masked_scaled == @least_SNPs_masked')[['isolate_T1']].iloc[0,0]
      isolate_T5 = df_trav_mlst.query('SNPs_masked_scaled == @least_SNPs_masked')[['isolate_T5']].iloc[0,0]
      SNP_verdict = 'Yes'
    else:
      SNP_verdict = 'No'
      isolate_T1 = 'NA'
      isolate_T5 = 'NA'

  else:
    MLST_verdict = 'No'
    MLST = 'NA'
    SNP_verdict = 'No'
    isolate_T1 = 'NA'
    isolate_T5 = 'NA'
    least_SNPs_masked = 'NA'

  print(traveler, MLST_verdict, MLST, SNP_verdict, least_SNPs_masked, isolate_T1, isolate_T5, sep = '\t')
