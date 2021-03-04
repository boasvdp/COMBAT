#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import stats
import sys
import argparse

# Parse argument
parser = argparse.ArgumentParser(description='Find within and between cluster comparisons based on SNP distances and fastBAPS clustering.')

parser.add_argument('--input', help="Input traveler conclusions", dest='input', required=True, type=str)
parser.add_argument('--expec', help="List of STs that represent ExPEC lineages", dest='expec', required=True, type=str)
parser.add_argument('--controls', help="Control metadata", dest='controls', required=True, type=str)
parser.add_argument('--mlst-controls', help="MLST typing of controls", dest='mlstcontrols', required=True, type=str)

args = parser.parse_args()

path_input = 'travelers_conclusions.tsv' # args.input
path_expec = 'Manges_ExPEC_lineages.tsv' # args.expec
path_controls = 'COMBAT_metadata_controls.tsv' # args.controls
path_mlst_controls = 'mlst_summary_controls.tsv' # args.mlstcontrols

df = pd.read_csv(path_input, sep = '\t')
df_controls = pd.read_csv(path_controls, sep ='\t')
df_mlst_controls = pd.read_csv(path_mlst_controls, sep ='\t', names = ['strain', 'MLST'], header=None)

e = open(path_expec, 'r+')
expec = e.readlines()
e.close()

expec = [x.strip('\n') for x in expec]

df_snps = df.query('SNP_verdict == "Yes"')

df_selected_controls = pd.DataFrame()

for traveler in df_snps['traveler']:
  df_tmp = df_controls.query('matched_to == @traveler')
  tmplist = [df_selected_controls, df_tmp]
  df_selected_controls = pd.concat(tmplist)

df_selected_controls_mlst = pd.merge(df_selected_controls, df_mlst_controls, how='left', on='strain')

df_snps = df_snps.drop(columns = ['MLST_verdict', 'SNP_verdict', 'least_SNPs_masked', 'isolate_T1', 'isolate_T5'])
df_selected_controls_mlst = df_selected_controls_mlst.drop(columns = ['strain', 'timepoint'])

df_snps = df_snps.assign(matched_to='NA')
df_snps = df_snps.assign(carriage='longterm')
df_selected_controls_mlst = df_selected_controls_mlst.assign(carriage='shortterm')

df_total = pd.concat([df_snps, df_selected_controls_mlst], sort=False, ignore_index=True)

df_total["ExPEC_coding"] = np.where(df_total["MLST"].isin(expec), 1, 0)

df_max = df_total.groupby('traveler')['carriage', 'ExPEC_coding'].max()

df_max['ExPEC_status'] = df_max['ExPEC_coding'].map({1: 'ExPEC', 0: 'Non-ExPEC'})



output = pd.crosstab(df_max.carriage, df_max.ExPEC_status)

output.to_csv(sys.stdout, sep = '\t')

print('P-value of Fisher\'s exact test is ' + str(stats.fisher_exact(np.array(output))[1]))
