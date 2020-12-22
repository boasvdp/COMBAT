#!/usr/bin/env python3

import pandas as pd
import sys
import os
import argparse

# Parse argument
parser = argparse.ArgumentParser(description='Find within and between cluster comparisons based on SNP distances and fastBAPS clustering.')

parser.add_argument('indir', help="Directory containing files with selected references", default='selected_references', type=str)
parser.add_argument('mlst', help="MLST summary to which the reference info will be added", default='mlst_summary.tsv', type=str)
parser.add_argument('outfile', help="Output file", default='mlst_references.tsv', type=str)

args = parser.parse_args()

#
dict = {}

for sel_ref in os.listdir(args.indir):
  ST = sel_ref.replace('reference_ST', '').replace('.tsv', '')
  c = open(args.indir + '/' + sel_ref)
  REF = c.readlines()[1].split('\t')[0]
  c.close()
  dict[ST] = REF

df_mlst = pd.read_csv(args.mlst, sep = '\t')

df_mlst['reference'] = df_mlst['st'].map(dict)

df_mlst = df_mlst.dropna()

df_mlst.to_csv(args.outfile, sep = '\t', index=False)
