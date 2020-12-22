#!/usr/bin/env python3

import pandas as pd
import argparse

# Parse argument
parser = argparse.ArgumentParser(description='Combine SNP typing tables.')

parser.add_argument('-l', '--list', nargs='+', dest='l', help='Tables containing SNP comparisons', required=True, type=str)
parser.add_argument('-o', '--output', dest='o', help='Output table', required=True, type=str)

args = parser.parse_args()

df_out = pd.DataFrame()

for i in args.l:
  tmpdf = pd.read_csv(i, sep = '\t')
  list_df = [df_out, tmpdf]
  df_out = pd.concat(list_df)

df_out.to_csv(args.o, sep = '\t', index=False)
