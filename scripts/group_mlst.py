#!/usr/bin/env python3

import pandas as pd
import os
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Group strains per ST and traveller and print to table.')

parser.add_argument('metadata', help="Input metadata file, containing info which isolate originates from which traveler at which timepoint", type=str)
parser.add_argument('mlst', help="MLST summary file of all isolates", type=str)
parser.add_argument('outdir', help="Output directory", default='tables_per_ST', type=str)

args = parser.parse_args()

# Read in metadat and ST data
df = pd.read_csv(args.metadata, sep = '\t', header=None)
mlst = pd.read_csv(args.mlst, sep = '\t')

# Rename metadata columns
df.columns = ['traveler', 'isolate', 'timepoint']

# Merge metadata and ST columns
df_mlst = pd.merge(df, mlst, left_on = 'isolate', right_on = 'strain')

# Make output dir if it doesn't exist already
if not os.path.exists(args.outdir):
  os.makedirs(args.outdir)

# Loop through STs present
for st in df_mlst['st'].unique():
  # Make temp df with only that ST
  df_mlst_st = df_mlst.query('st == @st')
  # Create empty output df and string for output table for this ST
  out_df = pd.DataFrame(columns = ['isolate', 'st', 'traveler', 'timepoint'])
  outstring = args.outdir + '/' + 'table_mapping_ST' + str(st) + '.tsv'

  # Loop through travelers in ST df and subset a df
  for traveler in df_mlst_st['traveler'].unique():    
    df_mlst_st_trav = df_mlst_st.query('traveler == @traveler')
    # Check whether two timepoints are there
    if len(df_mlst_st_trav['timepoint'].unique()) == 2:
      # Concatenate empty output dataframe and temp df containing data
      tmpdf = df_mlst_st_trav[['isolate', 'st', 'traveler', 'timepoint']]
      list_df = [out_df, tmpdf]
      out_df = pd.concat(list_df)
  # Write df to file if it contains data
  if out_df.shape[0] > 0:
    out_df.to_csv(outstring, sep = '\t', index=False)
