#!/usr/bin/env python3

import pandas as pd
import sys
import argparse

# Parse argument
parser = argparse.ArgumentParser(description='Select best reference genome for a set of genomes based on ReferenceSeeker output.')

parser.add_argument('table', help="Input table containing genomes to compare", type=str)
parser.add_argument('inputdir', help="Directory containing referenceseeker output tables", default='referenceseeker_out', type=str)

parser.add_argument('output', help="Output file", type=str)

args = parser.parse_args()

df_input_table = pd.read_csv(args.table, sep = '\t')

total_df = pd.DataFrame()

for file in df_input_table['isolate']:
  path_to_file = args.inputdir + '/' + file + '.tsv'
  tmpdf = pd.read_csv(path_to_file, sep = '\t')
  list_dfs = [total_df, tmpdf]
  total_df = pd.concat(list_dfs, ignore_index=True, sort=False)

mash_counts = total_df['#ID'].value_counts()
max_mash_counts = max(mash_counts)

list_mash_hits = list(mash_counts[mash_counts == max_mash_counts].index)

out_df = pd.DataFrame(columns = ['#ID', 'nr_mash_hits', 'mean_conserved_DNA', 'mean_ANI'])

if len(list_mash_hits) == 1:
  BEST_HIT = str(list_cons_dna_hits[0])
  max_cons_dna = total_df[total_df['#ID'] == BEST_HIT]['Cons. DNA']
  max_ani = total_df[total_df['#ID'] == BEST_HIT]['ANI']
  out_df = out_df.append({'#ID': BEST_HIT, 'nr_mash_hits': max_mash_counts, 'mean_conserved_DNA': max_cons_dna, 'mean_ANI': max_ani}, ignore_index=True)
else:
  cons_dna = total_df[total_df['#ID'].isin(list_mash_hits)].groupby('#ID')['Con. DNA'].mean()
  max_cons_dna = max(cons_dna)
  list_cons_dna_hits = list(cons_dna[cons_dna == max_cons_dna].index)
  if len(list_cons_dna_hits) == 1:
    BEST_HIT = str(list_cons_dna_hits[0])
    mean_ani = total_df[total_df['#ID'] == BEST_HIT]['ANI'].mean()
    out_df = out_df.append({'#ID': BEST_HIT, 'nr_mash_hits': max_mash_counts, 'mean_conserved_DNA': max_cons_dna, 'mean_ANI': mean_ani}, ignore_index=True)
  else:
    ani = total_df[total_df['#ID'].isin(list_cons_dna_hits)].groupby('#ID')['ANI'].mean()
    max_ani = max(ani)
    list_ani = list(ani[ani == max_ani].index)
    if len(list_ani) == 1:
      BEST_HIT = str(list_cons_dna_hits[0])
      out_df = out_df.append({'#ID': BEST_HIT, 'nr_mash_hits': max_mash_counts, 'mean_conserved_DNA': max_cons_dna, 'mean_ANI': max_ani}, ignore_index=True)
    else:
      for hit in list_ani:
        BEST_HIT = str(list_cons_dna_hits[0])
        out_df = out_df.append({'#ID': BEST_HIT, 'nr_mash_hits': max_mash_counts, 'mean_conserved_DNA': max_cons_dna, 'mean_ANI': max_ani}, ignore_index=True)

out_df.to_csv(args.output, sep = '\t', index=False)
