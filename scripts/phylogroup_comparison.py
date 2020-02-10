#!/usr/bin/env python3

import pandas as pd
import sys

### For long-term colonising strains
# Read snp_comparisons file, which will be used to find persistent strains which can be linked to phylogroups
snp_comparisons = pd.read_csv('snp_comparison/snp_comparisons.tsv', sep = '\t')

# Read tsv file which summarises which strain is which phylogroup
ezclermont = pd.read_csv('ezclermont_summary.tsv', sep  = '\t', names = ["strain", "phylogroup"])

# Set threshold for clonal persistence from argument passed from the command line. Defined in config.yaml.
T = int(sys.argv[1])

# Find strains which have persisted based on threshold defined above
snp_comparisons_clonal = snp_comparisons.query('comparison  == "within_traveler_between_timepoint" & SNPs_corrected <= @T')

# Merge dataframes: ezclermont phylogroups are added to persisting strains
snp_comparisons_clonal = pd.merge(left = snp_comparisons_clonal, right = ezclermont, left_on = 'strain2', right_on = 'strain')

# Take only the traveler and phylogroup columns, and remove duplicates to prevent inflated phylogroup numbers
snp_comparisons_clonal = snp_comparisons_clonal[['traveler1', 'phylogroup']].drop_duplicates()

### For short-term colonising strains
# Read the metadata file for controls
metadata_controls = pd.read_csv('COMBAT_metadata_controls.tsv', sep = '\t')

# Read the phylogroup summary for controls
ezclermont_controls = pd.read_csv('ezclermont_summary_controls.tsv', sep  = '\t', names = ["strain", "phylogroup"])

# For travelers in the metadata dataframe, get corresponding phylogroups
ezclermont_controls_travelers = pd.merge(left = ezclermont_controls, right = metadata_controls, on = 'strain')

# Take only the phylogroup counts from both dataframes, and add a column which specifies long/short term
phylogroups_longterm = pd.DataFrame(snp_comparisons_clonal["phylogroup"]).assign(type='longterm')
phylogroups_shortterm = pd.DataFrame(ezclermont_controls_travelers["phylogroup"]).assign(type='shortterm')

# Append dataframes
phylogroups = phylogroups_longterm.append(phylogroups_shortterm)

# Group the resulting dataframe by long/short term, calculate how often phylogroups occur, and format the table with unstack, fillna and transpose
out = phylogroups.groupby('type')['phylogroup'].value_counts().unstack().fillna(0).transpose()

# Write to a tsv file, and switch float to integer by formatting
out.to_csv('phylogroup_comparison.tsv', sep = '\t', float_format='%.f')
