#!/usr/bin/env python3

import pandas as pd
import sys

path_snp_comparisons = str(sys.argv[2])
path_ezclermont = str(sys.argv[3])
path_ezclermont_controls = str(sys.argv[4])
path_metadata_controls = str(sys.argv[5])
path_plot_data_out = str(sys.argv[6])
path_output = str(sys.argv[7])

### For long-term colonising strains
# Read snp_comparisons file, which will be used to find persistent strains which can be linked to phylogroups
snp_comparisons = pd.read_csv(path_snp_comparisons, sep = '\t')

# Read tsv file which summarises which strain is which phylogroup
ezclermont = pd.read_csv(path_ezclermont, sep  = '\t', names = ["strain", "phylogroup"])

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
metadata_controls = pd.read_csv(path_metadata_controls, sep = '\t')

# Read the phylogroup summary for controls
ezclermont_controls = pd.read_csv(path_ezclermont_controls, sep  = '\t', names = ["strain", "phylogroup"])

# For travelers in the metadata dataframe, get corresponding phylogroups
metadata_controls = pd.merge(left = metadata_controls, right = ezclermont_controls, left_on = 'strain', right_on = 'strain')

# Initiate empty dataframe to append short-term data to
phylogroups_shortterm = pd.DataFrame(columns = ['phylogroup'])

# Loop over short-term carriers and get phylogroups per traveler. If multiple phylogroups were acquired by a traveler, write 'Multiple'. Then append to empty dataframe.
for t in metadata_controls.traveler.unique():
    phylo = list(metadata_controls[metadata_controls['traveler'] == t].phylogroup.unique())
    if len(phylo) > 1:
        phylo = 'Multiple'
    else:
        phylo = str(phylo).strip('[]\'')
    phylogroups_shortterm = phylogroups_shortterm.append({'phylogroup': phylo}, ignore_index=True)

# Take only the phylogroup counts from both dataframes, and add a column which specifies long/short term
phylogroups_longterm = pd.DataFrame(snp_comparisons_clonal["phylogroup"]).assign(type='longterm')
phylogroups_shortterm = pd.DataFrame(phylogroups_shortterm["phylogroup"]).assign(type='shortterm')

# Append dataframes
phylogroups = phylogroups_longterm.append(phylogroups_shortterm)

# Change EzClermont output 'cryptic' to the more descriptive 'Cryptic clade'
phylogroups = phylogroups.replace('cryptic', 'Cryptic clade')

# Write data to plot to file
plotdata = pd.DataFrame(columns = ['type', 'phylogroup', 'count'])
t = phylogroups.type.unique()
p = phylogroups.phylogroup.unique()

for i in t:
	for j in p:
		s = sum((phylogroups['phylogroup'] == j) & (phylogroups['type'] == i))
		plotdata = plotdata.append({'type': i, 'phylogroup': j, 'count': s}, ignore_index=True)

plotdata.to_csv(path_plot_data_out, sep ='\t', index=False)

# Group the resulting dataframe by long/short term, calculate how often phylogroups occur, and format the table with unstack, fillna and transpose
out = phylogroups.groupby('type')['phylogroup'].value_counts().unstack().fillna(0).transpose()

# Write to a tsv file, and switch float to integer by formatting
out.to_csv(path_output, sep = '\t', float_format='%.f')
