#!/usr/bin/env python3

import pandas as pd
import sys
import argparse

# Parse argument
parser = argparse.ArgumentParser(description='Print comparison of phylogroups.')

parser.add_argument('--isolate-comparison', dest='isol', help='Isolate summary', required=True, type=str)
parser.add_argument('--threshold', dest='threshold', help='SNP threshold', required=True, type=int)
parser.add_argument('--ezclermont', dest='ezclermont', help='EzClermont summary', required=True, type=str)
parser.add_argument('--ezclermont-controls', dest='ezclermontcontrols', help='EzClermont summary of controls', required=True, type=str)
parser.add_argument('--metadata-controls', dest='metadatacontrols', help='Metadata of controls', required=True, type=str)
parser.add_argument('--plot-data-out', dest='plotdataout', help='Output file for plotting data', required=True, type=str)
parser.add_argument('--output', dest='output', help='Output file', required=True, type=str)

args = parser.parse_args()

path_ezclermont = args.ezclermont
path_ezclermont_controls = args.ezclermontcontrols
path_metadata_controls = args.metadatacontrols
path_plot_data_out = args.plotdataout
path_output = args.output

### For long-term colonising strains
# Read snp_comparisons file, which will be used to find persistent strains which can be linked to phylogroups
snp_comparisons = pd.read_csv(args.isol, sep = '\t')

# Read tsv file which summarises which strain is which phylogroup
ezclermont = pd.read_csv(path_ezclermont, sep  = '\t', names = ["strain", "phylogroup"])

# Set threshold for clonal persistence from argument passed from the command line. Defined in config.yaml.
T = args.threshold

# Find strains which have persisted based on threshold defined above
snp_comparisons_clonal = snp_comparisons.query('SNPs_masked_scaled != "NA"').query('SNPs_masked_scaled <= @T')

# Merge dataframes: ezclermont phylogroups are added to persisting strains
snp_comparisons_clonal = pd.merge(left = snp_comparisons_clonal, right = ezclermont, left_on = 'isolate_T5', right_on = 'strain')

# Take only the traveler and phylogroup columns, and remove duplicates to prevent inflated phylogroup numbers
snp_comparisons_clonal = snp_comparisons_clonal[['traveler', 'phylogroup']].drop_duplicates()

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
