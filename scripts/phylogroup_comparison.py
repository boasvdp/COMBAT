#!/usr/bin/env python3

import pandas as pd

snp_comparisons = pd.read_csv('snp_comparisons.tsv', sep = '\t')
metadata_controls = pd.read_csv('COMBAT_metadata_controls.tsv', sep = '\t')

ezclermont = pd.read_csv('ezclermont_summary.tsv', sep  = '\t', names = ["strain", "phylogroup"])
ezclermont_controls = pd.read_csv('ezclermont_summary_controls.tsv', sep  = '\t', names = ["strain", "phylogroup"])

snp_comparisons_clonal = snp_comparisons.query('comparison == "within_traveler_between_timepoint" & SNPs_per_mbp <= 17')

snp_comparisons_clonal = pd.merge(left = snp_comparisons_clonal, right = ezclermont, left_on = 'strain2', right_on = 'strain')

snp_comparisons_clonal = snp_comparisons_clonal[['traveler1', 'phylogroup']].drop_duplicates()

ezclermont_controls_travelers = pd.merge(left = ezclermont_controls, right = metadata_controls, on = 'strain')

phylogroups_longterm = pd.DataFrame(snp_comparisons_clonal["phylogroup"]).assign(type='longterm')
phylogroups_shortterm = pd.DataFrame(ezclermont_controls_travelers["phylogroup"]).assign(type='shortterm')

phylogroups = phylogroups_longterm.append(phylogroups_shortterm)

out = phylogroups.groupby('type')['phylogroup'].value_counts().unstack().fillna(0).transpose()

out.to_csv('phylogroup_comparison.tsv', sep = '\t', float_format='%.f')
