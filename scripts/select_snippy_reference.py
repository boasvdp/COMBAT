#!/usr/bin/env python3

import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Select and print Snippy samples for snippy-cre.')

parser.add_argument('info', help="MLST summary file with strain names and ST", type=str)
parser.add_argument('prefix', help="ST to select", type=str)
parser.add_argument('refdir', help="Directory containing Snippy references", default='references', type=str)

args = parser.parse_args()

# Open MLST summary file and read into lines object
c = open(args.info)
lines = c.readlines()
c.close()

# Initiate empty print string
list_ref = []

# Remove trailing slash (is present) from snippy directory
refdir = args.refdir.rstrip('/')

# Loop through MLST summary and select strains with right ST
for line in lines:
  if line.split('\t')[1] == args.prefix:
    add_to_list_ref = refdir + '/' + line.split('\t')[2].rstrip('\n') + '.fna'
    list_ref.append(add_to_list_ref)

# Remove leading space
list_ref = str(set(list_ref)).replace("'", '').replace('{', '').replace('}', '')

# Print selected samples
print(list_ref)
