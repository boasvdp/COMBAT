#!/usr/bin/env python3

import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Select and print Snippy samples for snippy-cre.')

parser.add_argument('info', help="MLST summary file with strain names and ST", type=str)
parser.add_argument('prefix', help="ST to select", type=str)
parser.add_argument('snippydir', help="Directory containing Snippy output", default='snippy_out', type=str)

args = parser.parse_args()

# Open MLST summary file and read into lines object
c = open(args.info)
lines = c.readlines()
c.close()

# Initiate empty print string
printstr = ''

# Remove trailing slash (is present) from snippy directory
snippydir = args.snippydir.rstrip('/')

# Loop through MLST summary and select strains with right ST
for line in lines:
  if line.split('\t')[1] == args.prefix:
    add_to_printstr = snippydir + '/' + line.split('\t')[0]
    printstr = str(printstr + ' ' + add_to_printstr)

# Remove leading space
printstr = printstr.lstrip(' ')

# Print selected samples
print(printstr)
