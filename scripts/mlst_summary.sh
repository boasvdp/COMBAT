#!/bin/bash

set -euo pipefail

echo -e "strain\tst"

for file in $@
do
	START=$(awk '{print $1}' "$file")
	NAME=$(basename "$START" .fasta)
	ST=$(awk '{print $3}' "$file")
	echo -e "$NAME\t$ST"
done
