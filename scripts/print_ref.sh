#!/usr/bin/bash

set -euo pipefail

if [[ $# -ne 2 ]]
then
	echo "Please supply two arguments. Usage:"
	echo "bash $0 [SAMPLE] [MLST_REFERENCES FILE]"
	exit 1
fi

SAMPLE=$1
MLST=$2

REFNAME=$(awk -F "\t" -v S="$SAMPLE" 'NR > 1 && $1 == S {print $3}' "$MLST")

ls references/"$REFNAME".fna
