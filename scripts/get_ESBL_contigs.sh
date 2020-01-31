#!/bin/bash

mkdir -p ESBL_contigs

for strain in $@
do
	NAME=$(basename $strain .tsv)
	awk -F "\t" '$6 ~ "blaCTX" {print $2}' amrfinder_nanopore_out/${NAME}.tsv > tmp.list
	seqtk subseq unicycler_out/${NAME}/assembly.fasta tmp.list > ESBL_contigs/${NAME}.fasta
	rm tmp.list
done
