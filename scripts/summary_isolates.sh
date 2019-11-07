#!/bin/bash

echo "Isolate,pct_classified_reads,pct_Enterobacteriaceae_reads,pct_E.coli_reads,Suspected_species,MLST,Phylogroup,Predicted_serotype,Contigs,N50,Assembly_length,Prokka_CDS,ABRicate_NCBI,ABRicate_ECOH"
for file in $@
do
	file=${file##*/}
	KRAKENROOT=`awk '$5 == "1" {print $1}' kraken_out/"$file"_kraken2_report.txt`
	KRAKENENTERO=`awk '$5 == "543" {print $1}' kraken_out/"$file"_kraken2_report.txt`
	KRAKENECOLI=`awk '$5 == "562" {print $1}' kraken_out/"$file"_kraken2_report.txt`
	KRAKENSPECIES=`awk '$4 == "S" {print $6" "$7}' kraken_out/"$file"_kraken2_report.txt | head -1`
	PHYLOGROUP=`awk '{print $2}' ezclermont_out/"$file".tsv`
	MLST=`awk '{print $3}' mlst/"$file".tsv`
	SEROTYPE=`awk 'NR>1 {print $5}' abricate_out/ecoh/"$file"_ecoh.tsv | cut -d '-' -f 2 | sort -r | uniq | tr "\n" ":" | sed 's/.$//'`
	NRCONTIGS=`grep '# contigs (>= 0 bp' quast_out/"$file"/report.tsv|awk '{print $6}'`
	N50=`grep '^N50' quast_out/"$file"/report.tsv|awk '{print $2}'`
	LENGTH=`grep '^Total length (>= 0 bp' quast_out/"$file"/report.tsv|awk '{print $6}'`
	ABRICATENCBI=`awk 'NR>1 {print $5}' abricate_out/ncbi/"$file"_ncbi.tsv | tr "\n" "|" | sed 's/.$//'`
	ABRICATEECOH=`awk 'NR>1 {print $5}' abricate_out/ecoh/"$file"_ecoh.tsv | tr "\n" "|" | sed 's/.$//'`
	PROKKA=`awk '/CDS/ {print $2}' prokka_out/"$file"/"$file".txt`
	echo "${file},${KRAKENROOT},${KRAKENENTERO},${KRAKENECOLI},${KRAKENSPECIES},${MLST},${PHYLOGROUP},${SEROTYPE},${NRCONTIGS},${N50},${LENGTH},${PROKKA},${ABRICATENCBI},${ABRICATEECOH}"
done
