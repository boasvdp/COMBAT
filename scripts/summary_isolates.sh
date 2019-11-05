#!/bin/bash

echo "Isolate,pct_classified_reads,pct_Enterobacteriaceae_reads,pct_E.coli_reads,Suspected_species,MLST,Phylogroup,Predicted_pathotype,Predicted_serotype,Contigs_SPAdes,N50_SPAdes,Assembly_length_SPAdes,Contigs_SKESA,N50_SKESA,Assembly_length_SKESA,Prokka_SPAdes,Prokka_SKESA,ABRicate_NCBI,ABRicate_ECOH"
for file in $@
do
	file=${file##*/}
	KRAKENROOT=`awk '$5 == "1" {print $1}' kraken_out/"$file"_kraken2_report.txt`
	KRAKENENTERO=`awk '$5 == "543" {print $1}' kraken_out/"$file"_kraken2_report.txt`
	KRAKENECOLI=`awk '$5 == "562" {print $1}' kraken_out/"$file"_kraken2_report.txt`
	KRAKENSPECIES=`awk '$4 == "S" {print $6" "$7}' kraken_out/"$file"_kraken2_report.txt | head -1`
	PHYLOGROUP=`awk -F "\t" '{print $5}' clermontyper_out/"$file"/*_phylogroups.txt | head -1`
	MLST=`awk '{print $3}' mlst/"$file".tsv`
	PATHOTYPE=`awk '{print $1}' patho_typing_out/"$file"/patho_typing.report.txt`
	SEROTYPE=`awk 'NR>1 {print $5}' abricate_out/ecoh/"$file"_ecoh.tsv | cut -d '-' -f 2 | sort -r | uniq | tr "\n" ":" | sed 's/.$//'`
	NRCONTIGSSPADES=`grep '# contigs (>= 0 bp' quast_out/spades/"$file"/report.tsv|awk '{print $6}'`
	N50SPADES=`grep '^N50' quast_out/spades/"$file"/report.tsv|awk '{print $2}'`
	LENGTHSPADES=`grep '^Total length (>= 0 bp' quast_out/spades/"$file"/report.tsv|awk '{print $6}'`
	NRCONTIGSSKESA=`grep '# contigs (>= 0 bp' quast_out/skesa/"$file"/report.tsv|awk '{print $6}'`
	N50SKESA=`grep '^N50' quast_out/skesa/"$file"/report.tsv|awk '{print $2}'`
	LENGTHSKESA=`grep '^Total length (>= 0 bp' quast_out/skesa/"$file"/report.tsv|awk '{print $6}'`
	ABRICATENCBI=`awk 'NR>1 {print $5}' abricate_out/ncbi/"$file"_ncbi.tsv | tr "\n" "|" | sed 's/.$//'`
	ABRICATEECOH=`awk 'NR>1 {print $5}' abricate_out/ecoh/"$file"_ecoh.tsv | tr "\n" "|" | sed 's/.$//'`
	PROKKASPADES=`awk '/CDS/ {print $2}' prokka_out/spades/"$file"_spades/"$file".txt`
	PROKKASKESA=`awk '/CDS/ {print $2}' prokka_out/skesa/"$file"_skesa/"$file".txt`
	echo "${file},${KRAKENROOT},${KRAKENENTERO},${KRAKENECOLI},${KRAKENSPECIES},${MLST},${PHYLOGROUP},${PATHOTYPE},${SEROTYPE},${NRCONTIGSSPADES},${N50SPADES},${LENGTHSPADES},${NRCONTIGSSKESA},${N50SKESA},${LENGTHSKESA},${PROKKASPADES},${PROKKASKESA},${ABRICATENCBI},${ABRICATEECOH}"
done
