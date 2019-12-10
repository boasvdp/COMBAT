#!/bin/bash

rm -rf tmp.list || true

echo "strain1	strain2	traveler_strain1	traveler_strain2	timepoint_strain1	timepoint_strain2	total_SNPs	aln_length	SNPs_per_mbp" > $1

while read traveler1 strain1 timepoint1
do
	while read traveler2 strain2 timepoint2
	do
		PRESENCE=$(awk -F "\t" -v ONE=$strain1 -v TWO=$strain2 '($1 == ONE || $2 == ONE) && ($1 == TWO || $2 == TWO) {t += 1} END {print t}' $1)
		if [[ $PRESENCE -eq 0 ]]
		then
			if [[ $strain1 != $strain2 ]]
			then
				echo "$strain1,$strain2" | tr "," "\n" > tmp.list
				seqtk subseq masked.aln tmp.list > "${strain1}_${strain2}.aln"
				SNP=$(snp-dists "${strain1}_${strain2}.aln" | awk -F "\t" 'NR == 2 {print $3}')
				ALNLENGTH=$(snp-sites -cb "${strain1}_${strain2}.aln" | wc -L)
				SNP_PER_MBP=$(echo "$SNP / $ALNLENGTH * 1000000" | bc -l)
				echo "$strain1	$strain2	$traveler1	$traveler2	$timepoint1	$timepoint2	$SNP	$ALNLENGTH	$SNP_PER_MBP" >> $1
				rm "${strain1}_${strain2}.aln"
			fi
		fi
	done < COMBAT_metadata.tsv
done < COMBAT_metadata.tsv

rm tmp.list || true
