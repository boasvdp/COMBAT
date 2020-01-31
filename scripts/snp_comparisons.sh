#!/bin/bash

rm -rf tmp.list || true

echo "strain1	strain2	traveler1	timepoint1	traveler2	timepoint2	comparison	total_SNPs	alignment_length	SNPs_per_mbp" > $1

awk 'NR > 1 {print $0}' comparisons_intermediate.tsv | while read strain1 strain2 traveler1 timepoint1 traveler2 timepoint2 comparison
do
	echo "$strain1,$strain2" | tr "," "\n" > tmp.list
	seqtk subseq masked.aln tmp.list > "${strain1}_${strain2}.aln"
	SNP=$(snp-dists "${strain1}_${strain2}.aln" | awk -F "\t" 'NR == 2 {print $3}')
	ALNLENGTH=$(snp-sites -cb "${strain1}_${strain2}.aln" | wc -L)
	SNP_PER_MBP=$(echo "$SNP / $ALNLENGTH * 1000000" | bc -l)
	echo "$strain1	$strain2	$traveler1	$timepoint1	$traveler2	$timepoint2	$comparison	$SNP	$ALNLENGTH	$SNP_PER_MBP" >> $1
	rm "${strain1}_${strain2}.aln"
done

rm tmp.list || true
rm comparisons_intermediate.tsv || true
