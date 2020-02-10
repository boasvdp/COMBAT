#!/bin/bash
set -eo pipefail

rm -rf tmp.list || true

echo "strain1	strain2	traveler1	timepoint1	traveler2	timepoint2	comparison	alignment_length	SNPs_no_gaps	SNPs_not_corrected	SNPs_corrected" > $1

awk 'NR > 1 {print $0}' snp_comparisons_intermediate.tsv | while read strain1 strain2 traveler1 timepoint1 traveler2 timepoint2 comparison SNPs_no_gaps SNPs_not_corrected
do
	echo "$strain1,$strain2" | tr "," "\n" > tmp.list
	seqtk subseq masked.aln tmp.list > "${strain1}_${strain2}.aln"
	ALNLENGTH=$(snp-sites -cb "${strain1}_${strain2}.aln" | wc -L)
	SNPs_corrected=$(echo "$SNPs_not_corrected / $ALNLENGTH * 1000000" | bc -l)
	echo "$strain1	$strain2	$traveler1	$timepoint1	$traveler2	$timepoint2	$comparison	$ALNLENGTH	$SNPs_no_gaps	$SNPs_not_corrected	$SNPs_corrected" >> $1
	rm "${strain1}_${strain2}.aln"
done

rm tmp.list || true
rm snp_comparisons_intermediate.tsv || true
