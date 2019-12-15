#!/bin/bash

THRESHOLD1=$1 # threshold marking upper limit of very likely persistence
THRESHOLD2=$2 # threshold marking upper limit of likely persistence

echo "Traveler	SNPs_per_mbp	Likelihood_of_persistence"

awk -F "\t" -v T1=$THRESHOLD1 '$7 == "within_traveler_between_timepoint" && $10 <= T1 {print $3}' snp_comparisons.tsv | sort | uniq | while read traveler
do
	SNP=$(awk -v T=$traveler '$3 == T && $7 == "within_traveler_between_timepoint" {print $10}' snp_comparisons.tsv | sort -n | head -n 1)
	echo "$traveler	$SNP	Very_likely_clonal"
done

awk -F "\t" -v T1=$THRESHOLD1 -v T2=$THRESHOLD2 '$7 == "within_traveler_between_timepoint" && $10 > T1 && $10 <= T2 {print $3}' snp_comparisons.tsv | sort | uniq | while read traveler
do
	SNP=$(awk -v T=$traveler '$3 == T && $7 == "within_traveler_between_timepoint" {print $10}' snp_comparisons.tsv | sort -n | head -n 1)
	if (( $(echo "$SNP > $THRESHOLD1" |bc -l) ))
	then
		echo "$traveler	$SNP	Likely_clonal"
	fi
done

awk -F "\t" -v T2=$THRESHOLD2 '$7 == "within_traveler_between_timepoint" && $10 > T2 {print $3}' snp_comparisons.tsv | sort | uniq | while read traveler
do
	SNP=$(awk -v T=$traveler '$3 == T && $7 == "within_traveler_between_timepoint" {print $10}' snp_comparisons.tsv | sort -n | head -n 1)
	if (( $(echo "$SNP > $THRESHOLD2" |bc -l) ))
	then
		echo "$traveler	$SNP	Not_clonal"
	fi
done
