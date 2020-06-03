#!/usr/bin/env Rscript

# Read command line arguments
args = commandArgs(trailingOnly=TRUE)


library(ggplot2)

snp_comparisons_withST38 <- read.delim(args[1])
snp_comparisons_withoutST38 <- read.delim(args[2])
output_basename <- args[3]

##### with ST38 #####

### Plot 0-50 SNPs
output_name = paste(output_basename, "_histogram_corrected_50_withST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withST38, aes(x = SNPs_corrected, fill = comparison)) +
	geom_histogram(bins = 50) +
	scale_x_continuous(limits = c(0,50), name = "SNPs per 1,000,000 bp") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot 0-500 SNPs
output_name = paste(output_basename, "_histogram_corrected_500_withST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withST38, aes(x = SNPs_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(limits = c(0,500), name = "SNPs per 1,000,000 bp") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot all SNPs
output_name = paste(output_basename, "_histogram_corrected_all_withST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 20)

ggplot(snp_comparisons_withST38, aes(x = SNPs_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(name = "SNPs per 1,000,000 bp") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot 0-50 SNPs
output_name = paste(output_basename, "_histogram_not_corrected_50_withST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withST38, aes(x = SNPs_not_corrected, fill = comparison)) +
	geom_histogram(bins = 50) +
	scale_x_continuous(limits = c(0,50), name = "SNPs") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot 0-500 SNPs
output_name = paste(output_basename, "_histogram_not_corrected_500_withST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withST38, aes(x = SNPs_not_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(limits = c(0,500), name = "SNPs") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot all SNPs
output_name = paste(output_basename, "_histogram_not_corrected_all_withST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 20)

ggplot(snp_comparisons_withST38, aes(x = SNPs_not_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(name = "SNPs") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

##### Without ST38 #####

### Plot 0-50 SNPs
output_name = paste(output_basename, "_histogram_corrected_50_withoutST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withoutST38, aes(x = SNPs_corrected, fill = comparison)) +
	geom_histogram(bins = 50) +
	scale_x_continuous(limits = c(0,50), name = "SNPs per 1,000,000 bp") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot 0-500 SNPs
output_name = paste(output_basename, "_histogram_corrected_500_withoutST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withoutST38, aes(x = SNPs_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(limits = c(0,500), name = "SNPs per 1,000,000 bp") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot all SNPs
output_name = paste(output_basename, "_histogram_corrected_all_withoutST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 20)

ggplot(snp_comparisons_withoutST38, aes(x = SNPs_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(name = "SNPs per 1,000,000 bp") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot 0-50 SNPs
output_name = paste(output_basename, "_histogram_not_corrected_50_withoutST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withoutST38, aes(x = SNPs_not_corrected, fill = comparison)) +
	geom_histogram(bins = 50) +
	scale_x_continuous(limits = c(0,50), name = "SNPs") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot 0-500 SNPs
output_name = paste(output_basename, "_histogram_not_corrected_500_withoutST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 10)

ggplot(snp_comparisons_withoutST38, aes(x = SNPs_not_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(limits = c(0,500), name = "SNPs") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot all SNPs
output_name = paste(output_basename, "_histogram_not_corrected_all_withoutST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 20)

ggplot(snp_comparisons_withoutST38, aes(x = SNPs_not_corrected, fill = comparison)) +
	geom_histogram(bins = 100) +
	scale_x_continuous(name = "SNPs") +
	labs(fill = "Comparison") +
	scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

