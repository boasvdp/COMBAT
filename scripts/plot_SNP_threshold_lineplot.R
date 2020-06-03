#!/usr/bin/env Rscript

# Read command line arguments
args = commandArgs(trailingOnly=TRUE)

library(reshape2)
library(ggplot2)

output_basename <- args[3]

##### withST38 #####
output_name = paste(output_basename, "_thresholds_lineplot_withST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 11)

df <- read.delim(args[1])

df_casted <- dcast(df, formula = SNP_threshold + Method ~ Comparison)

name_x <- "Number of isolate pairs from different travelers assigned as 'clonal'"
name_y <- "Number of isolate pairs from the same traveler and timepoint assigned as 'clonal'"

ggplot(df_casted, aes(x = between_traveler, y = within_traveler_within_timepoint, group = Method, color = Method, label = SNP_threshold)) +
	geom_line() +
	geom_point() +
	geom_text(hjust=0, vjust=0, color = "black") +
	scale_y_continuous(name = name_y) +
	scale_x_continuous(name = name_x, limits = c(0,100)) +
	scale_colour_discrete(name = "Method", labels = c("Core genome, no gaps", "Pairwise comparison, corrected for alignment length", "Pairwise comparison, uncorrected"))

dev.off()

##### withoutST38 #####
output_name = paste(output_basename, "_thresholds_lineplot_withoutST38.pdf", sep = "")
pdf(file = output_name, height = 6, width = 11)

df <- read.delim(args[2])

df_casted <- dcast(df, formula = SNP_threshold + Method ~ Comparison)

name_x <- "Number of isolate pairs from different travelers assigned as 'clonal'"
name_y <- "Number of isolate pairs from the same traveler and timepoint assigned as 'clonal'"

ggplot(df_casted, aes(x = between_traveler, y = within_traveler_within_timepoint, group = Method, color = Method, label = SNP_threshold)) +
	geom_line() +
	geom_point() +
	geom_text(hjust=0, vjust=0, color = "black") +
	scale_y_continuous(name = name_y) +
	scale_x_continuous(name = name_x, limits = c(0,100)) +
	scale_colour_discrete(name = "Method", labels = c("Core genome, no gaps", "Pairwise comparison, corrected for alignment length", "Pairwise comparison, uncorrected"))

dev.off()
