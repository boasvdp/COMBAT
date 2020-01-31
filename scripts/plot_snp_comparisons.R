library(ggplot2)

snp_comparisons <- read.delim("snp_comparisons.tsv")

### Plot 0-50 SNPs
pdf(file = "plots/snp_comparisons_50.pdf", height = 6, width = 10)

ggplot(snp_comparisons, aes(x = SNPs_per_mbp, fill = comparison)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(limits = c(0,50), name = "SNPs per 1,000,000 bp") +
  labs(fill = "Comparison") +
  scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot 0-500 SNPs
pdf(file = "plots/snp_comparisons_500.pdf", height = 6, width = 10)

ggplot(snp_comparisons, aes(x = SNPs_per_mbp, fill = comparison)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(limits = c(0,500), name = "SNPs per 1,000,000 bp") +
  labs(fill = "Comparison") +
  scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()

### Plot all SNPs
pdf(file = "plots/snp_comparisons_full.pdf", height = 6, width = 20)

ggplot(snp_comparisons, aes(x = SNPs_per_mbp, fill = comparison)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(name = "SNPs per 1,000,000 bp") +
  labs(fill = "Comparison") +
  scale_fill_discrete(labels = c("Between traveler", "Within traveler, between timepoint", "Within traveler, within timepoint"))

dev.off()
