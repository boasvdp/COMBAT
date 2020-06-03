#!/usr/bin/env Rscript

# Read command line arguments
args = commandArgs(trailingOnly=TRUE)

# library
library(ggplot2)
library(ggthemes)

# Read plot data
plot <- read.delim(args[1])

# Define labels
phylolabels <- c("A", "B1", "B2", "C", "Cryptic\nclade I", "D", "F", "Multiple")

# Define breaks


pdf(file = args[2], height = 4, width = 8)

# Grouped
ggplot(plot, aes(fill=type, y=count, x=phylogroup)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_clean() +
  scale_x_discrete(name = "Phylogroup/lineage", labels = phylolabels) +
  scale_y_continuous(name = "Number of strains", breaks = seq(0,10,2), labels = seq(0,10,2)) +
  scale_fill_discrete(name = "Carriage duration of \ntravel-acquired strains", labels = c("Long: >12 months (n = 14)", "Short: <1 month (n = 27)")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))

dev.off()
