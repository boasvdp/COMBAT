# COMBAT
Snakemake pipeline for the analysis of COMBAT data

## Introduction
This repository was created for the analysis of a COMBAT spin-off study. The original COMBAT study can be found [here](https://doi.org/10.1016/S1473-3099(16)30319-X). This study investigated acquisition and persistence of extended-spectrum beta-lactamase-positive Enterobacteriaceae (ESBL-E) in a population of 2001 Dutch international travellers. The original study investigated both acquisition and persistence of ESBL-E when returning from travel. 633 travellers acquired ESBL-E during international travel. Of these 633, 38 travellers harboured acquired ESBL-E (all E. coli, hence ESBL-Ec) for more than a year after return. In this analysis, using WGS we identify:

1) Which ESBL-Ec actually persisted;
2) How these ESBL genes persisted (did the whole strain persist, or only the plasmid?);
3) Associated bacterial characteristics for persistence.

## Pipeline
The pipeline is schematically picture below. Configuration can be managed through the `config.yaml` file.

![rulegraph](https://github.com/boasvdp/COMBAT/blob/master/rulegraph.svg "Rulegraph")

The input consists of Illumina and Nanopore sequencing data, which can be downloaded from ENA/SRA. The easiest way of running the pipeline is to use conda as follows:

```bash
snakemake --use-conda -j 16
```

Where the `-j` flag indicates how many cores the pipeline should use.


Default settings were used unless noted otherwise.

### Rules included in the pipeline

-	Fastp

Fastp (version xxx) was used to trim adapters and low quality bases from Illumina sequencing reads. Reads were not filtered based on length. Output were trimmed reads, together with a report in .html and .json formats.

-	Kraken2

Kraken2 (version xxx) was used to detect contamination in trimmed Illumina reads (from fastp). The Minikraken 8GB database without human sequences was used in this pipeline. Output consisted of a report.

-	Shovill

Shovill (version xxx) was used to assemble trimmed reads (from fastp) into draft assemblies, using SPAdes. Reads were subsampled to a maximum depth of 100X. Only contigs longer than 500 bp were outputted.

-	Quast

Quast (version xxx) was used to evaluate quality of draft genome assemblies constructed by Shovill.

-	Multiqc

MultiQC (version xxx) was used to aggregate reports from fastp and Quast.

-	Amrfinder

AMRFinderplus was used to detect resistance genes and mutations from draft genome assemblies (from Shovill). Specific annotation for the Escherichia genus was used.

-	Mlst

Mlst (version xxx) was used to identify the multi locus sequence type based on the draft genome assemblies from Shovill.

-	Prokka

Prokka (version xxx) was used to annotate draft genome assemblies from Shovill. Specific annotation for the Escherichia genus was used.

-	Ezclermont

Ezclermont (version xxx) was used to predict the phylogroup from the draft genome assemblies (from Shovill).

-	Snippy

Snippy (version xxx) was used to map trimmed reads (from fastp) onto a complete reference genome (E. coli ATCC25922, Genbank: CP009072). 

-	Snippycore

Snippy-core, implemented in snippy (version xxx) was used to obtain the core gene alignment from all Snippy runs together.

-	IQtree

IQtree (version xxx) was used to infer the phylogeny from the core genome alignment (from Snippy). Snp-sites (version xxx) was used to infer the number of constant sites per nucleotide which was inputted to IQtree via the –fconst flag. Modelfinder, implemented in IQtree, advised to use the “TVMe+ASC+R4” model.

-	ClonalFrameML

ClonalFrameML (version xxx, REF) was used to detect recombination events based on the core genome alignment from Snippy and the phylogeny from IQtree.

-	Maskrc

Maskrc-svg (version xxx) was used to mask recombination events identified by ClonalFrameML.

-	snp_comparisons

SNP comparisons were made using snp-dists (version xxx, REF) and custom Python and bash scripts available at https://github.com/boasvdp/COMBAT/blob/master/scripts/. SNPs were expressed as SNPs/1,000,000 aligned bases.

-	plot_snp_comparisons

SNP comparisons were plotted in R using ggplot2. Script available at https://github.com/boasvdp/COMBAT/blob/master/scripts/.

-	print_travelers

Script to output a list of travellers with data on how the ESBL genes persisted in their gut based on SNP data. Available at https://github.com/boasvdp/COMBAT/blob/master/scripts/

-	fastqc

FastQC (version xxx) was used to assess Nanopore read quality.

-	filtlong

Filtlong (version xxx) was used to filter ONT MinION reads based on read length and read identity, using trimmed Illumina reads (from fastp) as a reference. Depth of sequencing was capped at 100X.

-	unicycler

Unicycler (version xxx) was used to perform assembly of trimmed Illumina reads (from fastp) and trimmed ONT reads (from Filtlong).

-	amrfinder_nanopore

AMRFinderplus (version xxx) was used to detect ESBL genes in the complete genome assemblies (from Unicycler).

-	get_ESBL_contigs

Seqtk (version xxx, REF) was used to extract contigs from complete genome assemblies (from Unicycler) which contain ESBL genes (from amrfinder_nanopore).

-	prokka_ESBL_contigs

Prokka (version xxx) was used to annotate the ESBL plasmids (from get_ESBL_contigs) so that these can be used for the ANIcalculator webtool.

