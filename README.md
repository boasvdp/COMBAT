# COMBAT
Snakemake pipeline for the analysis of COMBAT data

## Introduction
This repository was created for the analysis of a COMBAT spin-off study. The original COMBAT study can be found [here](https://doi.org/10.1016/S1473-3099(16)30319-X). This study investigated acquisition and persistence of extended-spectrum beta-lactamase-positive Enterobacteriaceae (ESBL-E) in a population of 2001 Dutch international travellers. The original study investigated both acquisition and persistence of ESBL-E when returning from travel. 633 travellers acquired ESBL-E during international travel. Of these 633, 38 travellers harboured acquired ESBL-E (all E. coli, hence ESBL-Ec) for more than a year after return. In this analysis, using WGS we identify:

1) Which ESBL-Ec actually persisted;
2) How these ESBL genes persisted (did the whole strain persist, or only the plasmid?);
3) Associated bacterial characteristics for persistence.

## Pipeline
The pipeline is schematically pictured below. Configuration can be managed through the `config.yaml` file.

![rulegraph](https://github.com/boasvdp/COMBAT/blob/master/rulegraph.svg "Rulegraph")

The input consists of Illumina and Nanopore sequencing data, which can be downloaded from ENA/SRA. The necessary pseudonimised metadata is provided in this repository. The easiest way of running the standard pipeline is to use snakemake and conda as follows:

```bash
snakemake --use-conda -j 16
```

Where the `-j` flag indicates how many cores the pipeline should use.

We provide an additional analysis which repeats the main SNP calling and comparisons part of the study. This additional analysis is provided in a separate Snakefile `Snakefile_alternative_references`. To run, use:

```bash
snakemake --use-conda -j 16 -s Snakefile_alternative_references
```

Default settings were used unless noted otherwise.

## Summary of the pipeline

The pipeline consists of 4 parts:

1. Basic bacterial genomics workflow (mainly typing)
2. SNP calling and comparisons
3. Plasmid comparisons
4. Comparing findings with a matched control group

### 1. Basic bacterial genomics workflow (mainly typing)

We included all 38 long-term carriers (>= 12 months carriage after return from travel) originally identified in the COMBAT study. Long-term carriage in the COMBAT study was originally defined as travellers which carried the same ESBL gene **group** (as defined by microarray) at all timepoints after travel (0, 1, 3, 6 and 12 months after return). After employing WGS, we found that six travellers did not carry the same ESBL gene (but the same ESBL gene **group**) between timepoints. These travellers were excluded from the long-term carriers, as were four travellers that had missing samples. The final dataset comprised of 28 long-term carriers. Isolates from these 28 long-term carriers were sequencing using Illumina HiSeq and data was filtered using fastp (version 0.20.0, with `--disable_length_filtering`). Trimmed reads were subsequently processed using Kraken2 (version 2.0.8_beta, using the Minikraken 8GB database without human sequences) which showed all samples were *E. coli*. Shovill (version 1.0.9) using SPAdes (version 3.13.1) assembled contigs from the trimmed reads subsampled to a theoretical depth of 100X. Only contigs of 500 bp and longer were retained. Quast (version 4.6.3) was used to evaluate assembly quality. MultiQC (version 1.6) was finally used to aggregated reports from fastp and Quast.

For typing and annotation, we employed several tools. First, we used AMRfinderplus (version 3.2.3, with `--organism Escherichia`) on Shovill assemblies to find acquired resistance genes and *Escherichia* point mutations associated with resistance. We used prokka (version 1.14) for general annotation, again using an *Escherichia* specific database (`--genus Escherichia --usegenus`). mlst (version 2.17.6) and EzClermont (version 0.4.3) were used for sequence type and phylogroup prediction, respectively.

### 2. SNP calling and comparisons

To identify which strains have persisted in the long-term carriers, we compare isolates on single nucleotide polymorphisms (SNPs). SNPs are identified by mapping the trimmed reads on an *E. coli* reference genome (ATCC 25922, Genbank: CP009072) using Snippy (version 4.4.5). Mapped samples were combined using `snippy-core`. To mask recombinatory regions, maskrc-svg (version 0.5) is used based on ClonalFrameML (version 1.12) output. The input phylogeny for ClonalFrameML was inferred using IQtree (version 1.6.12), snp-sites (version 2.5.1, with option `-C` as input for IQtree's `-fconst` argument) and ModelFinder packaged with IQtree.

From the full alignment to ATCC 25922, SNP differences per 1,000,000 aligned bases were calculated using snp-dists (version 0.7.0) and a modified version of snp-dists (https://github.com/boasvdp/snp-dists). Python scripts (using pandas version 0.25.3, numpy version 1.17.3 and itertools version 8.0.2 libraries) available in this repo were used to parse results. Results were plotted using ggplot2 (version 3.1.1) and ggthemes (version 4.2.0). Scripts available in `scripts`.

We identified fourteen long-term carriers which harboured a single strain for 12 months. For nine travellers, we could not reliably differentiate between isolated ST38 clones. Five long-term carriers did not harbour a single persistent strain for 12 months (see point 3 below). We repated the SNP analysis using two other reference genomes (DSM 30083 and K-12) and arrived at the same conclusions. This analysis is provided in an extra Snakefile, `Snakefile_alternative_references`. 

### 3. Plasmid comparisons

In five long-term carriers, we did not identify any isolate pairs that could recently originate from a single strain. One traveller had an isolate pair which shared 84 SNPs/Mbp, while the other four had isolate pairs at least 900 SNPs/Mbp apart. Either the ESBL gene shifted between bacterial hosts, or the travellers lost the strain(s) acquired during travel and reacquired unrelated strains with the same ESBL gene. To check this, we employed Oxford Nanopore Technologies long-read sequencing which allows us to resolve plasmid structures. We sequenced one isolate per traveller per timepoint on the MinION platform and assessed read quality using fastQC (version 0.11.98). Subsequently we filtered reads using filtlong (version 0.2.0). Trimmed Illumina reads were used as reference. We retained a 100X coverage of Nanopore reads. Next, we assembled the genomes using trimmed Nanopore and trimmed Illumina data, using Unicycler (version 0.4.8). We identified contigs which carried resistance genes using AMRfinderplus (version 3.2.3) and extracted these contigs. Nine out of ten isolates harboured plasmid-located ESBL genes, and one isolate harboured a chromosomally-inserted ESBL gene. We annotated the ESBL contigs using prokka (version 1.14) and compared contigs using ANIcalculator and Easyfig (not in this pipeline).

### 4. Comparing findings with a matched control group

To find characteristics associated with long-term carriage, we matched 27 short-term carriers (<1 month carriage after acquisition during travel) with 14 long-term carriers on age, sex and travel destination. One long-term carrier could only be matched to one short-term carriers, instead of the intended two short-term carriers. We processed isolates from short-term carriers the same as the isolates from long-term carriers and compared these on sequence type (mlst), phylogroup (EzClermont) and ExPEC status (linking ST to known ExPEC lineages). We find ExPEC status is strongly associated to long-term carriage of travel-acquired ESBL+ *Escherichia coli*, mainly driven by phylogroup D (ST69, ST393, ST405) and ST131 specifically. 
