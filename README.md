# COMBAT
Snakemake pipeline for the analysis of COMBAT data

## Introduction
This repository was created for the analysis of a COMBAT spin-off study. The original COMBAT study can be found [here](https://doi.org/10.1016/S1473-3099(16)30319-X). This study investigated acquisition and persistence of extended-spectrum beta-lactamase-positive Enterobacteriaceae (ESBL-E) in a population of 2001 Dutch international travellers. The original study investigated both acquisition and persistence of ESBL-E when returning from travel. 633 travellers acquired ESBL-E during international travel. Of these 633, 38 travellers harboured acquired ESBL-E (all E. coli, hence ESBL-Ec) for more than a year after return. In this analysis, using WGS we identify:

1) Which ESBL-Ec actually persisted;
2) How these ESBL genes persisted (did the whole strain persist, or only the plasmid?);
3) Associated bacterial characteristics for persistence.

## Pipeline
The pipeline is schematically picture below. Configuration can be managed through the `config.yaml` file.

![rulegraph](https://github.com/boasvp/COMBAT/rulegraph.svg "Rulegraph")

The input consists of Illumina and Nanopore sequencing data, which can be downloaded from ENA/SRA. The easiest way of running the pipeline is to use conda as follows:

```bash
snakemake --use-conda -j 16
```

Where the `-j` flag indicates how many cores the pipeline should use.
