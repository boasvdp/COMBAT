configfile: "config.yaml"

import pandas as pd

# Read isolate list, ST38 excluded
df = pd.read_csv('lists/list_isolates_withoutST38.txt', sep = '\t')
WITHOUTST38 = tuple(df)

# Read control isolate list, ST38 excluded
df_controls = pd.read_csv('lists/list_isolates_controls_withoutST38.txt', sep = '\t')
CONTROLSWITHOUTST38 = tuple(df_controls)

# Other ID tuples can be globbed from input dirs
NANOPORE, = glob_wildcards("raw_nanopore/{id}.fastq.gz")
WITHST38, = glob_wildcards("raw_illumina/{id}_1.fastq.gz")
#CONTROLS, = glob_wildcards("controls/raw_illumina/{id}_1.fastq.gz")
CONTROLSWITHST38, = glob_wildcards("controls/raw_illumina/{id}_1.fastq.gz")

rule all:
	input:
		expand("kraken_out/{sample}_kraken2_report.txt", sample=WITHST38),
		expand("amrfinder_out/{sample}.tsv", sample=WITHST38),
		expand("mlst/{sample}.tsv", sample=WITHST38),
		expand("prokka_out/{sample}", sample=WITHST38),
		"multiqc_out/multiqc_fastp.html",
		"abricate_summary/summary_ncbi.tsv",
		expand("fastqc_out/{sample}", sample=NANOPORE),
		expand("ESBL_contigs_annotations/{sample}", sample=NANOPORE),
		expand("controls/kraken_out/{sample}_kraken2_report.txt", sample=CONTROLSWITHOUTST38),
		expand("controls/mlst/{sample}.tsv", sample=CONTROLSWITHOUTST38),
		"controls/multiqc_out/multiqc_fastp.html",
		"phylogroup_comparison_withST38.tsv",
		"phylogroup_comparison_withoutST38.tsv",
		"phylogroup_comparison_withST38.pdf",
		"phylogroup_comparison_withoutST38.pdf"

### STANDARD PIPELINE

rule fastp:
	input:
		fw = "raw_illumina/{sample}_1.fastq.gz",
		rv = "raw_illumina/{sample}_2.fastq.gz"
	output:
		fw = "trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "trimmed_illumina/{sample}_2_AT_QT.fastq.gz",
		html = "fastp_out/html/{sample}_AT_QT_fastp.html",
		json = "fastp_out/json/{sample}_AT_QT_fastp.json"
	conda:
		"envs/fastp.yaml"
	params:
		compression_level = config["fastp"]["compression_level"],
		general = config["fastp"]["general"]
	log:
		"logs/fastp/{sample}.log"
	threads: 8
	shell:
		"""
		fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}
		"""

rule kraken2:
	input:
		fw = "trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "trimmed_illumina/{sample}_2_AT_QT.fastq.gz"
	output:
		report = "kraken_out/{sample}_kraken2_report.txt"
	conda:
		"envs/kraken.yaml"
	params:
		general = config["kraken"]["general"],
		db = config["kraken"]["db"]
	log:
		"logs/kraken2/{sample}.log"
	threads: 8
	shell:
		"""
		kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}
		"""

rule shovill:
	input:
		fw = "trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "trimmed_illumina/{sample}_2_AT_QT.fastq.gz"
	output:
		assembly = "genomes/{sample}.fasta",
		shovill = directory("shovill_out/{sample}")
	params:
		minlen = config["shovill"]["minlen"],
		ram = config["shovill"]["ram"],
		depth = config["shovill"]["depth"],
		assembler = config["shovill"]["assembler"],
		tmpdir = config["shovill"]["tmpdir"]
	conda:
		"envs/shovill.yaml"
	log:
		"logs/shovill/{sample}.log"
	threads: 16
	shell:
		"""
		shovill --assembler {params.assembler} --outdir {output.shovill} --tmp {params.tmpdir} --depth {params.depth} --cpus {threads} --ram {params.ram} --minlen {params.minlen} --R1 {input.fw} --R2 {input.rv} 2>&1>{log}
		cp {output.shovill}/contigs.fa {output.assembly}
		"""

rule quast:
	input:
		"genomes/{sample}.fasta"
	output:
		directory("quast_out/{sample}")
	conda:
		"envs/quast.yaml"
	log:
		"logs/quast/{sample}.log"
	threads: 8
	shell:
		"""
		quast --threads {threads} -o {output} {input} 2>&1>{log}
		"""

rule abricate:
	input:
		assembly = "genomes/{sample}.fasta"
	output:	
		ncbi = "abricate_out/ncbi/{sample}_ncbi.tsv",
		vfdb = "abricate_out/vfdb/{sample}_vfdb.tsv",
		plasmidfinder = "abricate_out/plasmidfinder/{sample}_plasmidfinder.tsv",
		ecoh = "abricate_out/ecoh/{sample}_ecoh.tsv"
	conda:
		"envs/abricate.yaml"
	params:
		minid = config["abricate"]["minid"],
		mincov = config["abricate"]["mincov"],
		ncbi = config["abricate"]["ncbi"],
		vfdb = config["abricate"]["vfdb"],
		plasmidfinder = config["abricate"]["plasmidfinder"],
		ecoh = config["abricate"]["ecoh"]
	log:
		"logs/abricate/{sample}.log"
	threads: 8
	shell:
		"""
		abricate --db {params.ncbi} --nopath --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output.ncbi} 2>>{log}
		abricate --db {params.vfdb} --nopath --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output.vfdb} 2>>{log}
		abricate --db {params.plasmidfinder} --nopath --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output.plasmidfinder} 2>>{log}
		abricate --db {params.ecoh} --nopath --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output.ecoh} 2>>{log}
		"""

rule amrfinder:
	input:
		"genomes/{sample}.fasta"
	output:
		"amrfinder_out/{sample}.tsv"
	conda:
		"envs/amrfinder.yaml"
	params:
		organism = config["amrfinder"]["organism"]
	log:
		"logs/amrfinder/{sample}.log"
	threads: 16
	shell:
		"""
		amrfinder -u
		amrfinder --threads {threads} --nucleotide {input} --organism {params.organism} --output {output} 2>&1>{log}
		"""

rule mlst:
	input:
		"genomes/{sample}.fasta"
	output:
		"mlst/{sample}.tsv"
	conda:
		"envs/mlst.yaml"
	log:
		"logs/mlst/{sample}.log"
	shell:
		"""
		mlst {input} > {output} 2>>{log}
		"""

rule prokka:
	input:
		"genomes/{sample}.fasta"
	output:
		directory("prokka_out/{sample}")
	conda:
		"envs/prokka.yaml"
	params:
		general = config["prokka"]["general"],
		kingdom = config["prokka"]["kingdom"],
		genus = config["prokka"]["genus"],
		species = config["prokka"]["species"],
		prefix = "{sample}"
	log:
		"logs/prokka/{sample}.log"
	threads: 8
	shell:
		"""
		prokka {params.general} --force --outdir {output} --genus {params.genus} --species {params.species} --kingdom {params.kingdom} --cpus {threads} --prefix {params.prefix} {input} 2>&1>{log}
		if [ -f {output}/*.gff ]; then echo "{output} exists"; else exit 1; fi
		"""

rule ezclermont:
	input:
		"genomes/{sample}.fasta"
	output:
		"ezclermont_out/{sample}.tsv"
	conda:
		"envs/ezclermont.yaml"
	log:
		"logs/ezclermont/{sample}.log"
	shell:
		"""
		set +e
		ezclermont {input} 1>{output} 2>{log} || true
		"""

rule ezclermont_summary:
	input:
		expand("ezclermont_out/{sample}.tsv", sample=WITHST38)
	output:
		"ezclermont_summary.tsv"
	log:
		"logs/ezclermont_summary.log"
	shell:
		"""
		cat {input} 1> {output} 2>{log}
		"""

rule multiqc:
	input:
		fastp = expand("fastp_out/json/{sample}_AT_QT_fastp.json", sample=WITHST38),
		quast = expand("quast_out/{sample}", sample=WITHST38)
	output:
		fastp = "multiqc_out/multiqc_fastp.html",
		quast = "multiqc_out/multiqc_quast.html",
		fastp_data = directory("multiqc_out/multiqc_fastp_data"),
		quast_data = directory("multiqc_out/multiqc_quast_data")
	conda:
		"envs/multiqc.yaml"
	log:
		fastp = "logs/multiqc_fastp.log",
		quast = "logs/multiqc_quast.log"
	shell:
		"""
		OUTPUT={output.fastp}
		DIR=${{OUTPUT%/*}}
		NAME=${{OUTPUT##*/}}
		INPUT=$(for i in {input.fastp}; do echo ${{i%/*}}; done | sort | uniq)
		multiqc -n ${{NAME}} -o ${{DIR}} ${{INPUT}} 2>&1>{log.fastp}
		OUTPUT={output.quast}
		DIR=${{OUTPUT%/*}}
		NAME=${{OUTPUT##*/}}
		INPUT=$(for i in {input.quast}; do echo ${{i%/*}}; done | sort | uniq)
		multiqc -n ${{NAME}} -o ${{DIR}} ${{INPUT}} 2>&1>{log.quast}
		"""

rule abricate_summary:
	input:
		ncbi = expand("abricate_out/ncbi/{sample}_ncbi.tsv", sample=WITHST38),
		vfdb = expand("abricate_out/vfdb/{sample}_vfdb.tsv", sample=WITHST38),
		ecoh = expand("abricate_out/ecoh/{sample}_ecoh.tsv", sample=WITHST38),
		plasmidfinder = expand("abricate_out/plasmidfinder/{sample}_plasmidfinder.tsv", sample=WITHST38)
	output:
		ncbi = "abricate_summary/summary_ncbi.tsv",
		vfdb = "abricate_summary/summary_vfdb.tsv",
		ecoh = "abricate_summary/summary_ecoh.tsv",
		plasmidfinder = "abricate_summary/summary_plasmidfinder.tsv"
	log:
		"logs/abricate_summary.log"
	shell:
		"""
		mkdir -p abricate_summary
		abricate --summary {input.ncbi} > {output.ncbi} 2>>{log}
		abricate --summary {input.vfdb} > {output.vfdb} 2>>{log}
		abricate --summary {input.ecoh} > {output.ecoh} 2>>{log}
		abricate --summary {input.plasmidfinder} > {output.plasmidfinder} 2>>{log}
		"""

rule snippy:
	input:
		fw = "trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "trimmed_illumina/{sample}_2_AT_QT.fastq.gz",
		ATCC25922 = config["snippy"]["ATCC25922_ST69"]
	output:
		ATCC25922 = directory("snippy_out/{sample}")
	conda:
		"envs/snippy.yaml"
	params:
		general = config["snippy"]["general"]
	log:
		"logs/snippy/{sample}.log"
	threads: 8
	shell:
		"""
		snippy {params.general} --cpus {threads} --outdir {output.ATCC25922} --ref {input.ATCC25922} --pe1 {input.fw} --pe2 {input.rv} 2>{log}
		"""

### SNP COMPARISONS

rule snippycore_withST38:
	input:
		data = expand("snippy_out/{sample}", sample=WITHST38),
		ref = "references/ATCC25922.gbk"
	output:
		full = "snippy-core_withST38_out/core.full.aln",
		snps = "snippy-core_withST38_out/core.aln"
	conda:
		"envs/snippy.yaml"
	params:
		outdir = "snippy-core_withST38_out"
	log:
		"logs/snippycore.log"
	shell:
		"""
		mkdir -p {params.outdir}
		snippy-core --ref {input.ref} {input.data} 2>&1>{log}
		mv core.aln core.full.aln core.tab core.vcf core.txt core.ref.fa {params.outdir}
		"""

rule snippycore_withoutST38:
	input:
		data = expand("snippy_out/{sample}", sample=WITHOUTST38),
		ref = "references/ATCC25922.gbk",
		dummy = "snippy-core_withST38_out/core.full.aln"
	output:
		full = "snippy-core_withoutST38_out/core.full.aln",
		snps = "snippy-core_withoutST38_out/core.aln"
	conda:
		"envs/snippy.yaml"
	params:
		outdir = "snippy-core_withoutST38_out"
	log:
		"logs/snippycore_withoutST38.log"
	shell:
		"""
		mkdir -p {params.outdir}
		snippy-core --ref {input.ref} {input.data} 2>&1>{log}
		mv core.aln core.full.aln core.tab core.vcf core.txt core.ref.fa {params.outdir}
		"""

rule iqtree_withST38:
	input:
		core = "snippy-core_withST38_out/core.aln",
		fullcore = "snippy-core_withST38_out/core.full.aln"
	output:
		directory("iqtree_withST38_out")
	conda:
		"envs/iqtree_snp-sites.yaml"
	params:
		prefix = config["iqtree"]["prefix"]
	log:
		"logs/iqtree_withST38.log"
	threads: 16
	shell:
		"""
		mkdir -p {output} && cd {output}
		iqtree -fconst $(snp-sites -C ../{input.fullcore}) -nt AUTO -pre {params.prefix} -s ../{input.core} 2>&1>../{log}
		if [ -f {params.prefix}.treefile ]; then echo "{output}/{params.prefix}.treefile exists"; else exit 1; fi
		"""

rule iqtree_withoutST38:
	input:
		core = "snippy-core_withoutST38_out/core.aln",
		fullcore = "snippy-core_withoutST38_out/core.full.aln",
		dummy = "iqtree_withST38_out"
	output:
		directory("iqtree_withoutST38_out")
	conda:
		"envs/iqtree_snp-sites.yaml"
	params:
		prefix = config["iqtree"]["prefix"]
	log:
		"logs/iqtree_withoutST38.log"
	threads: 16
	shell:
		"""
		mkdir -p {output} && cd {output}
		iqtree -fconst $(snp-sites -C ../{input.fullcore}) -nt AUTO -pre {params.prefix} -s ../{input.core} 2>&1>../{log}
		if [ -f {params.prefix}.treefile ]; then echo "{output}/{params.prefix}.treefile exists"; else exit 1; fi
		"""

rule clonalframeml_withST38:
	input:
		tree = "iqtree_withST38_out",
		aln = "snippy-core_withST38_out/core.full.aln"
	output:
		directory("clonalframeml_withST38_out")
	conda:
		"envs/clonalframeml.yaml"
	params:
		prefix = "COMBAT_clonalframeml",
		iqtreeprefix = config["iqtree"]["prefix"]
	log:
		"logs/clonalframeml_withST38.log"
	threads: 16
	shell:
		"""
		mkdir -p {output} && cd {output}
		ClonalFrameML ../{input.tree}/{params.iqtreeprefix}.treefile ../{input.aln} {params.prefix} -show_progress true 2>&1>../{log}
		if [ -f {params.prefix}.labelled_tree.newick ]; then echo "CFML output exists"; else exit 1; fi
		"""

rule clonalframeml_withoutST38:
	input:
		tree = "iqtree_withoutST38_out",
		aln = "snippy-core_withoutST38_out/core.full.aln",
		dummy = "clonalframeml_withST38_out"
	output:
		directory("clonalframeml_withoutST38_out")
	conda:
		"envs/clonalframeml.yaml"
	params:
		prefix = "COMBAT_clonalframeml",
		iqtreeprefix = config["iqtree"]["prefix"]
	log:
		"logs/clonalframeml_withST38.log"
	threads: 16
	shell:
		"""
		mkdir -p {output} && cd {output}
		ClonalFrameML ../{input.tree}/{params.iqtreeprefix}.treefile ../{input.aln} {params.prefix} -show_progress true 2>&1>../{log}
		if [ -f {params.prefix}.labelled_tree.newick ]; then echo "CFML output exists"; else exit 1; fi
		"""

rule maskrc_withST38:
	input:
		aln = "snippy-core_withST38_out/core.full.aln",
		cfml = "clonalframeml_withST38_out"
	output:
		"masked_withST38.aln"
	conda:
		"envs/maskrc.yaml"
	params:
		prefix = "COMBAT_clonalframeml"
	log:
		"logs/maskrc.log"
	shell:
		"""
		bash scripts/download_maskrc.sh
		cd clonalframeml_withST38_out
		python3 ../maskrc-svg.py --aln ../{input.aln} --out ../{output} {params.prefix} 2>&1>../{log}
		"""

rule maskrc_withoutST38:
	input:
		aln = "snippy-core_withoutST38_out/core.full.aln",
		cfml = "clonalframeml_withoutST38_out",
		dummy = "masked_withST38.aln"
	output:
		"masked_withoutST38.aln"
	conda:
		"envs/maskrc.yaml"
	params:
		prefix = "COMBAT_clonalframeml"
	log:
		"logs/maskrc.log"
	shell:
		"""
		bash scripts/download_maskrc.sh
		cd clonalframeml_withoutST38_out
		python3 ../maskrc-svg.py --aln ../{input.aln} --out ../{output} {params.prefix} 2>&1>../{log}
		"""

rule snp_dists:
	input:
		ST38 = "masked_withST38.aln",
		noST38 = "masked_withoutST38.aln"
	output:
		snpmatstandard_withST38 = "snp_comparison/snp_distances_standard_withST38.tsv",
		snpmatnogaps_withST38 = "snp_comparison/snp_distances_no_gaps_withST38.tsv",
		snpmatstandard_withoutST38 = "snp_comparison/snp_distances_standard_withoutST38.tsv",
		snpmatnogaps_withoutST38 = "snp_comparison/snp_distances_no_gaps_withoutST38.tsv"
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/snp_dists.log"
	shell:
		"""
		snp-dists -m {input.ST38} > {output.snpmatstandard_withST38}
		snp-sites -cb {input.ST38} | snp-dists -m /dev/stdin > {output.snpmatnogaps_withST38}
		snp-dists -m {input.noST38} > {output.snpmatstandard_withoutST38}
		snp-sites -cb {input.noST38} | snp-dists -m /dev/stdin > {output.snpmatnogaps_withoutST38}
		"""

rule alnlengths:
	input:
		ST38 = "masked_withST38.aln",
		noST38 = "masked_withoutST38.aln"
	output:
		alnlengths_withST38 = "snp_comparison/alnlengths_withST38.tsv",
		alnlengths_withoutST38 = "snp_comparison/alnlengths_withoutST38.tsv"
	log:
		"logs/alnlengths.log"
	shell:
		"""
		bash scripts/download_snp-dists-alnlengths.sh
		scripts/snp-dists-alnlengths -l {input.ST38} > {output.alnlengths_withST38} 2>{log}
		scripts/snp-dists-alnlengths -l {input.noST38} > {output.alnlengths_withoutST38} 2>>{log}
		"""

rule snp_comparisons:
	input:
		ST38 = "masked_withST38.aln",
		noST38 = "masked_withoutST38.aln",
		snpmatstandard_withST38 = "snp_comparison/snp_distances_standard_withST38.tsv",
		snpmatnogaps_withST38 = "snp_comparison/snp_distances_no_gaps_withST38.tsv",
		snpmatstandard_withoutST38 = "snp_comparison/snp_distances_standard_withoutST38.tsv",
		snpmatnogaps_withoutST38 = "snp_comparison/snp_distances_no_gaps_withoutST38.tsv",
		metadata = "COMBAT_metadata.tsv",
		list_withST38 = "lists/list_isolates_withST38.txt",
		list_withoutST38 = "lists/list_isolates_withoutST38.txt",
		alnlengths_withST38 = "snp_comparison/alnlengths_withST38.tsv",
		alnlengths_withoutST38 = "snp_comparison/alnlengths_withoutST38.tsv"
	output:
		final_withST38 = "snp_comparison/snp_comparisons_withST38.tsv",
		final_withoutST38 = "snp_comparison/snp_comparisons_withoutST38.tsv"
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/snp_comparisons.log"	
	shell:
		"""
		python3 scripts/snp_comparisons.py {input.list_withST38} {input.snpmatstandard_withST38} {input.snpmatnogaps_withST38} {input.alnlengths_withST38} {output.final_withST38} 2>&1>{log}
		python3 scripts/snp_comparisons.py {input.list_withoutST38} {input.snpmatstandard_withoutST38} {input.snpmatnogaps_withoutST38} {input.alnlengths_withoutST38} {output.final_withoutST38} 2>&1>{log}
		"""

rule plot_snp_comparisons:
	input:
		withST38 = "snp_comparison/snp_comparisons_withST38.tsv",
		withoutST38 = "snp_comparison/snp_comparisons_withoutST38.tsv"
	output:
		data_withST38 = "snp_comparison/input_plot_SNP_threshold_withST38.tsv",
		data_withoutST38 = "snp_comparison/input_plot_SNP_threshold_withoutST38.tsv",
		lineplot_withST38 = "snp_comparison/snp_comparisons_thresholds_lineplot_withST38.pdf",
		histo_withST38 = "snp_comparison/snp_comparison_histogram_corrected_50_withST38.pdf",
		lineplot_withoutST38 = "snp_comparison/snp_comparisons_thresholds_lineplot_withoutST38.pdf",
		histo_withoutST38 = "snp_comparison/snp_comparison_histogram_corrected_50_withoutST38.pdf"
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/plot_snp_comparisons"
	shell:
		"""
		python3 scripts/prepare_input_plot_SNP_threshold.py {input.withST38} > {output.data_withST38} 2>{log}
		python3 scripts/prepare_input_plot_SNP_threshold.py {input.withoutST38} > {output.data_withoutST38} 2>>{log}
		Rscript scripts/plot_SNP_threshold_histogram.R 2>&1>>{log}
		Rscript scripts/plot_SNP_threshold_lineplot.R 2>&1>>{log}
		"""

rule print_travelers_withST38:
	input:
		lineplot = "snp_comparison/snp_comparisons_thresholds_lineplot_withST38.pdf",
		snpcomparisons = "snp_comparison/snp_comparisons_withST38.tsv"
	output:
		"traveler_persistence_types_withST38.tsv"
	params:
		threshold_verylikely  =  config["print_travelers_withST38"]["threshold_verylikely"],
		threshold_likely =  config["print_travelers_withST38"]["threshold_likely"]
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/print_travelers_withST38.log"
	shell:
		"""
		python3 scripts/print_travelers.py {params.threshold_verylikely} {params.threshold_likely} {input.snpcomparisons} > {output} 2>{log}
		"""

rule print_travelers_withoutST38:
	input:
		lineplot = "snp_comparison/snp_comparisons_thresholds_lineplot_withoutST38.pdf",
		snpcomparisons = "snp_comparison/snp_comparisons_withoutST38.tsv"
	output:
		"traveler_persistence_types_withoutST38.tsv"
	params:
		threshold_verylikely  =  config["print_travelers_withoutST38"]["threshold_verylikely"],
		threshold_likely =  config["print_travelers_withoutST38"]["threshold_likely"]
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/print_travelers_withoutST38.log"
	shell:
		"""
		python3 scripts/print_travelers.py {params.threshold_verylikely} {params.threshold_likely} {input.snpcomparisons} > {output} 2>{log}
		"""


### NANOPORE DATA

rule filtlong:
	input:
		nanopore = "raw_nanopore/{sample}.fastq.gz",
		fw = "trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "trimmed_illumina/{sample}_2_AT_QT.fastq.gz"
	output:
		"trimmed_nanopore/{sample}.fastq.gz"
	conda:
		"envs/filtlong.yaml"
	params:
		target_bases = config["filtlong"]["target_bases"],
		keep_percent = config["filtlong"]["keep_percent"]
	log:
		"logs/filtlong/{sample}.log"
	shell:
		"""
		filtlong --target_bases {params.target_bases} --illumina_1 {input.fw} --illumina_2 {input.rv} --trim {input.nanopore} | gzip > {output} 2>{log}
		"""

rule fastqc:
	input:
		nanopore = "raw_nanopore/{sample}.fastq.gz"
	output:
		directory("fastqc_out/{sample}")
	conda:
		"envs/fastqc.yaml"
	log:
		"logs/fastqc/{sample}.log"
	threads: 6
	shell:
		"""
		mkdir -p {output}
		fastqc -t {threads} --outdir {output} {input} 2>&1>{log}
		"""

rule unicycler:
	input:
		nanopore = "trimmed_nanopore/{sample}.fastq.gz",
		fw = "trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "trimmed_illumina/{sample}_2_AT_QT.fastq.gz"
	output:
		"unicycler_out/{sample}/assembly.fasta"
	conda:
		"envs/unicycler.yaml"
	params:
		outdir = "unicycler_out/{sample}"
	log:
		"logs/unicycler/{sample}.log"
	threads: 6
	shell:
		"""
		unicycler -1 {input.fw} -2 {input.rv} --long {input.nanopore} -o {params.outdir} --threads {threads} 2>&1>{log}
		"""

rule amrfinder_nanopore:
	input:
		assembly = "unicycler_out/{sample}/assembly.fasta"
	output:
		"amrfinder_nanopore_out/{sample}.tsv"
	conda:
		"envs/amrfinder.yaml"
	params:
		organism = config["amrfinder"]["organism"]
	log:
		"logs/amrfinder_nanopore/{sample}.log"
	threads: 16
	shell:
		"""
		amrfinder -u
		amrfinder --threads {threads} --nucleotide {input.assembly} --organism {params.organism} --output {output} 2>&1>{log}
		"""

rule get_ESBL_contigs:
	input:
		assembly = "unicycler_out/{sample}/assembly.fasta",
		amrfinder = "amrfinder_nanopore_out/{sample}.tsv"		
	output:
		ESBLcontig = "ESBL_contigs/{sample}.fasta"
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/get_ESBL_contigs/{sample}.log"
	shell:
		"""
		bash scripts/get_ESBL_contigs.sh {input.amrfinder} 2>&1>{log}
		"""

rule prokka_ESBL_contigs:
	input:
		"ESBL_contigs/{sample}.fasta"
	output:
		directory("ESBL_contigs_annotations/{sample}")
	conda:
		"envs/prokka.yaml"
	params:
		general = config["prokka"]["general"],
		kingdom = config["prokka"]["kingdom"],
		genus = config["prokka"]["genus"],
		species = config["prokka"]["species"],
		prefix = "{sample}"
	log:
		"logs/prokka_ESBL_contigs/{sample}.log"
	threads: 8
	shell:
		"""
		prokka {params.general} --force --outdir {output} --genus {params.genus} --species {params.species} --kingdom {params.kingdom} --cpus {threads} --prefix {params.prefix} {input} 2>&1>{log}
		if [ -f {output}/*.ffn ]; then echo "{output} exists"; else exit 1; fi
		"""

### MATCHED CONTROLS

rule fastp_controls:
	input:
		fw = "controls/raw_illumina/{sample}_1.fastq.gz",
		rv = "controls/raw_illumina/{sample}_2.fastq.gz"
	output:
		fw = "controls/trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "controls/trimmed_illumina/{sample}_2_AT_QT.fastq.gz",
		html = "controls/fastp_out/html/{sample}_AT_QT_fastp.html",
		json = "controls/fastp_out/json/{sample}_AT_QT_fastp.json"
	conda:
		"envs/fastp.yaml"
	params:
		compression_level = config["fastp"]["compression_level"],
		general = config["fastp"]["general"]
	log:
		"logs/fastp/{sample}.log"
	threads: 8
	shell:
		"""
		fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}
		"""

rule kraken2_controls:
	input:
		fw = "controls/trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "controls/trimmed_illumina/{sample}_2_AT_QT.fastq.gz"
	output:
		report = "controls/kraken_out/{sample}_kraken2_report.txt"
	conda:
		"envs/kraken.yaml"
	params:
		general = config["kraken"]["general"],
		db = config["kraken"]["db"]
	log:
		"logs/kraken2/{sample}.log"
	threads: 8
	shell:
		"""
		kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}
		"""

rule shovill_controls:
	input:
		fw = "controls/trimmed_illumina/{sample}_1_AT_QT.fastq.gz",
		rv = "controls/trimmed_illumina/{sample}_2_AT_QT.fastq.gz"
	output:
		assembly = "controls/genomes/{sample}.fasta",
		shovill = directory("controls/shovill_out/{sample}")
	params:
		minlen = config["shovill"]["minlen"],
		ram = config["shovill"]["ram"],
		depth = config["shovill"]["depth"],
		assembler = config["shovill"]["assembler"],
		tmpdir = config["shovill"]["tmpdir"]
	conda:
		"envs/shovill.yaml"
	log:
		"logs/shovill/{sample}.log"
	threads: 16
	shell:
		"""
		shovill --assembler {params.assembler} --outdir {output.shovill} --tmp {params.tmpdir} --depth {params.depth} --cpus {threads} --ram {params.ram} --minlen {params.minlen} --R1 {input.fw} --R2 {input.rv} 2>&1>{log}
		cp {output.shovill}/contigs.fa {output.assembly}
		"""

rule quast_controls:
	input:
		"controls/genomes/{sample}.fasta"
	output:
		directory("controls/quast_out/{sample}")
	conda:
		"envs/quast.yaml"
	log:
		"logs/quast/{sample}.log"
	threads: 8
	shell:
		"""
		quast --threads {threads} -o {output} {input} 2>&1>{log}
		"""

rule mlst_controls:
	input:
		"controls/genomes/{sample}.fasta"
	output:
		"controls/mlst/{sample}.tsv"
	conda:
		"envs/mlst.yaml"
	log:
		"logs/mlst/{sample}.log"
	shell:
		"""
		mlst {input} > {output} 2>>{log}
		"""

rule ezclermont_controls:
	input:
		"controls/genomes/{sample}.fasta"
	output:
		"controls/ezclermont_out/{sample}.tsv"
	conda:
		"envs/ezclermont.yaml"
	log:
		"logs/ezclermont/{sample}.log"
	shell:
		"""
		set +e
		ezclermont {input} 1>{output} 2>{log} || true
		"""

rule multiqc_controls:
	input:
		fastp = expand("controls/fastp_out/json/{sample}_AT_QT_fastp.json", sample=CONTROLSWITHOUTST38),
		quast = expand("controls/quast_out/{sample}", sample=CONTROLSWITHOUTST38)
	output:
		fastp = "controls/multiqc_out/multiqc_fastp.html",
		quast = "controls/multiqc_out/multiqc_quast.html",
		fastp_data = directory("controls/multiqc_out/multiqc_fastp_data"),
		quast_data = directory("controls/multiqc_out/multiqc_quast_data")
	log:
		fastp = "logs/multiqc_fastp.log",
		quast = "logs/multiqc_quast.log"
	shell:
		"""
		OUTPUT={output.fastp}
		DIR=${{OUTPUT%/*}}
		NAME=${{OUTPUT##*/}}
		INPUT=$(for i in {input.fastp}; do echo ${{i%/*}}; done | sort | uniq)
		multiqc -n ${{NAME}} -o ${{DIR}} ${{INPUT}} 2>&1>{log.fastp}
		OUTPUT={output.quast}
		DIR=${{OUTPUT%/*}}
		NAME=${{OUTPUT##*/}}
		INPUT=$(for i in {input.quast}; do echo ${{i%/*}}; done | sort | uniq)
		multiqc -n ${{NAME}} -o ${{DIR}} ${{INPUT}} 2>&1>{log.quast}
		"""

rule ezclermont_summary_controls:
	input:
		expand("controls/ezclermont_out/{sample}.tsv", sample=CONTROLSWITHOUTST38)
	output:
		"ezclermont_summary_controls.tsv"
	log:
		"logs/ezclermont_summary_controls.log"
	shell:
		"""
		cat {input} 1> {output} 2>{log}
		"""

rule ezclermont_summary_controls_withST38:
	input:
		expand("controls/ezclermont_out/{sample}.tsv", sample=CONTROLSWITHST38)
	output:
		"ezclermont_summary_controls_withST38.tsv"
	log:
		"logs/ezclermont_summary_controls_withST38.log"
	shell:
		"""
		cat {input} 1> {output} 2>{log}
		"""


rule comparison_phylogroups_withST38:
	input:
		snpcomparisons = "snp_comparison/snp_comparisons_withST38.tsv",
		metadata = "COMBAT_metadata.tsv",
		metadata_controls = "COMBAT_metadata_controls.tsv",
		phylo_summary = "ezclermont_summary.tsv",
		phylo_summary_controls = "ezclermont_summary_controls_withST38.tsv",
		travelers = "traveler_persistence_types_withST38.tsv"
	output:
		plotdata = "phylogroup_comparison_plotdata_withST38.tsv",
		final = "phylogroup_comparison_withST38.tsv"
	conda:
		"envs/snp_comparisons.yaml"
	params:
		threshold =  config["print_travelers_withST38"]["threshold_likely"]
	log:
		"logs/phylogroup_comparison.log"
	shell:
		"""
		python3 scripts/phylogroup_comparison.py {params.threshold} {input.snpcomparisons} {input.phylo_summary} {input.phylo_summary_controls} {input.metadata_controls} {output.plotdata} {output.final} 2>&1>{log}
		"""

rule comparison_phylogroups_withoutST38:
	input:
		snpcomparisons = "snp_comparison/snp_comparisons_withoutST38.tsv",
		metadata = "COMBAT_metadata.tsv",
		metadata_controls = "COMBAT_metadata_controls.tsv",
		phylo_summary = "ezclermont_summary.tsv",
		phylo_summary_controls = "ezclermont_summary_controls.tsv",
		travelers = "traveler_persistence_types_withoutST38.tsv"
	output:
		plotdata = "phylogroup_comparison_plotdata_withoutST38.tsv",
		final = "phylogroup_comparison_withoutST38.tsv"
	conda:
		"envs/snp_comparisons.yaml"
	params:
		threshold =  config["print_travelers_withoutST38"]["threshold_likely"]
	log:
		"logs/phylogroup_comparison.log"
	shell:
		"""
		python3 scripts/phylogroup_comparison.py {params.threshold} {input.snpcomparisons} {input.phylo_summary} {input.phylo_summary_controls} {input.metadata_controls} {output.plotdata} {output.final} 2>&1>{log}
		"""

rule plot_phylogroup_comparison_withoutST38:
	input:
		"phylogroup_comparison_plotdata_withoutST38.tsv"
	output:
		"phylogroup_comparison_withoutST38.pdf"
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/plot_phylogroup_comparison.log"
	shell:
		"""
		Rscript scripts/plot_phylogroup_comparison.R {input} {output}
		"""

rule plot_phylogroup_comparison_withST38:
	input:
		"phylogroup_comparison_plotdata_withST38.tsv"
	output:
		"phylogroup_comparison_withST38.pdf"
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/plot_phylogroup_comparison.log"
	shell:
		"""
		Rscript scripts/plot_phylogroup_comparison_withST38.R {input} {output}
		"""
