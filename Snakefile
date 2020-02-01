configfile: "config.yaml"

NANOPORE, = glob_wildcards("raw_nanopore/{id}.fastq.gz")
ILLUMINA, = glob_wildcards("raw_illumina/{id}_1.fastq.gz")
CONTROLS, = glob_wildcards("controls/raw_illumina/{id}_1.fastq.gz")

#ILLUMINA = [ "COMB0108", "COMB0109" ]
#NANOPORE = [ "COMB0108" ]

rule all:
	input:
		expand("kraken_out/{sample}_kraken2_report.txt", sample=ILLUMINA),
		expand("quast_out/{sample}", sample=ILLUMINA),
		expand("abricate_out/ncbi/{sample}_ncbi.tsv", sample=ILLUMINA),
		expand("amrfinder_out/{sample}.tsv", sample=ILLUMINA),
		expand("mlst/{sample}.tsv", sample=ILLUMINA),
		expand("prokka_out/{sample}", sample=ILLUMINA),
		expand("ezclermont_out/{sample}.tsv", sample=ILLUMINA),
		"multiqc_out/multiqc_fastp.html",
		"summary/summary.csv",
		"abricate_summary/summary_ncbi.tsv",
		expand("snippy_out/{sample}", sample=ILLUMINA),
		"masked.aln",
		expand("fastqc_out/{sample}", sample=NANOPORE),
		expand("unicycler_out/{sample}/assembly.fasta", sample=NANOPORE),
		"snp_comparisons.tsv",
		"plots/snp_comparisons_500.pdf",
		"traveler_persistence_types.tsv",
		expand("ESBL_contigs/{sample}.fasta", sample=NANOPORE),
		expand("ESBL_contigs_annotations/{sample}", sample=NANOPORE),
		"ezclermont_summary.tsv",
		expand("controls/kraken_out/{sample}_kraken2_report.txt", sample=CONTROLS),
		expand("controls/quast_out/{sample}", sample=CONTROLS),
		expand("controls/mlst/{sample}.tsv", sample=CONTROLS),
		expand("controls/ezclermont_out/{sample}.tsv", sample=CONTROLS),
		"controls/multiqc_out/multiqc_fastp.html",
		"ezclermont_summary_controls.tsv",
		"phylogroup_comparison.tsv"


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
		expand("ezclermont_out/{sample}.tsv", sample=ILLUMINA)
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
		fastp = expand("fastp_out/json/{sample}_AT_QT_fastp.json", sample=ILLUMINA),
		quast = expand("quast_out/{sample}", sample=ILLUMINA)
	output:
		fastp = "multiqc_out/multiqc_fastp.html",
		quast = "multiqc_out/multiqc_quast.html",
		fastp_data = directory("multiqc_out/multiqc_fastp_data"),
		quast_data = directory("multiqc_out/multiqc_quast_data")
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

rule summary:
	input:
		kraken = expand("kraken_out/{sample}_kraken2_report.txt", sample=ILLUMINA),
		quast = expand("quast_out/{sample}", sample=ILLUMINA),
		abricate = expand("abricate_out/ncbi/{sample}_ncbi.tsv", sample=ILLUMINA),
		mlst = expand("mlst/{sample}.tsv", sample=ILLUMINA),
		prokka = expand("prokka_out/{sample}", sample=ILLUMINA),
		ezclermont = expand("ezclermont_out/{sample}.tsv", sample=ILLUMINA)
	output:
		"summary/summary.csv"
	log:
		"logs/summary.log"
	shell:
		"""
		bash scripts/summary_isolates.sh {input.quast} > {output} 2>{log}
		"""

rule abricate_summary:
	input:
		ncbi = expand("abricate_out/ncbi/{sample}_ncbi.tsv", sample=ILLUMINA),
		vfdb = expand("abricate_out/vfdb/{sample}_vfdb.tsv", sample=ILLUMINA),
		ecoh = expand("abricate_out/ecoh/{sample}_ecoh.tsv", sample=ILLUMINA),
		plasmidfinder = expand("abricate_out/plasmidfinder/{sample}_plasmidfinder.tsv", sample=ILLUMINA)
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
		ref = "references/ATCC25922.gbk"
	output:
		directory("snippy_out/{sample}")
	conda:
		"envs/snippy.yaml"
	params:
		general = config["snippy"]["general"]
	log:
		"logs/snippy/{sample}.log"
	threads: 8
	shell:
		"""
		snippy {params.general} --cpus {threads} --outdir {output} --ref {input.ref} --pe1 {input.fw} --pe2 {input.rv} 2>{log}
		"""

rule snippycore:
	input:
		data = expand("snippy_out/{sample}", sample=ILLUMINA),
		ref = "references/ATCC25922.gbk"
	output:
		full = "snippy-core_out/core.full.aln",
		snps = "snippy-core_out/core.aln"
	conda:
		"envs/snippy.yaml"
	params:
		outdir = "snippy-core_out"
	log:
		"logs/snippycore.log"
	shell:
		"""
		mkdir -p {params.outdir}
		snippy-core --ref {input.ref} {input.data} 2>&1>{log}
		mv core.aln core.full.aln core.tab core.vcf core.txt core.ref.fa {params.outdir}
		"""

rule iqtree:
	input:
		core = "snippy-core_out/core.aln",
		fullcore = "snippy-core_out/core.full.aln"
	output:
		directory("iqtree_out")
	conda:
		"envs/iqtree_snp-sites.yaml"
	params:
		prefix = config["iqtree"]["prefix"]
	log:
		"logs/iqtree.log"
	shell:
		"""
		mkdir -p {output} && cd {output}
		iqtree -fconst $(snp-sites -C ../{input.fullcore}) -nt AUTO -pre {params.prefix} -s ../{input.core} 2>&1>../{log}
		if [ -f {params.prefix}.treefile ]; then echo "{output}/{params.prefix}.treefile exists"; else exit 1; fi
		"""

rule clonalframeml:
	input:
		tree = "iqtree_out",
		aln = "snippy-core_out/core.full.aln"
	output:
		directory("clonalframeml_out")
	conda:
		"envs/clonalframeml.yaml"
	params:
		prefix = "COMBAT_clonalframeml",
		iqtreeprefix = config["iqtree"]["prefix"]
	log:
		"logs/clonalframeml.log"
	shell:
		"""
		mkdir -p {output} && cd {output}
		ClonalFrameML ../{input.tree}/{params.iqtreeprefix}.treefile ../{input.aln} {params.prefix} 2>&1>../{log}
		if [ -f {params.prefix}.labelled_tree.newick ]; then echo "CFML output exists"; else exit 1; fi
		"""

rule maskrc:
	input:
		aln = "snippy-core_out/core.full.aln",
		cfml = "clonalframeml_out"
	output:
		"masked.aln"
	conda:
		"envs/maskrc.yaml"
	params:
		prefix = "COMBAT_clonalframeml"
	log:
		"logs/maskrc.log"
	shell:
		"""
		bash scripts/download_maskrc.sh
		cd clonalframeml_out
		python3 ../maskrc-svg.py --aln ../{input.aln} --out ../{output} {params.prefix} 2>&1>../{log}
		"""

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

rule snp_comparisons:
	input:
		"masked.aln",
		"COMBAT_metadata.tsv",
		"list_strains.txt"
	output:
		"snp_comparisons.tsv"
	conda:
		"envs/snp_comparisons.yaml"
	log:
		"logs/snp_comparisons"
	shell:
		"""
		python3 scripts/snp_comparisons.py 2>&1>{log}
		bash scripts/snp_comparisons.sh {output} 2>&1>{log}
		"""

rule plot_snp_comparisons:
	input:
		"snp_comparisons.tsv"
	output:
		"plots/snp_comparisons_50.pdf",
		"plots/snp_comparisons_500.pdf",
		"plots/snp_comparisons_full.pdf"
	conda:
		"envs/plot_snp_comparisons.yaml"
	log:
		"logs/plot_snp_comparisons"
	shell:
		"""
		Rscript scripts/plot_snp_comparisons.R 2>&1>{log}
		"""

rule print_travelers:
	input:
		"plots/snp_comparisons_50.pdf",
		"snp_comparisons.tsv"
	output:
		"traveler_persistence_types.tsv"
	params:
		threshold_verylikely  =  config["print_travelers"]["threshold_verylikely"],
		threshold_likely =  config["print_travelers"]["threshold_likely"]
	log:
		"logs/print_travelers.log"
	shell:
		"""
		bash scripts/print_travelers.py {params.threshold_verylikely} {params.threshold_likely} > {output} 2>{log}
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
		fastp = expand("controls/fastp_out/json/{sample}_AT_QT_fastp.json", sample=CONTROLS),
		quast = expand("controls/quast_out/{sample}", sample=CONTROLS)
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
		expand("controls/ezclermont_out/{sample}.tsv", sample=CONTROLS)
	output:
		"ezclermont_summary_controls.tsv"
	log:
		"logs/ezclermont_summary_controls.log"
	shell:
		"""
		cat {input} 1> {output} 2>{log}
		"""

rule comparison_phylogroups:
	input:
		"snp_comparisons.tsv",
		"COMBAT_metadata.tsv",
		"COMBAT_metadata_controls.tsv",
		"ezclermont_summary.tsv",
		"ezclermont_summary_controls.tsv"
	output:
		"phylogroup_comparison.tsv"
	conda:
		"envs/snp_comparisons.yaml"
	params:
		threshold =  config["print_travelers"]["threshold_likely"]
	log:
		"logs/phylogroup_comparison.log"
	shell:
		"""
		python3 scripts/phylogroup_comparison.py {params.threshold} 2>&1>{log}
		"""
