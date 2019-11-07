configfile: "config.yaml"

#IDS, = glob_wildcards("raw_reads/{id}_1.fastq.gz")
IDS = [ "COMB0108", "COMB0109" ]

rule all:
	input:
		expand("kraken_out/{sample}_kraken2_report.txt", sample=IDS),
		expand("quast_out/{sample}", sample=IDS),
		expand("abricate_out/ncbi/{sample}_ncbi.tsv", sample=IDS),
		expand("mlst/{sample}.tsv", sample=IDS),
		expand("prokka_out/{sample}", sample=IDS),
		expand("ezclermont_out/{sample}.tsv", sample=IDS),
		"multiqc_out/multiqc_fastp.html",
		"summary/summary.csv",
		"abricate_summary/summary_ncbi.tsv",
		expand("snippy_out/{sample}", sample=IDS),
		"masked.aln"

rule fastp:
	input:
		fw = "raw_reads/{sample}_1.fastq.gz",
		rv = "raw_reads/{sample}_2.fastq.gz"
	output:
		fw = "AT_QT_reads/{sample}_1_AT_QT.fastq.gz",
		rv = "AT_QT_reads/{sample}_2_AT_QT.fastq.gz",
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
		"fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}"

rule kraken2:
	input:
		fw = "AT_QT_reads/{sample}_1_AT_QT.fastq.gz",
		rv = "AT_QT_reads/{sample}_2_AT_QT.fastq.gz"
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
		"kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}"

rule spades:
	input:
		fw = "AT_QT_reads/{sample}_1_AT_QT.fastq.gz",
		rv = "AT_QT_reads/{sample}_2_AT_QT.fastq.gz"
	output:
		assembly = "genomes/{sample}.fasta",
		tar_archive = "spades_out/{sample}.tar.gz"
	params:
		intermediate = "{sample}",
		spades = config["spades"]["spades"],
		removesmalls_length = config["spades"]["removesmalls_length"]
	log:
		"logs/spades/{sample}.log"
	threads: 16
	shell:
		"""
		spades.py -1 {input.fw} -2 {input.rv} -o ${{TMPDIR}}/{params.intermediate} {params.spades} --threads {threads} 2>&1>{log}
		perl scripts/removesmalls.pl {params.removesmalls_length} ${{TMPDIR}}/{params.intermediate}/contigs.fasta > {output.assembly}
		tar zcvf {output.tar_archive} ${{TMPDIR}}/{params.intermediate}
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
		quast.py --threads {threads} -o {output} {input} 2>&1>{log}
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
		ezclermont {input} 1>{output} 2>{log}
		"""

rule multiqc:
	input:
		fastp = expand("fastp_out/json/{sample}_AT_QT_fastp.json", sample=IDS),
		quast = expand("quast_out/{sample}", sample=IDS)
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
		kraken = expand("kraken_out/{sample}_kraken2_report.txt", sample=IDS),
		quast = expand("quast_out/{sample}", sample=IDS),
		abricate = expand("abricate_out/ncbi/{sample}_ncbi.tsv", sample=IDS),
		mlst = expand("mlst/{sample}.tsv", sample=IDS),
		prokka = expand("prokka_out/{sample}", sample=IDS),
		ezclermont = expand("ezclermont_out/{sample}.tsv", sample=IDS)
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
		ncbi = expand("abricate_out/ncbi/{sample}_ncbi.tsv", sample=IDS),
		vfdb = expand("abricate_out/vfdb/{sample}_vfdb.tsv", sample=IDS),
		ecoh = expand("abricate_out/ecoh/{sample}_ecoh.tsv", sample=IDS),
		plasmidfinder = expand("abricate_out/plasmidfinder/{sample}_plasmidfinder.tsv", sample=IDS)
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
		fw = "AT_QT_reads/{sample}_1_AT_QT.fastq.gz",
		rv = "AT_QT_reads/{sample}_2_AT_QT.fastq.gz",
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
		data = expand("snippy_out/{sample}", sample=IDS),
		ref = "references/ATCC25922.gbk"
	output:
		"snippy-core_out/core.full.aln"
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
		"snippy-core_out/core.full.aln"
	output:
		"iqtree_out/COMBAT.treefile"
	conda:
		"envs/iqtree_snp-sites.yaml"
	params:
		outdir = "iqtree_out",
		prefix = "COMBAT"
	log:
		"logs/iqtree.log"
	shell:
		"""
		mkdir -p {params.outdir} && cd {params.outdir}
		iqtree -fconst $(snp-sites -C ../{input}) -nt AUTO -pre {params.prefix} -s <(snp-sites -c ../{input})
		"""

rule clonalframeml:
	input:
		tree = "iqtree_out/COMBAT.treefile",
		aln = "snippy-core_out/core.full.aln"
	output:
		"clonalframeml_out"
	conda:
		"envs/clonalframeml.yaml"
	params:
		prefix = "COMBAT_clonalframeml"
	log:
		"logs/clonalframeml.log"
	shell:
		"""
		mkdir -p {output} && cd {output}
		ClonalFrameML {input.tree} {input.aln} {params.prefix} 2>&1>{log}
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
		cd clonalframeml_out
		maskrc-svg.py --aln ../{input.aln} --out {output} {params.prefix} 2>&1>{log}
		"""
