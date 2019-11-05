IDS, = glob_wildcards("raw_reads/{id}_1.fastq.gz")
#IDS = [ "COMB0052", "COMB0053" ]

rule all:
	input:
		expand("kraken_out/{sample}_kraken2_report.txt", sample=IDS),
		expand("quast_out/{sample}", sample=IDS),
		expand("abricate_out/ncbi/{sample}_card.tsv", sample=IDS),
		expand("mlst/{sample}.tsv", sample=IDS),
		expand("prokka_out/{sample}", sample=IDS),
		expand("ezclermont/{sample}.tsv", sample=IDS),
		"multiqc_out/multiqc_fastp.html",
		"summary/summary.csv",
		"abricate_summary/summary_ncbi.tsv",
		expand("snippy_out/{sample}", sample=IDS),
		"snippy-core_out"

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
		path = "/home/vdputten/.conda/envs/my_root/bin",
		compression_level = "9",
		general = "--disable_length_filtering"
	log:
		"logs/fastp/{sample}.log"
	threads: 8
	shell:
		"{params.path}/fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}"

rule kraken2:
	input:
		fw = "AT_QT_reads/{sample}_1_AT_QT.fastq.gz",
		rv = "AT_QT_reads/{sample}_2_AT_QT.fastq.gz"
	output:
		report = "kraken_out/{sample}_kraken2_report.txt"
	conda:
		"envs/kraken.yaml"
	params:
		path = "/home/vdputten/.conda/envs/my_root/bin",
		general = "--output - --fastq-input --gzip-compressed --paired",
		db = "/home/vdputten/kraken_db_nohuman"
	log:
		"logs/kraken2/{sample}.log"
	threads: 8
	shell:
		"{params.path}/kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}"

rule spades:
	input:
		fw = "AT_QT_reads/{sample}_1_AT_QT.fastq.gz",
		rv = "AT_QT_reads/{sample}_2_AT_QT.fastq.gz"
	output:
		assembly = "genomes/{sample}.fasta",
		tar_archive = "spades_out/{sample}.tar.gz"
	params:
		intermediate = "{sample}"
		spades = "--only-assembler --careful",
		removesmalls_length = "500"
	log:
		"logs/spades/{sample}.log"
	threads: 16
	shell:
		"""
		{params.path}/spades.py -1 {input.fw} -2 {input.rv} -o ${{TMPDIR}}/{params.intermediate} {params.spades} --threads {threads} 2>&1>{log}
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
		minid = "90",
		mincov = "60",
		ncbi = "ncbi",
		vfdb = "vfdb",
		plasmidfinder = "plasmidfinder",
		ecoh = "ecoh"
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
		directory("prokka_out/spades/{sample}")
	params:
		general = "--usegenus",
		kingdom = "Bacteria",
		genus = "Escherichia",
		species = "coli",
		outdir_spades = "prokka_out/{sample}",
		prefix = "{sample}"
	log:
		"logs/prokka/{sample}.log"
	threads: 8
	shell:
		"""
		prokka {params.general} --force --outdir {params.outdir} --genus {params.genus} --species {params.species} --kingdom {params.kingdom} --cpus {threads} --prefix {params.prefix} {input} 2>&1>{log}
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
		quast_spades = "logs/multiqc_quast.log"
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
		quastspades = expand("quast_out/spades/{sample}", sample=IDS),
		quastskesa = expand("quast_out/skesa/{sample}", sample=IDS),
		abricate = expand("abricate_out/card/{sample}_card.tsv", sample=IDS),
		mlst = expand("mlst/{sample}.tsv", sample=IDS),
		prokka = expand("prokka_out/spades/{sample}_spades", sample=IDS),
		clermontyper = expand("clermontyper_out/{sample}", sample=IDS)
	output:
		"summary/summary.csv"
	shell:
		"""
		bash scripts/summary_isolates.sh {input.quastspades} > {output}
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
	shell:
		"""
		mkdir -p abricate_summary
		abricate --summary {input.ncbi} > {output.ncbi}
		abricate --summary {input.vfdb} > {output.vfdb}
		abricate --summary {input.ecoh} > {output.ecoh}
		abricate --summary {input.plasmidfinder} > {output.plasmidfinder}
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
	log:
		"logs/snippy/{sample}.log"
	threads: 8
	shell:
		"""
		snippy --force --cpus {threads} --outdir {output} --ref {input.ref} --pe1 {input.fw} --pe2 {input.rv} 2>{log}
		"""

rule snippycore:
	input:
		data = expand("snippy_out/{sample}", sample=IDS),
		ref = "references/ATCC25922.gbk"
	output:
		directory("snippy-core_out")
	conda:
		"envs/snippy.yaml"
	params:
		prefix = "COMBAT"
	shell:
		"""
		snippy-core --ref {input.ref} --prefix {params.prefix} {input.data}
		"""

rule clonalframeml

rule snippy_three_groups


