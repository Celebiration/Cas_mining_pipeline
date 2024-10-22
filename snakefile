import os
from snakemake.logging import logger
#dry-run指令：snakemake --snakefile /Yol_Data/GSA/workflow/snakefile -d . -np
#DAG画图指令：snakemake --dag --snakefile /Yol_Data/GSA/workflow/snakefile -d .|dot -Tsvg > DAG.svg
#指令：nohup snakemake --snakefile /Yol_Data/GSA/workflow/snakefile -p --cores 126 --latency-wait 10 -d . > snakemake.nohup.out 2>&1 &

configfile: "/Yol_Data/GSA/workflow/config.yaml"
contig_suffix=config["contig_suffix"]
run_type=config["run_type"]
run_cctyper=config["run_cctyper"]
l_filter="_"+str(config["filter"])
if l_filter=="_0":
	l_filter=""

if run_type=="GSA":
	sample_name=str(subprocess.check_output(config["scripts"]["get_sample_name"],shell=True).decode('utf-8').strip())
elif run_type=="contig" or run_type=="JGI" or run_type=="SRA":
	sample_name=config["base"] # 需在命令行指定
#print('sample: '+str(sample_name))
logger.info('sample: '+str(sample_name))
def GSA_single_or_pair_input():
	files=os.listdir()
	if any([i in files for i in [sample_name+'_f1.fastq', sample_name+'_f1.fastq.gz', sample_name+'_f1.fq', sample_name+'_f1.fq.gz']]) and \
	any([i in files for i in [sample_name+'_r2.fastq', sample_name+'_r2.fastq.gz', sample_name+'_r2.fq', sample_name+'_r2.fq.gz']]):
		#print("This is a paired library.")
		targets=[sample_name+".contigs"+l_filter+".fa_out/all_summary.csv", "fastqc_p/"+sample_name+"_f1_fastqc.html", "fastqc_p/"+sample_name+"_r2_fastqc.html",sample_name+".contigs"+l_filter+".fa_out/hitted_contigs.uniq.cctyper/hmmer.tab"]
		logger.info("This is a paired library.\ntarget(s): "+str(targets))
		return(targets if run_cctyper else targets[:-1])
	elif sample_name+'.fq' in files or sample_name+'.fastq' in files:
		#print("This is a single library.")
		targets=[sample_name+".contigs"+l_filter+".fa_out/all_summary.csv", "fastqc_s/"+sample_name+"_fastqc.html",sample_name+".contigs"+l_filter+".fa_out/hitted_contigs.uniq.cctyper/hmmer.tab"]
		logger.info("This is a single library.\ntarget(s): "+str(targets))
		return(targets if run_cctyper else targets[:-1])
	else:
		raise Error('no appropriate fastq/gz files.')

if run_type=="GSA":
	rule all:
		input:
			GSA_single_or_pair_input()
elif run_type=="contig":
	if run_cctyper:
		rule all:
			input:
				config['base']+l_filter+".fa_out/all_summary.csv",
				config['base']+l_filter+".fa_out/hitted_contigs.uniq.cctyper/hmmer.tab"
	else:
		rule all:
			input:
				config['base']+l_filter+".fa_out/all_summary.csv"
elif run_type=="JGI":
	if run_cctyper:
		rule all:
			input:
				config['base']+l_filter+"_contigs.fa_out/all_summary.csv",
				config['base']+l_filter+"_contigs.fa_out/hitted_contigs.uniq.cctyper/hmmer.tab"
	else:
		rule all:
			input:
				config['base']+l_filter+"_contigs.fa_out/all_summary.csv"
elif run_type=="SRA":
	if run_cctyper:
		rule all:
			input:
				config['base']+".contigs"+l_filter+".fa_out/all_summary.csv",
				config['base']+".contigs"+l_filter+".fa_out/hitted_contigs.uniq.cctyper/hmmer.tab"
	else:
		rule all:
			input:
				config['base']+".contigs"+l_filter+".fa_out/all_summary.csv"

if run_type=="JGI":
	if os.path.isfile(config['base']+"_proteins.faa"):
		config['protein_file']=config['base']+"_proteins.faa"
	else:
		config['protein_file']="not_empty"


def get_suffix():
	suffix=str(subprocess.check_output(config["scripts"]["get_suffix"],shell=True).decode('utf-8').strip())
	return(suffix)

rule gunzip_file:
	input:
		"{file}.gz"
	output:
		"{file}"
	shell:
		"gunzip {input}"

# rule unzip_fastq_single:
#	 input:
#		 "{sample}."+get_suffix()+".gz"
#	 output:
#		 "{sample}.fastq"
#	 shell:
#		 "gunzip -c {input} > {output};rm {input}"

rule unzip_fastq_pair1:
	input:
		"{sample}_f1."+get_suffix()+".gz"
	output:
		"{sample}_f1.fastq"
	shell:
		"gunzip -c {input} > {output};rm {input}"

rule unzip_fastq_pair2:
	input:
		"{sample}_r2."+get_suffix()+".gz"
	output:
		"{sample}_r2.fastq"
	shell:
		"gunzip -c {input} > {output};rm {input}"

rule fasterq_dump_file:
	input:
		"{sample}.sra"
	output:
		fq1="{sample}_1.fastq",
		fq2="{sample}_2.fastq"
	shell:
		"fasterq-dump {input}"

# rule fastqc_single:
#	 input:
#		 "{sample}.fastq"
#	 output:
#		 "fastqc_s/{sample}_fastqc.html"
#	 threads: 6
#	 shell:
#		 "fastqc -o fastqc_s {input} -t {threads}"

rule fastqc_pair:
	input:
		"{sample}_1.fastq",
		"{sample}_2.fastq"
	output:
		"fastqc_p/{sample}_1_fastqc.html",
		"fastqc_p/{sample}_2_fastqc.html"
	threads: 6
	shell:
		"fastqc -o fastqc_p {input} -t {threads}"

rule fastqc_pair_GSA:
	input:
		"{sample}_f1.fastq",
		"{sample}_r2.fastq"
	output:
		"fastqc_p/{sample}_f1_fastqc.html",
		"fastqc_p/{sample}_r2_fastqc.html"
	threads: 6
	shell:
		"fastqc -o fastqc_p {input} -t {threads}"

# rule fastp_single:
#	 input:
#		 "{sample}.fastq"
#	 output:
#		 "{sample}.fastp_s.fastq"
#	 log:
#		 "log/{sample}_fastp.log.txt"
#	 threads: 16
#	 shell:
#		 "fastp -i {input} -o {output} -w {threads} > {log} 2>&1"

rule fastp_pair:
	input:
		fq1="{sample}_1.fastq",
		fq2="{sample}_2.fastq"
	output:
		fq1="{sample}_1.fastp_p.fastq",
		fq2="{sample}_2.fastp_p.fastq"
	log:
		"log/{sample}_fastp.log.txt"
	threads: 16
	shell:
		"fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} -w {threads} > {log} 2>&1"

rule fastp_pair_GSA:
	input:
		fq1="{sample}_f1.fastq",
		fq2="{sample}_r2.fastq"
	output:
		fq1="{sample}_f1.fastp_p.fastq",
		fq2="{sample}_r2.fastp_p.fastq"
	log:
		"log/{sample}_fastp.log.txt"
	threads: 16
	shell:
		"fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} -w {threads} > {log} 2>&1"

if config["run_fastp"]:
	megahit_ind=".fastp_p"
else:
	megahit_ind=""

rule megahit_paired_GSA:
	input:
		fq1="{sample}_f1"+megahit_ind+".fastq",
		fq2="{sample}_r2"+megahit_ind+".fastq"
	output:
		"{sample}.contigs.fa"
	log:
		"log/{sample}_megahit.log.txt"
	params:
		presets=config["megahit"]["presets"],
		mincount=config["megahit"]["mincount"],
		klist=config["megahit"]["k-list"]
	threads: threads=config["megahit"]["threads"]
	shell:
		"megahit -1 {input.fq1} -2 {input.fq2} --k-list {params.klist} -f --out-prefix {wildcards.sample} --min-contig-len 1000 -t {threads} > {log} 2>&1;mv megahit_out/{output} ./"
		#"megahit -1 {input.fq1} -2 {input.fq2} --presets {params.presets} --min-count {params.mincount} --k-list {params.klist} -f --out-prefix {wildcards.sample} --min-contig-len 1000 -t {params.threads} > {log} 2>&1"

rule megahit_paired:
	input:
		fq1="{sample}_1"+megahit_ind+".fastq",
		fq2="{sample}_2"+megahit_ind+".fastq"
	output:
		"{sample}.contigs.fa"
	log:
		"log/{sample}_megahit.log.txt"
	params:
		presets=config["megahit"]["presets"],
		mincount=config["megahit"]["mincount"],
		klist=config["megahit"]["k-list"]
	threads: config["megahit"]["threads"]
	shell:
		#"megahit -1 {input.fq1} -2 {input.fq2} --k-list {params.klist} -f --out-prefix {wildcards.sample} --min-contig-len 200 -t {threads} --continue > {log} 2>&1;mv megahit_out/{output} ./"
		"megahit -1 {input.fq1} -2 {input.fq2} --presets {params.presets} --min-count {params.mincount} -f --out-prefix {wildcards.sample} --min-contig-len 150 -t {threads} --continue >> {log} 2>&1;mv megahit_out/{output} ./"

# rule megahit_single:
#	 input:
#		 "{sample}.fastp_s.fastq"
#	 output:
#		 "{sample}.contigs.fa"
#	 log:
#		 "log/{sample}_megahit.log.txt"
#	 params:
#		 presets=config["megahit"]["presets"],
#		 mincount=config["megahit"]["mincount"],
#		 klist=config["megahit"]["k-list"],
#		 threads=config["megahit"]["threads"]
#	 threads: 118
#	 shell:
#		 "megahit -r {input} --k-list {params.klist} -f --out-prefix {wildcards.sample} --min-contig-len 1000 -t {params.threads} > {log} 2>&1;mv megahit_out/{output} ./"
		#"megahit -r {input} --presets {params.presets} --min-count {params.mincount} --k-list {params.klist} -f --out-prefix {wildcards.sample} --min-contig-len 1000 -t {params.threads} > {log} 2>&1"

if contig_suffix!="fa":
	rule generalize_contig_suffix:
		input:
			"{contig}."+contig_suffix
		output:
			"{contig}.fa"
		shell:
			"mv {input} {output}"
if l_filter!="":
	rule filter_long_seqs:
		input:
			"{contig}.fa"
		output:
			"{contig}"+l_filter+".fa"
		threads: 8
		shell:
			"seqkit seq -m "+str(config["filter"])+" {input} -o {output} -j {threads}"

protein0=config["protein_file"]
if protein0=="not_empty":
	protein="contig_prodigal"
else:
	protein1=protein0.split("/")[-1]
	protein_suffix=protein1.split(".")[-1]
	if len(protein_suffix)==len(protein1):
		protein=protein1
		protein_suffix=""
	else:
		protein=protein1[:-(len(protein_suffix)+1)]
if protein=="contig_prodigal":
	rule find_genes:
		input:
			"{contig}.fa"
		output:
			"{contig}.fa_out/contig_prodigal.fa"
		params:
			genome_type=config["genome_type"]
		log:
			"log/{contig}_prodigal.log.txt"
		shell:
			"prodigal -i {input} -o {wildcards.contig}.fa_out/contig_prodigal.gbk -m -a {output} -p {params.genome_type} > {log} 2>&1"
else:
	rule mv_protein_file:
		input:
			protein0
		output:
			"{contig}.fa_out/"+protein+".fa"
		shell:
			"cp {input} {output}"

rule find_cas12:
	input:
		gene_file="{contig}.fa_out/"+protein+".fa",
		hmm_file=config["hmm_file"]
	output:
		"{contig}.fa_out/"+protein+"_per_seq.txt"
	log:
		"log/{contig}_hmmsearch.log.txt"
	params:
		hmmsearch_cut_command=config["hmmsearch_cut_command"]
	threads: 16
	shell:
		"hmmsearch -o {wildcards.contig}.fa_out/"+protein+"_hmmout.txt --tblout {output} --noali {params.hmmsearch_cut_command} --cpu {threads} {input.hmm_file} {input.gene_file} > {log} 2>&1"
		#"hmmsearch -o {wildcards.contig}.fa_out/"+protein+"_hmmout.txt -A {wildcards.contig}.fa_out/"+protein+"_ali.txt --tblout {output} --noali {params.hmmsearch_cut_command} --cpu {threads} "+config["hmm_file"]+" {input.gene_file} > {log} 2>&1"

rule process_per_seq:
	input:
		contig="{contig}.fa",
		perseq="{contig}.fa_out/"+protein+"_per_seq.txt"
	output:
		cas12i_tmp="{contig}.fa_out/"+protein+"_cas12i_tmp.txt",
		hitted_contig="{contig}.fa_out/hitted_contigs.fa",
		uniq_hitted_contig="{contig}.fa_out/hitted_contigs.uniq.fa"
	params:
		nametype=protein
	shell:
		config["scripts"]["process_per_seq"]+" -f {input.contig} -i {input.perseq} -r {params.nametype}"

rule find_crispr:
	input:
		"{contig}.fa_out/hitted_contigs.fa"
	output:
		"{contig}.fa_out/hitted_CRISPR_out/result.json"
	log:
		"log/{contig}_find_crispr.log.txt"
	shell:
		"""
			if [ -d {wildcards.contig}.fa_out/hitted_CRISPR_out ];then rm -rf {wildcards.contig}.fa_out/hitted_CRISPR_out;fi
			cd {wildcards.contig}.fa_out
			sed 's/|/_fuck_/g' hitted_contigs.fa > tmp.fa
			sed -i 's/\./_shit_/g' tmp.fa
			seqkit rmdup -n tmp.fa -o tmp_uniq.fa
			rm tmp.fa
			CRISPRCasFinder.pl \
				--in tmp_uniq.fa \
				-out hitted_CRISPR_out \
				-so /home/fengchr/software/CRISPRCasFinder/sel392v2.so \
				-DBcrispr /home/fengchr/software/CRISPRCasFinder/supplementary_files/CRISPR_crisprdb.csv \
				-repeats /home/fengchr/software/CRISPRCasFinder/supplementary_files/Repeat_List.csv \
				-DIRrepeat /home/fengchr/software/CRISPRCasFinder/supplementary_files/repeatDirection.tsv \
				-q > ../{log} 2>&1
			sed -i 's/_fuck_/|/g' hitted_CRISPR_out/result.json
			sed -i 's/_shit_/\./g' hitted_CRISPR_out/result.json
			rm tmp_uniq.fa
			cd ../
		"""

rule final_process:
	input:
		cas12i_tmp="{contig}.fa_out/"+protein+"_cas12i_tmp.txt",
		crispr_out="{contig}.fa_out/hitted_CRISPR_out/result.json"
	output:
		"{contig}.fa_out/all_summary.csv"
	params:
		work_dir="{contig}.fa_out",
		cutoff=config["params"]["cutoff"]
	log:
		"log/{contig}_final_process.txt"
	shell:
		"process_hmmoutput_cas12bij.R {params.work_dir} {input.cas12i_tmp} {params.cutoff} > {log} 2>&1"

rule cctyper:
	input:
		"{contig}.fa_out/hitted_contigs.uniq.fa"
	output:
		"{contig}.fa_out/hitted_contigs.uniq.cctyper/hmmer.tab"
	log:
		"log/{contig}_cctyper.log.txt"
	threads: 4
	shell:
		"""
			if [ -d {wildcards.contig}.fa_out/hitted_contigs.uniq.cctyper ];then rm -rf {wildcards.contig}.fa_out/hitted_contigs.uniq.cctyper;fi
			bash -c 'source ~/.bashrc
			source activate cctyper
			cctyper --prodigal meta {input} {wildcards.contig}.fa_out/hitted_contigs.uniq.cctyper -t {threads} > {log} 2>&1
			conda deactivate'
		"""

# rule cctyper:
#	 input:
#		 "{contig}.fa_out/hitted_contigs.uniq.fa"
#	 output:
#		 "{contig}.fa_out/hitted_contigs.uniq.cctyper/hmmer.tab"
#	 log:
#		 "log/{contig}_cctyper.log.txt"
#	 threads: 4
#	 conda:
#		 config['envs']['cctyper']
#	 shell:
#		 """
#			 bash -c 'source ~/.bashrc
#			 conda activate cctyper
#			 cctyper --prodigal meta {input} hitted_contigs.uniq.cctyper -t {threads} > {log} 2>&1
#			 conda deactivate'
#		 """

# rule process_contig:
#	 input:
#		 expand("megahit_out/{sample}.contigs.fa", sample=sample_name)
#	 output:
#		 "megahit_out/all_summary.csv"
#	 log:
#		 expand("log/{sample}_process_contig.log.txt", sample=sample_name)
#	 shell:
#		 "process_contig.sh -i {input} > {log} 2>&1"