modelDNA_dir = ""
genome_dir = "/gpfs1/chengqiyi_pkuhpc/hebo/hebo_luster1/ref/hg38"
dataset = "batch9"
sample_ids = open(f'sample_id/{dataset}.txt').read().strip().split('\n')
sample_ids = [ele for ele in sample_ids if 'input' not in ele] #remove input
#sample format:
# <platin>_<second_label>_[1|2].fq.gz <-> <platin>_<second_label>


#data_dir = f"data/{dataset}"
#out_dir = f"output/{dataset}"
#log_dir = f"{out_dir}/{dataset}/log"




rule all:
	input:
		trim1 = expand("output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",dataset = dataset, sample_id = sample_ids),
		trim2 = expand("output/{dataset}/trim/{sample_id}_2_val_2.fq.gz",dataset = dataset, sample_id = sample_ids),
		report1 = expand("output/{dataset}/trim/{sample_id}_1_val_1_fastqc.html",dataset = dataset, sample_id = sample_ids),
		report2 = expand("output/{dataset}/trim/{sample_id}_2_val_2_fastqc.html", dataset = dataset, sample_id = sample_ids),
		rmSbfI1 = expand("output/{dataset}/trim/{sample_id}_rmSbfI_1.fq.gz",dataset = dataset, sample_id = sample_ids),
		rmSbfI2 = expand("output/{dataset}/trim/{sample_id}_rmSbfI_2.fq.gz",dataset = dataset, sample_id = sample_ids),
#		unfastq1 = expand("output/{dataset}/align/model/{sample_id}_unmapped.fq.1.gz",dataset = dataset, sample_id = sample_ids),
#		unfastq2 = expand("output/{dataset}/align/model/{sample_id}_unmapped.fq.2.gz",dataset = dataset, sample_id= sample_ids),
#		spikebam = expand("output/{dataset}/align/model/{sample_id}_model.bam", dataset = dataset , sample_id = sample_ids),
#		spikebam1 = expand("output/{dataset}/align/model1/{sample_id}_model.bam", dataset = dataset, sample_id = sample_ids),
		mipfwd = expand("output/{dataset}/align/model/{sample_id}_fwd.mip", dataset=dataset, sample_id = sample_ids),
		miprvs = expand("output/{dataset}/align/model/{sample_id}_rvs.mip", dataset=dataset, sample_id = sample_ids),
#		mipfwd1 = expand("output/{dataset}/align/model1/{sample_id}_fwd.mip", dataset=dataset, sample_id = sample_ids),
#		miprvs1 = expand("output/{dataset}/align/model1/{sample_id}_rvs.mip", dataset=dataset, sample_id = sample_ids),
		genome_bam = expand("output/{dataset}/align/genome/{sample_id}_genome_sorted.bam",dataset = dataset, sample_id = sample_ids),
		genome_bamstat = expand("output/{dataset}/align/genome/{sample_id}_genome_sorted.bam.stat",dataset = dataset, sample_id = sample_ids),
		rmdupbam = expand("output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam", dataset = dataset, sample_id = sample_ids),
		rmdupbamstat = expand("output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam.stat", dataset = dataset, sample_id = sample_ids),
		#cleanbam = expand("output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_clean.bam", dataset = dataset, sample_id = sample_ids),
		covinfo = expand("output/{dataset}/align/genome/{sample_id}_cover_genome_covered.info", dataset = dataset, sample_id = sample_ids),
		bam_pos = expand("output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam", dataset = dataset, sample_id = sample_ids),
		#mip_pos = expand("output/{dataset}/align/genome/{sample_id}_fwd.mip", dataset = dataset, sample_id = sample_ids),
		bam_neg = expand("output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam", dataset = dataset, sample_id = sample_ids),
		#mip_neg = expand("output/{dataset}/align/genome/{sample_id}_rvs.mip", dataset = dataset, sample_id = sample_ids)
		#rmdupbam_model = expand("output/{dataset}/align/model/{sample_id}_model_sorted_rmdup.bam", dataset = dataset, sample_id = sample_ids),
		#bam_pos_model = expand("output/{dataset}/align/model/{sample_id}_model_sorted_rmdup_fwd.bam", dataset = dataset, sample_id = sample_ids),
		#bam_neg_model = expand("output/{dataset}/align/model/{sample_id}_model_sorted_rmdup_rvs.bam", dataset = dataset, sample_id = sample_ids)
#

rule trim_galore:
	input:
		fastq1 = "data/{dataset}/{sample_id}_1.fq.gz",
		fastq2 = "data/{dataset}/{sample_id}_2.fq.gz"
	output:
		fastq1 = "output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_2_val_2.fq.gz"
	threads: 10
	log:
		"output/{dataset}/log/trim/{sample_id}.log"
	shell:
		"""
		trim_galore -j 10 -q 30 --phred33 --length 25 --stringency 3 --paired -o output/{dataset}/trim/ {input.fastq1} {input.fastq2}
		"""

rule fastqc:
	input:
		fastq1 = "output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_2_val_2.fq.gz"
	output:
		report1 = "output/{dataset}/trim/{sample_id}_1_val_1_fastqc.html",
		report2 = "output/{dataset}/trim/{sample_id}_2_val_2_fastqc.html"
	threads: 3
	shell:
		"""
		fastqc -o output/{wildcards.dataset}/trim/ -t {threads} {input.fastq1} {input.fastq2}
		"""
#		
#
#
rule rm_SbfI:
	input:
		report1 = "output/{dataset}/trim/{sample_id}_1_val_1_fastqc.html", #wait fastqc
		report2 = "output/{dataset}/trim/{sample_id}_2_val_2_fastqc.html",
		fastq1 = "output/{dataset}/trim/{sample_id}_1_val_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_2_val_2.fq.gz"
	output:
		fastq1 = "output/{dataset}/trim/{sample_id}_rmSbfI_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_rmSbfI_2.fq.gz"
	threads: 5
	log:
		"output/{dataset}/log/trim/{sample_id}_rmSbfI.log"
	shell:
		"""
		pigz -d -p {threads} {input.fastq1} {input.fastq2}
		python3 filter_reads_withSbfI.py -read1 output/{wildcards.dataset}/trim/{wildcards.sample_id}_1_val_1.fq -read2 output/{wildcards.dataset}/trim/{wildcards.sample_id}_2_val_2.fq -outname_prx output/{wildcards.dataset}/trim/{wildcards.sample_id}
		pigz -p {threads} output/{wildcards.dataset}/trim/{wildcards.sample_id}_1_val_1.fq output/{wildcards.dataset}/trim/{wildcards.sample_id}_2_val_2.fq
		"""


rule align_model:
	input:
		fastq1 = "output/{dataset}/trim/{sample_id}_rmSbfI_1.fq.gz",
		fastq2 = "output/{dataset}/trim/{sample_id}_rmSbfI_2.fq.gz"
	output:
		unfastq1 = "output/{dataset}/align/model/{sample_id}_unmapped.fq.1.gz",
		spikebam = "output/{dataset}/align/model/{sample_id}_model.bam",
		unfastq2 = "output/{dataset}/align/model/{sample_id}_unmapped.fq.2.gz"
	threads: 20
	log:
		"output/{dataset}/log/align/model/{sample_id}.log"
	shell:
		"""
		bowtie2  --sensitive --no-discordant --end-to-end -N 1 -p 20 -x ./spikein_new/spikein -1 {input.fastq1} -2 {input.fastq2} --un-conc-gz output/{wildcards.dataset}/align/model/{wildcards.sample_id}_unmapped.fq.gz | samtools view -Sbhu -@ 20 -F 4 - | samtools sort -@ 20 -o {output.spikebam} - > {log} 2>&1
		samtools index -@ 20 {output.spikebam}
		"""

rule mip_modelbam:
	input:
		spikebam = "output/{dataset}/align/model/{sample_id}_model.bam"
	output:
		mipfwd = "output/{dataset}/align/model/{sample_id}_fwd.mip",
		miprvs = "output/{dataset}/align/model/{sample_id}_rvs.mip"
	threads: 10
	log:
		"output/{dataset}/log/align/model/{sample_id}.mip"
	shell:
		"""
		samtools view -b -F 276 -@ {threads} -o output/{wildcards.dataset}/align/model/{wildcards.sample_id}_model_fwd.bam {input.spikebam}
		samtools view -b -F 260 -f 16 -@ {threads} -o output/{wildcards.dataset}/align/model/{wildcards.sample_id}_model_rvs.bam {input.spikebam}
		samtools index -@ {threads} output/{wildcards.dataset}/align/model/{wildcards.sample_id}_model_fwd.bam
		samtools index -@ {threads} output/{wildcards.dataset}/align/model/{wildcards.sample_id}_model_rvs.bam
		samtools mpileup -f ./spikein_new/spikein.fa output/{wildcards.dataset}/align/model/{wildcards.sample_id}_model_fwd.bam -o {output.mipfwd} -B -d 10000000 -q 0 -Q 0
		samtools mpileup -f ./spikein_new/spikein.fa output/{wildcards.dataset}/align/model/{wildcards.sample_id}_model_rvs.bam -o {output.miprvs} -B -d 10000000 -q 0 -Q 0
		"""
#
#rule mip_modelbam1:
#	input:
#		spikebam = "output/{dataset}/align/model1/{sample_id}_model.bam"
#	output:
#		mipfwd = "output/{dataset}/align/model1/{sample_id}_fwd.mip",
#		miprvs = "output/{dataset}/align/model1/{sample_id}_rvs.mip"
#	threads: 20
#	log:
#		"output/{dataset}/log/align/model1/{sample_id}.mip"
#	shell:
#		"""
#		samtools view -b -F 276 -@ {threads} -o output/{wildcards.dataset}/align/model1/{wildcards.sample_id}_model_fwd.bam {input.spikebam}
#		samtools view -b -F 260 -f 16 -@ {threads} -o output/{wildcards.dataset}/align/model1/{wildcards.sample_id}_model_rvs.bam {input.spikebam}
#		samtools index -@ {threads} output/{wildcards.dataset}/align/model1/{wildcards.sample_id}_model_fwd.bam
#		samtools index -@ {threads} output/{wildcards.dataset}/align/model1/{wildcards.sample_id}_model_rvs.bam
#		samtools mpileup -f ./spikein_new/spikein.fa output/{wildcards.dataset}/align/model1/{wildcards.sample_id}_model_fwd.bam -o {output.mipfwd} -B -d 10000000 -q 0 -Q 0
#		samtools mpileup -f ./spikein_new/spikein.fa output/{wildcards.dataset}/align/model1/{wildcards.sample_id}_model_rvs.bam -o {output.miprvs} -B -d 10000000 -q 0 -Q 0
#		"""


rule align_genome:
	input:
		unfastq1 = "output/{dataset}/align/model/{sample_id}_unmapped.fq.1.gz",
		unfastq2 = "output/{dataset}/align/model/{sample_id}_unmapped.fq.2.gz"
	output:
		bam = "output/{dataset}/align/genome/{sample_id}_genome_sorted.bam",
	threads: 20
	log:
		"output/{dataset}/log/align/genome/{sample_id}.log"
	shell:
		"""
		bowtie2 --sensitive -N 1 --no-discordant --end-to-end -p 20 -x /gpfs1/chengqiyi_pkuhpc/MB_project/hg38/hg38_only_chromsomes -1 {input.unfastq1} -2 {input.unfastq2} | samtools view -Sbhu -@ 20 -F 780 -f 2 -q 0 - | samtools sort -@ 20 -o {output.bam} - > {log} 2>&1 
		samtools index -@ 20 {output.bam}
		"""


rule rmduplicate:
	input:
		bam = "output/{dataset}/align/genome/{sample_id}_genome_sorted.bam"
	output:
		rmdupbam = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam",
	threads: 5
	log:
		"output/{dataset}/log/align/genome/{sample_id}_rmdup.log"
	shell:
		"""
		java -jar picard.jar MarkDuplicates -I {input.bam} -O {output.rmdupbam} -REMOVE_DUPLICATES true -METRICS_FILE output/{wildcards.dataset}/align/genome/{wildcards.sample_id}_rmdup.txt -CREATE_INDEX true -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
		"""

rule stat:
	input:
		rmdupbam = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam",
		bam = "output/{dataset}/align/genome/{sample_id}_genome_sorted.bam"
	output:
		bamstat = "output/{dataset}/align/genome/{sample_id}_genome_sorted.bam.stat",
		rmdupbamstat = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam.stat"
	threads: 10
	shell:
		"""
		samtools stats -@ 20 {input.bam} > {output.bamstat} &
		samtools stats -@ 20 {input.rmdupbam} > {output.rmdupbamstat} &
		wait
		"""

rule rmdup_coverage:
	input:
		rmdupbam = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam"
	output:
		covinfo = "output/{dataset}/align/genome/{sample_id}_cover_genome_covered.info"
	threads: 5
	log:
		"output/{dataset}/log/align/genome/{sample_id}_cov.log"
	shell:
		"""
		python3 cal_covered_genome.py -file1 {input.rmdupbam} -outname_prx output/{wildcards.dataset}/align/genome/{wildcards.sample_id}_cover_genome
		rm -rf output/{wildcards.dataset}/align/genome/{wildcards.sample_id}_cover_genome.depth
		"""



rule divide_by_strand:
	input:
		rmdupbam = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam"
	output:
		bam_pos = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam",
		bam_neg = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam"
	threads: 20
	log:
		"output/{dataset}/log/align/genome/{sample_id}_devide.log"
	shell:
		"""
		samtools view -b -F 276 -@ 20 -o {output.bam_pos} {input.rmdupbam}
		samtools view -b -F 260 -f 16 -@ 20 -o {output.bam_neg} {input.rmdupbam}
		samtools index -@ 20 {output.bam_pos}
		samtools index -@ 20 {output.bam_neg}
		"""
#rm -rf output/{wildcards.dataset}/align/genome/{wildcards.sample_id}_genome_sorted.bam

rule bam_mip:
	input:
		bam_pos = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam",
		bam_neg = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam"
	output:
		mip_pos = "output/{dataset}/align/genome/{sample_id}_fwd.mip",
		mip_neg = "output/{dataset}/align/genome/{sample_id}_rvs.mip"
	threads: 10
	log:
		"output/{dataset}/log/align/genome/{sample_id}_mip.log"
	shell:
		"""
		samtools mpileup -f hg38/hg38_only_chromsomes.fa {input.bam_pos} -o {output.mip_pos} -B -d 10000000 -q 0 -Q 0 &
		samtools mpileup -f hg38/hg38_only_chromsomes.fa {input.bam_neg} -o {output.mip_neg} -B -d 10000000 -q 0 -Q 0 &
		wait
		"""
		
		

