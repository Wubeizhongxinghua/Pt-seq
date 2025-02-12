'''    
Copyright (c) 2024-02-12 by LiMingyang, YiLab, Peking University.    
     
Author: Li Mingyang (limingyang200101@gmail.com)    
     
Institute: AAIS, Peking University    
     
File Name: PTseq_analysis/xx.py                                                                              
     
Permission is hereby granted, free of charge, to any person obtaining a copy    
of this software and associated documentation files (the "Software"), to deal    
in the Software without restriction, including without limitation the rights    
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell    
copies of the Software, and to permit persons to whom the Software is    
furnished to do so, subject to the following conditions:    
     
The above copyright notice and this permission notice shall be included in all    
copies or substantial portions of the Software.    
''' 

genome_dir = ""
genome_fasta_file = ""
genome_index_basename = ""

datasets = []

sample_id_dict = {}

for dataset in datasets:
#sample format:
# <platin>_<second_label>_[1|2].fq.gz <-> <platin>_<second_label>
    sample_id_dict[dataset] = open(f'sample_id/{dataset}.txt').read().strip().split('\n')


def gen_input():
    required_file = []
    for dataset in datasets:
        for sample_id in sample_id_dict[dataset]:
            required_file.append(f"output/{dataset}/trim/{sample_id}_1_val_1.fq.gz")
    		required_file.append(f"output/{dataset}/trim/{sample_id}_2_val_2.fq.gz")
    		required_file.append(f"output/{dataset}/trim/{sample_id}_1_val_1_fastqc.html")
    		required_file.append(f"output/{dataset}/trim/{sample_id}_2_val_2_fastqc.html")
    		required_file.append(f"output/{dataset}/trim/{sample_id}_rmSbfI_1.fq.gz")
    		required_file.append(f"output/{dataset}/trim/{sample_id}_rmSbfI_2.fq.gz")
    		required_file.append(f"output/{dataset}/align/model/{sample_id}_fwd.mip")
    		required_file.append(f"output/{dataset}/align/model/{sample_id}_rvs.mip")
    		required_file.append(f"output/{dataset}/align/genome/{sample_id}_genome_sorted.bam")
    		required_file.append(f"output/{dataset}/align/genome/{sample_id}_genome_sorted.bam.stat")
    		required_file.append(f"output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam")
    		required_file.append(f"output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam.stat")
    		required_file.append(f"output/{dataset}/align/genome/{sample_id}_cover_genome_covered.info")
    		required_file.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam")
    		required_file.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam")
            required_file.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd_sig.txt")
            required_file.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs_sig.txt")
    return required_file


rule all:
	input:
        gen_input()

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

rule rm_SbfI:
	input:
		report1 = "output/{dataset}/trim/{sample_id}_1_val_1_fastqc.html",
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

rule align_genome:
	input:
		unfastq1 = "output/{dataset}/align/model/{sample_id}_rmSbfI_1.fq.gz",
		unfastq2 = "output/{dataset}/align/model/{sample_id}_rmSbfI_2.fq.gz"
	output:
		bam = "output/{dataset}/align/genome/{sample_id}_genome_sorted.bam",
	threads: 20
	log:
		"output/{dataset}/log/align/genome/{sample_id}.log"
	shell:
		"""
		bowtie2 --sensitive -N 1 --no-discordant --end-to-end -p 20 -x {genome_dir}/{genome_index_basename} -1 {input.unfastq1} -2 {input.unfastq2} | samtools view -Sbhu -@ 20 -F 780 -f 2 -q 0 - | samtools sort -@ 20 -o {output.bam} - > {log} 2>&1 
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


rule ext_bam_all:
	input:
		fwdbam = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam",
	output:
		fwdbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam",
	threads: 3
	shell:
		"""
		sambamba view -t {threads} -f bam --filter "mapping_quality>=40 or (mapping_quality>=1 and mapping_quality<=10)" {input.fwdbam} > {output.fwdbam}
		"""


rule divide_by_strand:
	input:
		rmdupbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam"
	output:
		bam_pos = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam",
		bam_neg = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam"
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

rule to_mip:
	input:
		fwdbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam",
		rvsbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam"
	output:
		fwdmip = temp("output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd.mip"),
		rvsmip = temp("output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs.mip")
	threads: 3
	shell:
		"""
		samtools mpileup -f {genome_dir}/{genome_fasta_file} -o {output.fwdmip} -B -d 10000000 -q 0 -Q 0 {input.fwdbam} &
		samtools mpileup -f {genome_dir}/{genome_fasta_file} -o {output.rvsmip} -B -d 10000000 -q 0 -Q 0 {input.rvsbam} &
		wait
		"""

rule to_sig:
	input:
		fwdmip = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd.mip",
		rvsmip = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs.mip"
	output:
		fwdsig = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd_sig.txt",
		rvssig = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs_sig.txt"
	threads: 3
	shell:
		"""
		python3 tools/0.all_signal.py -i {input.fwdmip} -o {output.fwdsig} -s fwd -p {threads} &
		python3 tools/0.all_signal.py -i {input.rvsmip} -o {output.rvssig} -s rvs -p {threads} &
		wait
		"""



################################# Process Model Sequence ##########################################

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
		bowtie2  --sensitive --no-discordant --end-to-end -N 1 -p 20 -x ./spikein/spikein -1 {input.fastq1} -2 {input.fastq2} --un-conc-gz output/{wildcards.dataset}/align/model/{wildcards.sample_id}_unmapped.fq.gz | samtools view -Sbhu -@ 20 -F 4 - | samtools sort -@ 20 -o {output.spikebam} - > {log} 2>&1
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

