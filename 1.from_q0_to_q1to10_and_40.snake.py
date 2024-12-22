'''
Copyright (c) 2024-03-26 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/from_q0_to_q1_10_40.snake.py

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
'''
from rich import print, pretty
from rich.traceback import install
pretty.install()
install(show_locals=True)

#datasets = ['batch3','batch4','batch5','pnas2016','batch7','batch8','batch9','batch10']
datasets = ['batch11_spike']

dataset_sample = {}
for dataset in datasets:
	sample_ids = open(f'sample_id/{dataset}.txt').read().strip().split('\n')
	#sample_ids = [ele for ele in sample_ids if 'input' not in ele]
	dataset_sample[dataset] = sample_ids
	if dataset in ['pnas2016', 'batch7']:
		dataset_sample[dataset].append('ctrl')
	if dataset in ['batch8']:
		mask_ctrls = [f'ctrl-{i}' for i in [1,2,3]]
		dataset_sample[dataset] = [ele for ele in dataset_sample[dataset] if not any(f in ele for f in mask_ctrls)]
		dataset_sample[dataset].append('Dctrl')
		dataset_sample[dataset].append('Wctrl')

def gen_input():
	all_input = []
	for dataset in datasets:
		for sample_id in dataset_sample[dataset]:
			#all_input.append(f'output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd.mip')
			all_input.append(f'output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd_sig.txt')
			#all_input.append(f'output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs.mip')
			all_input.append(f'output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs_sig.txt')
			all_input.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam")
			all_input.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam")
			all_input.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam")
			all_input.append(f"output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_removed_mip.flag")
		#all_input.append(f"output_q1_10_40_unique/{dataset}/align/genome_down/down_finished.flag")
		#all_input.append(f"output_q1_10_40_unique/{dataset}/align/genome_down/downmip_finished.flag")
	return all_input

rule all:
	input:
		gen_input()



rule ext_bam:
	input:
		fwdbam = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam",
		rvsbam = "output/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam"
	output:
		fwdbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam",
		rvsbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam"
	threads: 3
	shell:
		"""
		sambamba view -t {threads} -f bam --filter "mapping_quality>=40 or (mapping_quality>=1 and mapping_quality<=10)" {input.fwdbam} > {output.fwdbam} &
		sambamba view -t {threads} -f bam --filter "mapping_quality>=40 or (mapping_quality>=1 and mapping_quality<=10)" {input.rvsbam} > {output.rvsbam} &
		wait
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

rule to_mip:
	input:
		fwdbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_fwd.bam",
		rvsbam = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup_rvs.bam"
	output:
		fwdmip = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd.mip",
		rvsmip = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs.mip"
	threads: 3
	shell:
		"""
		samtools mpileup -f hg38/hg38_only_chromsomes.fa -o {output.fwdmip} -B -d 10000000 -q 0 -Q 0 {input.fwdbam} &
		samtools mpileup -f hg38/hg38_only_chromsomes.fa -o {output.rvsmip} -B -d 10000000 -q 0 -Q 0 {input.rvsbam} &
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
rule rm_mip:
	input:
		fwdsig = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_fwd_sig.txt",
		rvssig = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_rvs_sig.txt"
	output:
		flag = "output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_removed_mip.flag"
	threads: 3
	shell:
		"""
		rm -rf output_q1_10_40_unique/{wildcards.dataset}/align/genome/{wildcards.sample_id}_fwd.mip output_q1_10_40_unique/{wildcards.dataset}/align/genome/{wildcards.sample_id}_rvs.mip && touch {output.flag}
		"""


#rule downsample:
#	input:
#		fwdbam = lambda wildcards: expand("output_q1_10_40_unique/{dataset}/align/genome/{sample_id}_genome_sorted_rmdup.bam", dataset=wildcards.dataset, sample_id=dataset_sample[wildcards.dataset])
#	output:
#		downbam = "output_q1_10_40_unique/{dataset}/align/genome_down/down_finished.flag"
#	threads: 20
#	shell:
#		"""
#		thedataset="{wildcards.dataset}"
#		if [[ ${{thedataset}} == "batch7" ]]; then
#			~/Applications/zsh-5.9/bin/zsh tools/3.down_sample.sh output_q1_10_40_unique/{wildcards.dataset}/align/genome/ ctrl-1 ctrl-2
#		elif [[ ${{thedataset}} == "pnas2016" ]]; then
#			~/Applications/zsh-5.9/bin/zsh tools/3.down_sample.sh output_q1_10_40_unique/{wildcards.dataset}/align/genome/ ctrl-rep1 ctrl-rep2
#		else
#			~/Applications/zsh-5.9/bin/zsh tools/3.down_sample.sh output_q1_10_40_unique/{wildcards.dataset}/align/genome/
#		fi
#
#		mkdir -p /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/output_q1_10_40_unique/{wildcards.dataset}/align/genome_down/
#		cd /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/output_q1_10_40_unique/{wildcards.dataset}/align/genome_down/
#		rm -rf *.bam *.bai
#		ln -sf ../genome/*_down_genome_sorted_rmdup.bam ./
#		~/miniforge3/bin/rename 's/_down_genome/_genome/g' *
#		samtools index -M -@ {threads} *_genome_sorted_rmdup.bam
#		touch down_finished.flag
#		"""
#
#rule divide_mip_downsample:
#	input:
#		downbam = "output_q1_10_40_unique/{dataset}/align/genome_down/down_finished.flag"
#	output:
#		downmip = "output_q1_10_40_unique/{dataset}/align/genome_down/downmip_finished.flag"
#	threads: 20
#	shell:
#		"""
#		cd output_q1_10_40_unique/{wildcards.dataset}/align/genome_down/
#		for sample in `ls -1 *_genome_sorted_rmdup.bam | cut -d '_' -f 1`
#		do
#			samtools view -b -F 276 -@ 20 -o ${{sample}}_genome_sorted_rmdup_fwd.bam ${{sample}}_genome_sorted_rmdup.bam && samtools index -@ 20 ${{sample}}_genome_sorted_rmdup_fwd.bam && samtools mpileup -f /gpfs1/chengqiyi_pkuhpc/limingyang/hg38/hg38_only_chromsomes.fa ${{sample}}_genome_sorted_rmdup_fwd.bam -o ${{sample}}_fwd.mip -B -d 10000000 -q 0 -Q 0 &
#			samtools view -b -F 260 -f 16 -@ 20 -o ${{sample}}_genome_sorted_rmdup_rvs.bam ${{sample}}_genome_sorted_rmdup.bam && samtools index -@ 20 ${{sample}}_genome_sorted_rmdup_rvs.bam && samtools mpileup -f /gpfs1/chengqiyi_pkuhpc/limingyang/hg38/hg38_only_chromsomes.fa ${{sample}}_genome_sorted_rmdup_rvs.bam -o ${{sample}}_rvs.mip -B -d 10000000 -q 0 -Q 0 &
#		done
#		wait
#		touch downmip_finished.flag
#		"""
