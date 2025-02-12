'''
Copyright (c) 2024-03-02 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: main_downs.py

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
'''

"""
Generate stop sites (ctrl, treat, treat.final) and seq of stop sites.
"""

# You need to set:
trts = ['sample1']
ctrltrt = 'sample2'



strands = ['fwd', 'rvs']
trttypes = ["treat.final", "treat", "ctrl"]

hg38 = 'hg38/hg38_only_chromsomes.fa'

rates = [10, 50]
cuts = [1, 2, 3, 4, 5, 6, 7]
rateseq = ' '.join([str(ele) for ele in rates])
cutseq = ' '.join([str(ele) for ele in cuts])
trtseq = ' '.join(trts)
strandseq = ' '.join(strands)



rule all:
	input:
		outsig = expand("{trt}_sig_proba.txt", trt = trts),
		treat_stop = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat", cut = cuts, rate = rates, trt = trts, strand = strands),
		ctrl_stop = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl", cut = cuts, rate = rates, trt = trts, strand = strands),
		outtrt_final = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final", cut = cuts, rate = rates, trt = trts, strand = strands),
		stopsite_bedlist = expand("2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt", cut = cuts, rate = rates),
		cluster_flag = expand("2.flank_region/cut{cut}rate{rate}/cluster_bed_finish.flag", cut = cuts, rate = rates)

rule get_sig_bg:
	input:	
		mipfwd = "{trt}_fwd_sig.txt",
		miprvs = "{trt}_rvs_sig.txt"
	output:
		outsig = "{trt}_sig_proba.txt"
	threads: 1
	shell:
		"""
		python3 snk_script/0.0.cal_sig_bg.py -if {input.mipfwd} -ir {input.miprvs} -o {output.outsig}
		"""

rule get_stop:
	input:
		outsig = "{trt}_sig_proba.txt",
		miptrt = "{trt}_{strand}_sig.txt",
		mipctrl = ctrltrt+"_{strand}_sig.txt" #modify?
	output:
		outtrt = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat",
		outctrl = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl"
	threads: 1
	log: "logs/2.get_stop/{trt}_{strand}_{cut}_{rate}.log"
	shell:
		"""
		sigproba=$(cat {input.outsig})
		python3 snk_script/cisplatin_get_stop.py -ctrl {input.mipctrl} -treat {input.miptrt} -O 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt -stop_reads_cut {wildcards.cut} -stop_rate_cut {wildcards.rate} -cov_cut {wildcards.cut} -p_cut 0.05 -proba ${{sigproba}}
		"""

rule treat_rm_ctrl:
	input:
		outtrt = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat",
		outctrl = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl"
	output:
		outtrt_final = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final"
	threads: 1
	log: "logs/3.rm_ctrl/{trt}_{strand}_{cut}_{rate}.log"
	shell:
		"""
		bash snk_script/3.treat_rm_ctrl.sh {wildcards.trt} {wildcards.strand} {wildcards.rate} {wildcards.cut} 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}
		"""


rule chr_enrich:
	input:
		outfa_trt = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands),
		outfa_ctrl = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands)
	output:
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt",
	threads: 2
	shell:
		"""
		bash snk_script/5.get_bedlist.sh 2.flank_region/cut{wildcards.cut}rate{wildcards.rate} {wildcards.cut} {wildcards.rate} "{trtseq}"
		"""

rule cluster:
	input:
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt"
	output:
		cluster_flag = "2.flank_region/cut{cut}rate{rate}/cluster_bed_finish.flag"
	threads: 2
	shell:
		r"""
		ctrl=0
		for pl in $(ls -1 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/*_stopsite.bed | grep -v GG)
		do
			if [[ $(echo ${{pl}} | grep -oP "(ctrl|treat.final)") == "ctrl" ]] && [[ ${{ctrl}} == 0 ]]; then
				ctrl=1
				ctrltrt=$(echo ${{pl}} | grep -oP "(?<={wildcards.rate}\/).*(?=_ctrl)")
			elif [[ $(echo ${{pl}} | grep -oP "(ctrl|treat.final)") == "ctrl" ]] && [[ ${{ctrl}} == 1  ]]; then
				continue
			fi

			bash snk_script/12.cluster.sh ${{pl}} ${{pl%%.bed}}_cluster.bed
			#awk '{{print $1"\t"$2"\t"$3}}' ${{pl%%.bed}}_cluster.bed | awk '!a[$0]++' > ${{pl%%.bed}}_clusterpos.bed 
		done
		cd 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/ && ln -sf ${{ctrltrt}}_ctrl_stopsite_cluster_pos.bed ctrl_treat.final_stopsite_cluster_pos.bed && touch cluster_bed_finish.flag
		cd ../../
		"""

