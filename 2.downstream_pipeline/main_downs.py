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
		outfa_trt = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat_expand_flank10_cisplatin_site.txt", cut = cuts, rate = rates, trt = trts, strand = strands),
		outfa_ctrl = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl_expand_flank10_cisplatin_site.txt", cut = cuts, rate = rates, trt = trts, strand = strands),
		sitenum = "2.flank_region/cpg_num.txt",
		outtrt_finalbed = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final.bed", cut = cuts, rate = rates, trt = trts, strand = strands),
		outctrlbed = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl.bed", cut = cuts, rate = rates, trt = trts, strand = strands),
		g_frac_flag = expand("2.flank_region/cut{cut}rate{rate}/g_frac_single_finished.flag", cut = cuts, rate = rates),
		gg_fraction = expand("2.flank_region/cut{cut}rate{rate}/gg_fraction.txt", cut = cuts, rate = rates), 
		treat_seqlogo = expand("2.flank_region/cut{cut}rate{rate}/treat_seqlogo.pdf", cut = cuts, rate = rates),
		ctrl_seqlogo = expand("2.flank_region/cut{cut}rate{rate}/ctrl_seqlogo.pdf", cut = cuts, rate = rates),
		pattern_fraction = expand("2.flank_region/cut{cut}rate{rate}/pattern_fraction.pdf", cut = cuts, rate = rates),
		gfrac_ATCG_flag = expand("2.flank_region/cut{cut}rate{rate}/g_frac_single_ATCG_finished.flag", cut = cuts, rate = rates),
		G_frac_all = "2.flank_region/G_frac_all.txt",
		stopsite_bedlist = expand("2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt", cut = cuts, rate = rates),
		stopsite_enrichment = expand("2.flank_region/cut{cut}rate{rate}/stopsite_enrichment.txt", cut = cuts, rate = rates),
		stopsite_enrichment_pdf = expand("2.flank_region/cut{cut}rate{rate}/stopsite_enrichment.pdf", cut = cuts, rate = rates),
		distribution = expand("2.flank_region/cut{cut}rate{rate}/distribution/annotation.pdf", cut = cuts, rate = rates),
		outflag = expand("2.flank_region/cut{cut}rate{rate}/pattern_finish.flag", cut = cuts, rate = rates),
		outposflag = expand("2.flank_region/cut{cut}rate{rate}/patternpos_finish.flag", cut = cuts, rate = rates),
		bamcov_rpkm = expand("bedGraph/{trt}_genome_sorted_rmdup_BPM_1mb.bedGraph", trt = trts),
		bamcov = expand("bedGraph/{trt}_genome_sorted_rmdup_1mb.bedGraph", trt = trts),
		box = expand("2.flank_region/cut{cut}rate{rate}/distance_box.pdf", cut = cuts, rate = rates),
		boxlog = expand("2.flank_region/cut{cut}rate{rate}/distance_box_log.pdf", cut = cuts, rate = rates),
		stopsite_bedlist_onlyGG = expand("2.flank_region/cut{cut}rate{rate}/stopsite_bedlist_onlyGG.txt", cut = cuts, rate = rates),
		stopsite_enrichment_onlyGG = expand("2.flank_region/cut{cut}rate{rate}/stopsite_enrichment_onlyGG.txt", cut = cuts, rate = rates),
		stopsite_enrichment_pdf_onlyGG = expand("2.flank_region/cut{cut}rate{rate}/stopsite_enrichment_onlyGG.pdf", cut = cuts, rate = rates),
		distribution_onlyGG = expand("2.flank_region/cut{cut}rate{rate}/distribution_onlyGG/annotation.pdf", cut = cuts, rate = rates),
		ptrn_signal_flag = expand("2.flank_region/cut{cut}rate{rate}/pattern_signal_finish.flag",cut = cuts, rate = rates),
		cluster_flag = expand("2.flank_region/cut{cut}rate{rate}/cluster_bed_finish.flag", cut = cuts, rate = rates),
		cluster_flaggg = expand("2.flank_region/cut{cut}rate{rate}/clustergg_bed_finish.flag", cut = cuts, rate = rates),
		clusterpdf_flag = expand("2.flank_region/cut{cut}rate{rate}/cluster_pdf_finished.flag", cut = cuts, rate = rates),
		clusterpdfgg_flag = expand("2.flank_region/cut{cut}rate{rate}/clustergg_pdf_finished.flag", cut = cuts, rate = rates)

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

rule get_seq_of_stopsite:
	input:
		outtrt_final = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final"
	output:
		outfa_trt = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat_expand_flank10_cisplatin_site.txt",
		outfa_ctrl = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl_expand_flank10_cisplatin_site.txt"
	threads: 2
	log: "logs/4.getseq/{trt}_{strand}_{cut}_{rate}.log"
	shell:
		r"""
		thestrand={wildcards.strand}
		if [[ $thestrand == 'fwd' ]]; then
			strandi='+'
		else
			strandi='-'
		fi

		awk -v strand="$strandi" '{{print $1"\t"$2-11"\t"$2+10"\t"$2"\t"$3"\t"strand}}' 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.treat.final > ./2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.treat_expand_flank10.bed
		cut -f 1-6 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.treat_expand_flank10.bed | awk '$2>0' | sort -u > 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.treat_expand_flank10_cisplatin_site.bed
		bedtools getfasta -fi {hg38} -bed 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.treat_expand_flank10_cisplatin_site.bed -s -fo 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.treat_expand_flank10_cisplatin_site.txt -tab

		awk -v strand="$strandi" '{{print $1"\t"$2-11"\t"$2+10"\t"$2"\t"$3"\t"strand}}' 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.ctrl > 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.ctrl_expand_flank10.bed
		cut -f 1-6 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.ctrl_expand_flank10.bed | awk '$2>0' | sort -u > 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.ctrl_expand_flank10_cisplatin_site.bed
		bedtools getfasta -fi {hg38} -bed 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.ctrl_expand_flank10_cisplatin_site.bed -s -fo 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/{wildcards.trt}_{wildcards.strand}_stop_articut_{wildcards.cut}_{wildcards.rate}.txt.ctrl_expand_flank10_cisplatin_site.txt -tab

		"""

rule stat_sitenum:
	input:
		outtrt = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat", cut = cuts, rate = rates, trt = trts, strand = strands),
		outctrl = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl", cut = cuts, rate = rates, trt = trts, strand = strands),
		outtrt_final = expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final", cut = cuts, rate = rates, trt = trts, strand = strands)
	output:
		sitenum = "2.flank_region/cpg_num.txt"
	threads: 5
	log: "logs/21.stat_sitenum/sitenum.log"
	shell:
		"""
		if [ -e cpg_num.txt ]; then
			rm -rf  cpg_num.txt
		fi

		for cut in {cutseq}
		do
			for rate in {rateseq}
			do
				for tp in {trtseq}
				do
					for strand in fwd rvs
					do
						for trt in treat.final treat ctrl
						do
							sample=${{tp}}_${{strand}}
							bash snk_script/1.stat_cpgnum.sh $sample $cut $trt $rate &
						done
					done
				done
			done
		done
		wait

		Rscript snk_script/1.stat_cpgnum.r
		"""

rule gen_bed_of_stopsites:
	input:
		outtrt_final = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final",
		outctrl = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl"
	output:
		outtrt_finalbed = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final.bed",
		outctrlbed = "2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl.bed"
	log:
		"logs/20.genbed/{trt}_{strand}_{cut}_{rate}.log"
	threads: 1
	shell:
		"""
		bash snk_script/0.gen_site_bed.sh {wildcards.cut} {wildcards.rate} {wildcards.trt} {wildcards.strand}
		"""


rule cal_g_frac:
	input:
		outfa_trt = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat_expand_flank10_cisplatin_site.txt", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands),
		outfa_ctrl = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl_expand_flank10_cisplatin_site.txt", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands)
	output:
		g_frac_flag = "2.flank_region/cut{cut}rate{rate}/g_frac_single_finished.flag"
	log:
		"logs/22.g_frac/{cut}_{rate}.log"
	threads: 10
	shell:
		"""
		bash snk_script/2.cal_g_frac_singlefile.sh 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}
		"""

rule cal_g_frac_ATCG:
	input:
		g_frac_flag = "2.flank_region/cut{cut}rate{rate}/g_frac_single_finished.flag"
	output:
		gg_fraction = "2.flank_region/cut{cut}rate{rate}/gg_fraction.txt",
		treat_seqlogo = "2.flank_region/cut{cut}rate{rate}/treat_seqlogo.pdf",
		ctrl_seqlogo = "2.flank_region/cut{cut}rate{rate}/ctrl_seqlogo.pdf",
		pattern_fraction = "2.flank_region/cut{cut}rate{rate}/pattern_fraction.pdf",
		gfrac_ATCG_flag = "2.flank_region/cut{cut}rate{rate}/g_frac_single_ATCG_finished.flag"
	log:
		"logs/22.g_frac_ATCG/{cut}_{rate}.log"
	threads: 10
	shell:
		"""
		bash snk_script/2.cal_g_frac_singlefile_atcg.sh 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}
		"""

rule merge_g_frac:
	input:
		g_frac_flag = expand("2.flank_region/cut{cut}rate{rate}/g_frac_single_finished.flag", cut = cuts, rate = rates),
		gg_fraction = expand("2.flank_region/cut{cut}rate{rate}/gg_fraction.txt", cut = cuts, rate = rates), 
		treat_seqlogo = expand("2.flank_region/cut{cut}rate{rate}/treat_seqlogo.pdf", cut = cuts, rate = rates),
		ctrl_seqlogo = expand("2.flank_region/cut{cut}rate{rate}/ctrl_seqlogo.pdf", cut = cuts, rate = rates),
		pattern_fraction = expand("2.flank_region/cut{cut}rate{rate}/pattern_fraction.pdf", cut = cuts, rate = rates),
		gfrac_ATCG_flag = expand("2.flank_region/cut{cut}rate{rate}/g_frac_single_ATCG_finished.flag", cut = cuts, rate = rates)
	output:
		G_frac_all = "2.flank_region/G_frac_all.txt"
	threads: 2
	shell:
		"""
		bash snk_script/3.merge_gfrac.sh "{cutseq}" "{rateseq}"
		"""


### stopsite analysis
rule chr_enrich:
	input:
		outfa_trt = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat.final", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands),
		outfa_ctrl = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands)
	output:
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt",
		stopsite_enrichment = "2.flank_region/cut{cut}rate{rate}/stopsite_enrichment.txt",
		stopsite_enrichment_pdf = "2.flank_region/cut{cut}rate{rate}/stopsite_enrichment.pdf",
	threads: 7
	shell:
		"""
		bash snk_script/5.chr_enrich.sh 2.flank_region/cut{wildcards.cut}rate{wildcards.rate} {wildcards.cut} {wildcards.rate} "{trtseq}"
		"""
	   
rule element_enrich:
	input:
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt",
	output:
		distribution = "2.flank_region/cut{cut}rate{rate}/distribution/annotation.pdf"
	threads: 5
	shell:
		"""
		bash snk_script/7genomic_element_distribution_annl_6mA.sh 2.flank_region/cut{wildcards.cut}rate{wildcards.rate} {trtseq}
		"""

rule divide_gg_not:
	input:
		outfa_trt = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.treat_expand_flank10_cisplatin_site.txt", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands),
		outfa_ctrl = lambda wildcards: expand("2.flank_region/cut{cut}rate{rate}/{trt}_{strand}_stop_articut_{cut}_{rate}.txt.ctrl_expand_flank10_cisplatin_site.txt", cut = wildcards.cut, rate = wildcards.rate, trt = trts, strand = strands)
	output:
		outflag = "2.flank_region/cut{cut}rate{rate}/pattern_finish.flag"
	threads: 15
	shell:
		"""
		cd 2.flank_region
		for i in {trtseq}
		do
			python3 ../snk_script/13.divide_gg_not.py -t $i -c {wildcards.cut} -r {wildcards.rate} &
		done

		touch cut{wildcards.cut}rate{wildcards.rate}/pattern_finish.flag
		wait
		"""

rule pattern_position:
	input:
		outflag = "2.flank_region/cut{cut}rate{rate}/pattern_finish.flag"
	output:
		outposflag = "2.flank_region/cut{cut}rate{rate}/patternpos_finish.flag"
	threads: 10
	shell:
		"""
		cd 2.flank_region
		for i in {trtseq}
		do
			Rscript ../snk_script/site_position.r -b cut{wildcards.cut}rate{wildcards.rate}/${{i}}_pattern.bed &
		done

		wait
		touch cut{wildcards.cut}rate{wildcards.rate}/patternpos_finish.flag
		"""

rule site_distance:
	input:
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt"
	output:
		box = "2.flank_region/cut{cut}rate{rate}/distance_box.pdf",
		boxlog = "2.flank_region/cut{cut}rate{rate}/distance_box_log.pdf" 
	threads: 4
	shell:
		"""
		Rscript snk_script/site_distance.r -b {input.stopsite_bedlist} -t "{trtseq} {ctrltrt}"
		"""



rule bam_cov:
	input:
		bam = "{trt}_genome_sorted_rmdup.bam"
	output:
		bamcov_rpkm = "bedGraph/{trt}_genome_sorted_rmdup_BPM_1mb.bedGraph",
		bamcov = "bedGraph/{trt}_genome_sorted_rmdup_1mb.bedGraph"
	threads: 20
	shell:
		"""
		samtools index -@ 20 {input.bam}
		bash snk_script/0.1.bam_to_bedGraph.sh {input.bam}
		"""


##### only GG analysis

rule chr_enrich_onlyGG:
	input:
		outflag = "2.flank_region/cut{cut}rate{rate}/pattern_finish.flag"
	output:
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist_onlyGG.txt",
		stopsite_enrichment = "2.flank_region/cut{cut}rate{rate}/stopsite_enrichment_onlyGG.txt",
		stopsite_enrichment_pdf = "2.flank_region/cut{cut}rate{rate}/stopsite_enrichment_onlyGG.pdf",
	threads: 7
	shell:
		"""
		echo "Here generate GG bed"
		bash snk_script/5.1.chr_enrich_onlyGG.sh 2.flank_region/cut{wildcards.cut}rate{wildcards.rate} {wildcards.cut} {wildcards.rate} "{trtseq}"
		"""
	   
rule element_enrich_onlyGG:
	input:
		outflag = "2.flank_region/cut{cut}rate{rate}/pattern_finish.flag",
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist_onlyGG.txt"
	output:
		distribution = "2.flank_region/cut{cut}rate{rate}/distribution_onlyGG/annotation.pdf"
	threads: 5
	shell:
		"""
		bash snk_script/7.1.genomic_element_distribution_annl_6mA_onlyGG.sh 2.flank_region/cut{wildcards.cut}rate{wildcards.rate} {trtseq}
		"""

rule pattern_signal:
	input:
		outflag = "2.flank_region/cut{cut}rate{rate}/pattern_finish.flag"
	output:
		ptrn_signal_flag = "2.flank_region/cut{cut}rate{rate}/pattern_signal_finish.flag"
	threads: 10
	shell:
		"""
		for i in {trtseq}
		do
			python3 snk_script/16.pattern_signal.py -t $i -c {wildcards.cut} -r {wildcards.rate} &
		done
		touch {output.ptrn_signal_flag}
		wait
		"""

rule cluster:
	input:
		stopsite_bedlist = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist.txt",
		stopsite_bedlistgg = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist_onlyGG.txt",
	output:
		cluster_flag = "2.flank_region/cut{cut}rate{rate}/cluster_bed_finish.flag",
		cluster_flaggg = "2.flank_region/cut{cut}rate{rate}/clustergg_bed_finish.flag",
		clusterpdf_flag = "2.flank_region/cut{cut}rate{rate}/cluster_pdf_finished.flag",
		clusterpdfgg_flag = "2.flank_region/cut{cut}rate{rate}/clustergg_pdf_finished.flag",
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

		# process for non GG
		for pl in $(ls -1 2.flank_region/cut3rate10/*_cluster_pos.bed | grep -v GG | grep -v ctrl)
		do
			Rscript snk_script/cluster_position.r -b ${{pl}} -c 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/ctrl_treat.final_stopsite_cluster_pos.bed &
		done
		wait
		touch 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/cluster_pdf_finished.flag



		ctrl=0
		for pl in $(ls -1 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/*_stopsite.bed | grep GG)
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
		cd 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/ && ln -sf ${{ctrltrt}}_ctrl_stopsite_cluster_pos.bed ctrl_GG_treat.final_stopsite_cluster_pos.bed && touch clustergg_bed_finish.flag
		cd ../../

		# process for GG
		for pl in $(ls -1 2.flank_region/cut3rate10/*_cluster_pos.bed | grep GG | grep -v ctrl)
		do
			Rscript snk_script/cluster_position.r -b ${{pl}} -c 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/ctrl_GG_treat.final_stopsite_cluster_pos.bed &
		done
		wait
		touch 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/clustergg_pdf_finished.flag
		"""

#rule clusterpos:
#	input:
#		cluster_flag = "2.flank_region/cut{cut}rate{rate}/cluster_bed_finish.flag"
#	output:
#		cluster_pospdf = "2.flank_region/cut{cut}rate{rate}/{trt}_treat.final_stopsite_cluster_pos.bed_distribution.pdf"
#	params:
#		cluster_pos = "2.flank_region/cut{cut}rate{rate}/{trt}_treat.final_stopsite_cluster_pos.bed"
#	threads: 2
#	shell:
#		r"""
#		Rscript snk_script/cluster_position.r -b {params.cluster_pos} -c 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/ctrl_treat.final_stopsite_cluster_pos.bed
#		"""
#rule cluster_GG:
#	input:
#		stopsite_bedlistgg = "2.flank_region/cut{cut}rate{rate}/stopsite_bedlist_onlyGG.txt"
#	output:
#		cluster_flaggg = "2.flank_region/cut{cut}rate{rate}/clusterGG/clusterGG_bed_finish.flag"
#	threads: 2
#	shell:
#		r"""
#		ctrl=0
#		for pl in $(cat {input.stopsite_bedlistgg})
#		do
#			if [[ $(echo ${{pl}} | grep -oP "(ctrl|treat.final)") == "ctrl" ]] && [[ ${{ctrl}} == 0 ]]; then
#				ctrl=1
#				ctrltrt=$(echo ${{pl}} | grep -oP "(?<={wildcards.rate}\/).*(?=_ctrl)")
#			elif [[ $(echo ${{pl}} | grep -oP "(ctrl|treat.final)") == "ctrl" ]] && [[ ${{ctrl}} == 1  ]]; then
#				continue
#			fi
#			plbase=echo "${pl}" | xargs -n 1 basenaem
#			bash snk_script/12.cluster.sh ${{pl}} 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/clusterGG/${{plbase%%.bed}}_cluster.bed
#		done
#		cd 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/clusterGG && ln -sf ${{ctrltrt}}_ctrl_stopsite_cluster_pos.bed ctrl_treat.final_stopsite_cluster_pos.bed && touch clusterGG_bed_finish.flag
#		"""
#
#rule clusterpos_GG:
#	input:
#		cluster_flaggg = "2.flank_region/cut{cut}rate{rate}/clusterGG/clusterGG_bed_finish.flag"
#	output:
#		cluster_pospdfgg = "2.flank_region/cut{cut}rate{rate}/clusterGG/{trt}_GG_treat.final_stopsite_cluster_pos.bed_distribution.pdf"
#	params:
#		cluster_pos = "2.flank_region/cut{cut}rate{rate}/clusterGG/{trt}_GG_treat.final_stopsite_cluster_pos.bed"
#	threads: 2
#	shell:
#		r"""
#		Rscript snk_script/cluster_position.r -b {params.cluster_pos} -c 2.flank_region/cut{wildcards.cut}rate{wildcards.rate}/clusterGG/ctrl_treat.final_stopsite_cluster_pos.bed
#		"""
