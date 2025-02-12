#!/bin/bash
workdir=$1
cut=$2
rate=$3
trts=$4

for trt in ${trts[@]}
do
	#treat
	awk -v con="treat" '{if($5=="GG" && $10==con) print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' ${workdir}/${trt}_pattern.bed > ${workdir}/${trt}_GG_treat.final_stopsite.bed
	#ctrl, actually all site not only GG
	awk -v con="ctrl" '{if($10==con) print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' ${workdir}/${trt}_pattern.bed > ${workdir}/${trt}_GG_ctrl_stopsite.bed
done

ls -1 ${workdir}/*GG*_stopsite.bed > ${workdir}/stopsite_bedlist_onlyGG.txt


python3 snk_script/5.1.cal_chr_distribution_onlyGG.py -d ${workdir} && \
Rscript snk_script/5.1.plot_bar_chrEnrich_onlyGG.r -i ${workdir}/stopsite_enrichment_onlyGG.txt -o ${workdir}/stopsite_enrichment_onlyGG.pdf
#python3 snk_script/5.1.cal_chr_distribution_onlyGG.py -d ${workdir} && \
