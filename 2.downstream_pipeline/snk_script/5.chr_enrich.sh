#!/bin/bash
workdir=$1
cut=$2
rate=$3
trts=$4

for trt in ${trts[@]}
do                                                                                                             
		for con in treat.final ctrl
		do
			for strand in fwd rvs 
			do
				awk '{print $1"\t"$2"\t"$2+1}' ${workdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.${con} > ${workdir}/${trt}_${strand}_${con}_stopsite.bed
			done
			cat ${workdir}/${trt}_fwd_${con}_stopsite.bed ${workdir}/${trt}_rvs_${con}_stopsite.bed | sort -k1,1 -k2,2n > ${workdir}/${trt}_${con}_stopsite.bed 
			rm -rf ${workdir}/${trt}_fwd_${con}_stopsite.bed ${workdir}/${trt}_rvs_${con}_stopsite.bed
		done
done  

ls -1 ${workdir}/*_stopsite.bed | grep -v GG > ${workdir}/stopsite_bedlist.txt


python3 snk_script/5.cal_chr_distribution.py -d ${workdir} && \
Rscript snk_script/5.plot_bar_chrEnrich.r -i ${workdir}/stopsite_enrichment.txt -o ${workdir}/stopsite_enrichment.pdf
