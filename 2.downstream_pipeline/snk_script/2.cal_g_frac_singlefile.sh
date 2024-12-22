#!/bin/bash
workdir=$1

ls ${workdir}/*_site.txt | xargs -n 1 basename > ${workdir}/list.txt
for thefile in `cat ${workdir}/list.txt`
do
	python3 snk_script/2.cal_g_frac_single.py -i ${thefile} -d ${workdir} &
	sleep 0.5
done
wait
#Rscript frequency_plot.r 
touch ${workdir}/g_frac_single_finished.flag
