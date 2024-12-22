#!/bin/bash

cutseq=$1
rateseq=$2
basedir="2.flank_region"
if [ -f ${basedir}/G_frac_all.txt ]; then
	rm -rf ${basedir}/G_frac_all.txt 
fi

for cut in ${cutseq[@]}
do
	for rate in ${rateseq[@]}
	do
		workdir="2.flank_region/cut${cut}rate${rate}"
		cat ${workdir}/G_frac*.txt.txt | sed '/sample/d' > ${workdir}/G_frac_all.txt
		sed "s/$/\t${cut}\t${rate}/g" ${workdir}/G_frac_all.txt >> ${basedir}/G_frac_all.txt

		#cat ${workdir}/G_frac*.txt_atcg.txt | sed '/sample/d' > ${workdir}/G_frac_all_atcg.txt
        #sed "s/$/\t${cut}\t${rate}/g" ${workdir}/G_frac_all_atcg.txt >> ${basedir}/G_frac_all_atcg.txt
	done
done


sed -i "1isample\tpos\tfraction\tcut\trate" ${basedir}/G_frac_all.txt
#sed -i "1isample\tpos\tfraction\tbase\tcut\trate" ${basedir}/G_frac_all_atcg.txt
Rscript snk_script/frequency_plot.r
