#!/bin/bash
workdir=$1

ls ${workdir}/*_site.txt | xargs -n 1 basename > ${workdir}/list.txt

Rscript snk_script/2.cal_g_frac_singlefile_atcg.r -d ${workdir}
#pkurun-cnlong 1 4 Rscript 2.cal_g_frac_singlefile_atcg_bit.r -d ${workdir}
#sleep 1
touch ${workdir}/g_frac_single_ATCG_finished.flag
