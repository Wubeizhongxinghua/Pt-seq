#!/bin/bash

cut=$1
rate=$2
t=$3
strand=$4
workdir="2.flank_region/cut${cut}rate${rate}"

rmdup () {
    awk '!a[$0]++' $@
}

#for t in "car-1" "car-2" "cis-1" "cis-5" "lo-1" "lo-2" "pt-1" "pt-2"
#do
#    for strand in fwd rvs
#    do
        awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ${workdir}/${t}_${strand}_stop_articut_${cut}_${rate}.txt.treat.final > ${workdir}/${t}_${strand}_stop_articut_${cut}_${rate}.txt.treat.final.bed
        awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ${workdir}/${t}_${strand}_stop_articut_${cut}_${rate}.txt.ctrl > ${workdir}/${t}_${strand}_stop_articut_${cut}_${rate}.txt.ctrl.bed
#    done
#done
