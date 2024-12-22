#!/bin/bash
trt=$1
strand=$2
rate=$3
cut=$4
outdir=$5

awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat | sort --parallel=10 -k1,1 -k2,2n > ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat.bed
awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.ctrl | sort --parallel=10 -k1,1 -k2,2n > ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.ctrl.bed
bedtools intersect -sorted -v -wa -a ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat.bed -b ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.ctrl.bed > ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat.final.bed
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat.final.bed > ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat.final
rm -rf ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat.bed ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.ctrl.bed ${outdir}/${trt}_${strand}_stop_articut_${cut}_${rate}.txt.treat.final.bed
