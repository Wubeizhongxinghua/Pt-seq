#!/bin/bash
sample=$1
cutoff=$2
trt=$3
rate=$4
num=$(wc -l 2.flank_region/cut${cutoff}rate${rate}/${sample}_stop_articut_${cutoff}_${rate}.txt.${trt} | cut -d ' ' -f 1)
echo -e "${sample}\t${trt}\t${cutoff}\t${rate}\t${num}" >> 2.flank_region/cpg_num.txt
