: <<'END'
Copyright (c) 2024-06-17 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/figs/fig3/ana_sites/12.cluster.sh

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
END

source ~/env_backup/.functions
#!/bin/bash

bed=$1 #input
clusterbed=$2 #output

sort -k1,1 -k2,2n ${bed} | bedtools merge -i - -d 10000 |\
    bedtools intersect -a - -b ${bed} -wa |\
    uniq -c - |\
    awk '$1>5 {$1=""; print substr($0, 2)}' - |\
    sed "s/ /\t/g" |\
    bedtools intersect -a - -b ${bed} -wa -wb > ${clusterbed}


sort -k1,1 -k2,2n ${bed} | bedtools merge -i - -d 10000 |\
    bedtools intersect -a - -b ${bed} -wa |\
    uniq -c - |\
    awk '$1>5 {$1=""; print substr($0, 2)}' - |\
    sed "s/ /\t/g" |\
    bedtools intersect -a - -b ${bed} -c > ${clusterbed%%.bed}_pos.bed
