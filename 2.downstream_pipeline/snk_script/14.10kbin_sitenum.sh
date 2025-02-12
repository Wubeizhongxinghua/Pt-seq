: <<'END'
Copyright (c) 2024-03-08 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/tools/downstream_pipeline/snk_script/14.10kbin_sitenum.sh

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
END

#source ~/env_backup/.functions
#!/usr/bin/bash


cut=$1
rate=$2
trtsinput=${@:3}
IFS=' ' read -r -a trts <<< "$trtsinput"


workdir="2.flank_region/cut${cut}rate${rate}"



mkdir -p ${workdir}/bin10k

for i in ${trts[@]}
do
	awk '{if($10=="treat") print $1"\t"$3"\t"$3+1"\t"$6"\t1\t"$5}' ${workdir}/${i}_pattern.bed | sort -k1,1 -k2,2n - | bedtools map -a snk_script/hg38_10K_bin_sort.bed -b - -o count > ${workdir}/bin10k/${i}_10k.bed &
done

awk '{if($10=="ctrl") print $1"\t"$3"\t"$3+1"\t"$6"\t1\t"$5}' ${workdir}/${trts[1]}_pattern.bed | sort -k1,1 -k2,2n - | bedtools map -a snk_script/hg38_10K_bin_sort.bed -b - -o count > ${workdir}/bin10k/ctrl_10k.bed &

wait
