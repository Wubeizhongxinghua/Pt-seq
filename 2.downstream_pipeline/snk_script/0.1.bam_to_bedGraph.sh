: <<'END'
Copyright (c) 2024-03-01 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/output/batch4/align/genome/0.1.bam_to_bedGraph.sh

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
#!/usr/bin/bash

bam=$1


mkdir -p bedGraph

/home/chengqiyi_pkuhpc/profiles/limingyang/miniforge3/envs/py310/bin/bamCoverage -b ${bam} -of bedgraph -o bedGraph/${bam%%.bam}_BPM_1mb.bedGraph -bs 1000000 --normalizeUsing BPM -p 10 &
/home/chengqiyi_pkuhpc/profiles/limingyang/miniforge3/envs/py310/bin/bamCoverage -b ${bam} -of bedgraph -o bedGraph/${bam%%.bam}_1mb.bedGraph -bs 1000000 -p 10 &
wait
