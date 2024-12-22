: <<'END'
Copyright (c) 2024-03-25 by LiMingyang, YiLab, Peking University.

Author: Li Mingyang (limingyang200101@gmail.com)

Institute: AAIS, Peking University

File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/tools/downstream_pipeline/snk_script/15.gg_fraction_by_cutoff.sh

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

IFS=' ' read -r -a cuts <<< "$1"
IFS=' ' read -r -a rates <<< "$2"

echo -e "treatment\tpattern\tcount\ttotal_count\tfraction\tcut\trate" > 2.flank_region/gg_fraction_all.txt

for cut in "${cuts[@]}"
do
	for rate in "${rates[@]}"
	do
		sed '1d' 2.flank_region/cut${cut}rate${rate}/gg_fraction.txt | awk -v cut=$cut -v rate=$rate '{print $0"\t"cut"\t"rate}' >> 2.flank_region/gg_fraction_all.txt
	done
done
