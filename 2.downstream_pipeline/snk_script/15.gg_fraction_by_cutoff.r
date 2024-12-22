#Copyright (c) 2024-03-25 by LiMingyang, YiLab, Peking University.
#
#Author: Li Mingyang (limingyang200101@gmail.com)
#
#Institute: AAIS, Peking University
#
#File Name: /gpfs1/chengqiyi_pkuhpc/limingyang/cisplatin/tools/downstream_pipeline/snk_script/15.gg_fraction_by_cutoff.r
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

default_library <- 1

if(default_library){
	libs <- c('ggplot2','dplyr','readr','tidyr','ggsci','patchwork','ggprism','argparse')
	lapply(
		libs,
		function(lib){
			suppressPackageStartupMessages(library(lib, character.only=TRUE))
		}
	)
	if ('package:argparse' %in% search()){
		parser <- ArgumentParser()
		parser$add_argument("-n","--number",help='explanation')
		args <- parser$parse_args()
	}
}

df <- read_tsv('2.flank_region/gg_fraction_all.txt')


fig1 <- ggplot(df, aes(x=rate, y=cut, fill=fraction)) +
	geom_raster()+
	facet_wrap(treatment~., ncol=4)
ggsave('2.flank_region/gg_fraction_by_cutoff.pdf', fig1)


dfall <- df |> tidyr::separate(col=treatment, c('platin','strand'),'_') |>
	group_by(platin, cut, rate) |>
	summarize(GGcount = sum(count), totalcount = sum(total_count)) |>
	mutate(fraction = GGcount/totalcount)
fig2 <- ggplot(dfall, aes(x=rate, y=cut, fill=fraction)) +
	geom_raster()+
	facet_wrap(platin~., ncol=5) + 
	geom_text(aes(label=round(fraction,2)))
ggsave('2.flank_region/ggall_fraction_by_cutoff.pdf', fig2, width=16,height=9)


