library(readr)
library(ggplot2)
#library(tidyverse)
library(dplyr)
library(stringr)
library(patchwork)
library(ggprism)
library(ggsci)
df <- read_tsv('2.flank_region/G_frac_all.txt')

df <- df %>% tidyr::separate(col=sample, c('con','type','strand'),'_', remove=FALSE)

mincut <- df$cut |> min()
maxcut <- df$cut |> max()


for(therate in c(10,50)){
	for(thecut in mincut:maxcut){
		dff <- df %>% dplyr::filter(cut == thecut, rate == therate)
		subfig <- ggplot(dff , aes(x=pos, y=fraction)) +
			geom_line(aes(linetype=con, color=type), linewidth=0.5, alpha=0.5) +
			geom_point(aes(shape=con, color=type), size=0.2, alpha=0.5) +
			theme_prism(base_size=6) +
			ylim(c(0,1)) +
			labs(subtitle = paste0('Cov>',thecut,' Rate>=',therate)) +
			geom_hline(yintercept=0.25, linetype='dashed')+
			facet_wrap(type ~strand) +
			scale_linetype_manual(values=c('treat'='solid','ctrl'='dashed'))
		if(thecut==mincut && therate==10){
			fig <- subfig
		} else {
			fig <- fig + subfig
		}
	}
}

fig <- fig + plot_layout(ncol=8, guides='collect')

#fig_fwd_oxa <- ggplot(df %>% filter(strand=='fwd', type=='oxa'), aes(x=pos, y=fraction)) +
#	geom_line(aes(linetype=con, color=anti)) +
#	geom_point(aes(shape=con, color=anti)) +
#	theme_prism() +
#	ylim(c(0,0.6)) +
#	geom_hline(yintercept=0.25, linetype='dashed')+
#	labs(subtitle='Oxa')
#
#fig_rvs_cis <- ggplot(df %>% filter(strand=='rvs', type=='cis'), aes(x=pos, y=fraction)) +
#	geom_line(aes(linetype=con, color=anti)) +
#	geom_point(aes(shape=con, color=anti)) +
#	theme_prism() +
#	ylim(c(0,0.6)) +
#	geom_hline(yintercept=0.25, linetype='dashed')+
#	labs(subtitle='Cis, Reverse Strand')
#
#fig_rvs_oxa <- ggplot(df %>% filter(strand=='rvs', type=='oxa'), aes(x=pos, y=fraction)) +
#	geom_line(aes(linetype=con, color=anti)) +
#	geom_point(aes(shape=con, color=anti)) +
#	theme_prism() +
#	ylim(c(0,0.6)) +
#	geom_hline(yintercept=0.25, linetype='dashed')+
#	labs(subtitle='Oxa')

#fig <- (fig_fwd_cis + fig_fwd_oxa) / (fig_rvs_cis + fig_rvs_oxa)
ggsave('2.flank_region/G_frequency.pdf', width=32, height=18)
