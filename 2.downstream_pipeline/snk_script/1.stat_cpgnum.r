library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(ggsci)
library(tidyverse)
library(patchwork)
library(ggprism)
library(rayshader)
df <- read_tsv('2.flank_region/cpg_num.txt', col_names=c('sample','treatment','cutoff','rate','number'))
df <- df %>% separate(col=sample, c('type','strand'), '_', remove=FALSE)
#df$number <- log(df$number)

types <- df$type |> unique()
i <- 0

#df <- df |> mutate(
#                   antibody = case_when(
#                                        treatment == 'ctrl' ~ "5ug",
#                                        TRUE ~ antibody
#                                        )
#                    
#                   )

df <- df |> distinct()


for(atype in types){
	
	fig_cis <- ggplot(data=df %>% filter(type==atype, cutoff>=1), aes(x=cutoff,y=rate, fill=log(number))) +
		geom_tile(color='white',lwd=0.5, linetype=0.5) + 
		geom_text(aes(label = number), color = "white", size = 2) +
		labs(
			subtitle=atype,
			x = 'Cutoff',
			y = 'Rate',
			fill='Number' ) + 
		scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8)) + 
		scale_y_continuous(breaks=c(10,30,50,70,90)) +
		facet_grid(treatment ~ strand, scales='fixed') +
		theme(legend.text = element_text(size = 4)) +
		theme_prism(base_size=6)

    

	if( i == 0){
		fig <- fig_cis
		i <- i + 1
	} else {
		fig <- fig / (fig_cis & labs(x=NULL, y=NULL))
	}

}

figp <- ggplot(data=df |> filter(cutoff>=1, rate==10), aes(x=cutoff, y=number, color=type)) +
    geom_line(alpha=0.6) +
    geom_point(aes(fill=type), size=2, shape=21, color='black') +
    facet_grid(treatment~strand) +
    theme_prism() +
    coord_cartesian(ylim=c(0,100000)) +
    scale_color_aaas() +
    scale_fill_aaas()

ggsave('2.flank_region/stat_stopsite.pdf', fig, height=16, width=16)
ggsave('2.flank_region/stat_stopsite_point.pdf', figp, height=9, width=16)
