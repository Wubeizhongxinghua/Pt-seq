library(dplyr)
library(forcats)
library(readr)
library(ggplot2)
library(tidyr)
library(argparse)
library(ggprism)
library(ggsci)
parser <- ArgumentParser() #创建参数解析对象

parser$add_argument("-i","--input",help='Input for plot file') #添加参数
parser$add_argument("-o","--output",help='Output fig file.') #添加参数
args <- parser$parse_args()

df <- read_tsv(args$input)

#df <- df %>% gather('HyperEnrichment','HypoEnrichment',key='Condition',value='log2(obs/exp)') %>% filter(chr != 'chrY')

plot <- ggplot(data=df,aes(x=factor(chr, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM')), y=log2(enrichment), fill=chr))+
	geom_col(color='black', show.legend=FALSE) + 
	labs(
		x='Chromosomes',
		y='Enrichment'
	) + 
	facet_wrap(treatment ~  condition, scales='free') + 
	theme_prism(axis_text_angle=45, base_size=7)

ggsave(args$output,width=16,height=9)

