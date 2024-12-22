library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggprism)
library(argparse)
parser <- ArgumentParser() #创建参数解析对象

parser$add_argument("-i","--input",help='Input for plot file') #添加参数
parser$add_argument("-o","--output",help='Output fig file.') #添加参数
args <- parser$parse_args()

df <- read_tsv(args$input)

elements <- c('3UTR', 'Retroposon', 'RC?', 'RNA', 'miRNA', 'ncRNA', 'TTS', 'LINE', 'srpRNA', 'SINE', 'RC', 'tRNA', 'DNA?', 'pseudo', 'DNA', 'Exon', 'Intron', 'Intergenic', 'Promoter', '5UTR', 'snoRNA', 'LTR?', 'scRNA', 'CpG-Island', 'Low_complexity', 'LTR', 'Simple_repeat', 'snRNA', 'Unknown', 'SINE?', 'Satellite', 'rRNA')

df <- df %>% gather('3UTR', 'Retroposon', 'RC?', 'RNA', 'miRNA', 'ncRNA', 'TTS', 'LINE', 'srpRNA', 'SINE', 'RC', 'tRNA', 'DNA?', 'pseudo', 'DNA', 'Exon', 'Intron', 'Intergenic', 'Promoter', '5UTR', 'snoRNA', 'LTR?', 'scRNA', 'CpG-Island', 'Low_complexity', 'LTR', 'Simple_repeat', 'snRNA', 'Unknown', 'SINE?', 'Satellite', 'rRNA', key='Element', value='log2(obs/exp)')

plot <- ggplot(data=df,aes(x=factor(Element,levels=elements), y=`log2(obs/exp)`, fill=Sample))+
	geom_col(position='dodge',color='black') + 
	labs(
		x='Elements',
		y='Enrichment'
	) +
	facet_grid(rows=vars(Con)) + theme_prism(base_size=14)

ggsave(args$output)

