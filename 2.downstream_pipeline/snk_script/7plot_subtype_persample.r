library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggprism)
library(stringr)
library(argparse)
parser <- ArgumentParser() #创建参数解析对象

parser$add_argument("-i","--input",help='Input dir') #添加参数
parser$add_argument("-o","--output",help='Output fig file.') #添加参数
args <- parser$parse_args()

files <- list.files(path=args$input, pattern="*_broadelement.txt")

dflist <- lapply(files, function(dfname){
					 df <- read_tsv(paste0(args$input,'/',dfname)) |>
					 	 select(Annotation, `Log2 Ratio (obs/exp)`) |>
					 	 rename(Element = Annotation)
					 df$Sample <- str_extract(dfname, ".*(?=_broadelement.txt)")
					 if(nrow(df) == 0){
					 } else {
					 	return(df)
					 }
})

df <- bind_rows(dflist)

elements <- c('3UTR', 'Retroposon', 'RC?', 'RNA', 'miRNA', 'ncRNA', 'TTS', 'LINE', 'srpRNA', 'SINE', 'RC', 'tRNA', 'DNA?', 'pseudo', 'DNA', 'Exon', 'Intron', 'Intergenic', 'Promoter', '5UTR', 'snoRNA', 'LTR?', 'scRNA', 'CpG-Island', 'Low_complexity', 'LTR', 'Simple_repeat', 'snRNA', 'Unknown', 'SINE?', 'Satellite', 'rRNA')

#df <- df %>% gather('3UTR', 'Retroposon', 'RC?', 'RNA', 'miRNA', 'ncRNA', 'TTS', 'LINE', 'srpRNA', 'SINE', 'RC', 'tRNA', 'DNA?', 'pseudo', 'DNA', 'Exon', 'Intron', 'Intergenic', 'Promoter', '5UTR', 'snoRNA', 'LTR?', 'scRNA', 'CpG-Island', 'Low_complexity', 'LTR', 'Simple_repeat', 'snRNA', 'Unknown', 'SINE?', 'Satellite', 'rRNA', key='Element', value='log2(obs/exp)')


plot <- ggplot(data=df,aes(x=factor(Element,levels=elements), y=`Log2 Ratio (obs/exp)`))+
	geom_col(aes(fill = ifelse(`Log2 Ratio (obs/exp)` > 1, 'red',
							  ifelse(`Log2 Ratio (obs/exp)` < -1, 'blue', 'grey'))),
			 position='dodge',color='black', show.legend = FALSE) +
	labs(
		x='Elements',
		y='Enrichment'
	) +
	facet_grid(rows=vars(Sample)) + theme_prism(base_size=8, axis_text_angle = 45) +
	scale_fill_manual(values = c('#008E9B', 'grey','#D65DB1'))

ggsave(args$output)

