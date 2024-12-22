library(GenomicRanges)
library(rtracklayer)
library(ggridges)
library(ggplot2)
library(ggprism)
library(ggsci)
library(argparse)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

parser <- ArgumentParser() 
parser$add_argument('-b','--bed',help='Bedlist file input')
parser$add_argument('-t','--treat',help='treat numbers (last is ctrl), space separated')

args <- parser$parse_args()

treats <- strsplit(args$treat, " ")[[1]]
print('Input treats are:')
print(treats)

thedir <- str_extract(args$bed, "^.*(?=/)")

bedlist <- read.table(args$bed, header=FALSE, sep='\t', col.names=c('list'))

ctk <- 0


print('Read beds....')
beds <- lapply(bedlist$list, function(bedfile){
    con <- str_extract(bedfile, "(ctrl|treat.final)")
    if((con == 'ctrl') && (ctk == 0)){
        ctk <<- 1
        bed_df <<- read.table(bedfile, sep = "\t", header = FALSE,col.names = c("chromosome", "start", "end")) 
        bed_df$type <- treats[length(treats)] 
        return(bed_df)
    } else if((con=='ctrl') && (ctk != 0)) {
        dff <- data.frame(matrix(ncol=4, nrow=0))
        colnames(dff) <- c("chromosome", "start", "end", "type")
        dff$chromosome <- as.character(dff$chromosome)
        dff$start <- as.numeric(dff$start)
        dff$end <- as.numeric(dff$end)
        dff$type <- as.character(dff$type)
        return(dff) 
    } else {
        thetype <- str_extract(bedfile, "(?<=[0-9]/)[^_]+")
        tryCatch(
			expr={
        		bed_df <<- read.table(bedfile, sep = "\t", header = FALSE,col.names = c("chromosome", "start", "end")) 
        		bed_df$type <- thetype
			},
			error = function(e){
				bed_df <- data.frame(chromosome=c(), start=c(), end=c(), type=c())
			}
        )
        return(bed_df)
    }
})
beds <- keep(beds, ~ nrow(.x) > 0) #filter null dfs
bed_df <- bind_rows(beds)
#print("Read Random...")
#dfran <- read.table('random_8000_GG.bed', header=FALSE, col.names = c('chromosome','start','end','strand','GG'))
#dfran <- dfran |> select(chromosome, start, end)
#dfran$type <- 'Random'
#
#bed_df <- bed_df |> bind_rows(dfran)


bed <- GRanges(
    seqnames = bed_df$chromosome,
    ranges = IRanges(start = bed_df$start, end = bed_df$end),
    type = bed_df$type
)

bed_by_type <- split(bed, bed$type)
print('Calculating distance...')
distances <- lapply(names(bed_by_type), function(x) {
    data <- bed_by_type[[x]]
    dis <- distanceToNearest(data)
    disv <- mcols(dis)$distance
    tryCatch(
    expr={
		return(
			data.frame(distance=disv, platin=x)
		)
    },
	error = function(e){
		return(
			   data.frame(distance=c(), platin=c())
		)
	}
    )
})

disdf <- bind_rows(distances)
print('Plotting...')
fig <- ggplot(disdf, aes(x=factor(platin, levels=treats), y=distance, fill=platin)) +
    geom_boxplot(color='black', show.legend=FALSE, outlier.alpha=0.1, outlier.size=0.5) +
    labs(x='', y='Distance') +
    coord_cartesian(ylim=c(0,1000000)) +
    scale_fill_aaas()+
    theme_prism(axis_text_angle=45)

ggsave(paste0(thedir,"/distance_box.pdf"))


fig <- ggplot(disdf, aes(x=factor(platin, levels=treats), y=distance+1, fill=platin)) +
    geom_boxplot(color='black', show.legend=FALSE, outlier.alpha=0.1, outlier.size=0.5) +
    labs(x='', y='log10(Distance+1)') +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_fill_aaas()+
    theme_prism(axis_text_angle=45) +
    annotation_logticks(sides = "l")  

ggsave(paste0(thedir,"/distance_box_log.pdf"))
