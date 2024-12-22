library(karyoploteR)
library(readr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('-b','--bed', help='Cluster bed files.')
parser$add_argument('-c','--ctbed', help='Ctrl cluster bed files.')
args <- parser$parse_args()
bedfile <- args$bed
ctbedfile <- args$ctbed

bed_df <- read.table(bedfile, sep = "\t", header = FALSE, 
    col.names = c("chromosome", "start", "end", "num"))

bed_ctdf <- read.table(ctbedfile, sep = "\t", header = FALSE, 
    col.names = c("chromosome", "start", "end", "num"))

bed <- GRanges(
    seqnames = bed_df$chromosome,
    ranges = IRanges(start = bed_df$start, end = bed_df$end)
)
ctbed <- GRanges(
    seqnames = bed_ctdf$chromosome,
    ranges = IRanges(start = bed_ctdf$start, end = bed_ctdf$end)
)

pdf(paste0(bedfile,'_distribution.pdf'),height=8, width=18)

pp <- getDefaultPlotParams(plot.type = 2)                                                       
pp$data1height <- 60                                                                            
pp$data2height <- 60  

kp <- plotKaryotype(genome="hg38", cex=1.6, plot.type=2, plot.params=pp)
kpDataBackground(kp, color = "#FFFFFFAA")

kpPlotRegions(kp, data=bed, col="#EEAAFF", r0=0, r1=1)
kpPlotDensity(kp, data=bed, col='#1FD086', r0=0, r1=1)

kpPlotRegions(kp, data=ctbed, col="#A8A7A8", r0=0, r1=1, data.panel = 2)
kpPlotDensity(kp, data=ctbed, col='#ACA58B', r0=0, r1=1, data.panel = 2)
#kpPoints(kp, data=bed, data.panel=1, y=bed$percentage, ymax=100, ymin=0, cex=0.3, col='#2DD4FA', border='#265FF2')
#kpAxis(kp, ymax=100, ymin=0, data.panel=1)

dev.off()

#png(paste0(bedfile,'_cov.png'), width=2400, height=3000, res=150)

#kp <- plotKaryotype(genome="hg38", cex=1.6, plot.type=2)
#kpDataBackground(kp, color = "#FFFFFFAA")


#dev.off()
