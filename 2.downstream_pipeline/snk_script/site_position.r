library(karyoploteR)
library(readr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('-b','--bed', help='Site bed files.')
args <- parser$parse_args()
bedfile <- args$bed

bed_df <- read.table(bedfile, sep = "\t", header = FALSE, 
    col.names = c("chromosome", "start", "end", "seq",'ptrn','strand','sample','cut','rate','type'))


df_tr_GG <- bed_df |> filter(ptrn=='GG', type == 'treat')
df_tr_nGG <- bed_df |> filter(ptrn!='GG', type == 'treat')
df_ct_GG <- bed_df |> filter(ptrn=='GG', type == 'ctrl')
df_ct_nGG <- bed_df |> filter(ptrn!='GG', type == 'ctrl')


bed_tr_GG <- GRanges(
    seqnames = df_tr_GG$chromosome,
    ranges = IRanges(start = df_tr_GG$start, end = df_tr_GG$end)
)

bed_tr_nGG <- GRanges(
    seqnames = df_tr_nGG$chromosome,
    ranges = IRanges(start = df_tr_nGG$start, end = df_tr_nGG$end)
)

bed_ct_GG <- GRanges(
    seqnames = df_ct_GG$chromosome,
    ranges = IRanges(start = df_ct_GG$start, end = df_ct_GG$end)
)

bed_ct_nGG <- GRanges(
    seqnames = df_ct_nGG$chromosome,
    ranges = IRanges(start = df_ct_nGG$start, end = df_ct_nGG$end)
)

pdf(paste0(bedfile,'_distribution.pdf'),height=19, width=9)

pp <- getDefaultPlotParams(plot.type = 2)
pp$data1height <- 60
pp$data2height <- 60
kp <- plotKaryotype(genome="hg38", plot.type=2, plot.params=pp)
kpDataBackground(kp, color = "#FFFFFFAA")

kpPlotRegions(kp, data=bed_tr_GG, col="#EEAAFF", r0=0, r1=0.5)
kpPlotDensity(kp, data=bed_tr_GG, col='#1FD086', r0=0, r1=0.5)
kpPlotRegions(kp, data=bed_tr_nGG, col="#EEAAFF", r0=0.5, r1=1)
kpPlotDensity(kp, data=bed_tr_nGG, col='#365EE9', r0=0.5, r1=1)

kpPlotRegions(kp, data=bed_ct_GG, col="#EEAAFF", r0=0, r1=0.5, data.panel=2)
kpPlotDensity(kp, data=bed_ct_GG, col='#1FD086', r0=0, r1=0.5, data.panel=2)
kpPlotRegions(kp, data=bed_ct_nGG, col="#EEAAFF", r0=0.5, r1=1, data.panel=2)
kpPlotDensity(kp, data=bed_ct_nGG, col='#365EE9', r0=0.5, r1=1, data.panel=2)

#kpPoints(kp, data=bed, data.panel=1, y=bed$percentage, ymax=100, ymin=0, cex=0.3, col='#2DD4FA', border='#265FF2')
#kpAxis(kp, ymax=100, ymin=0, data.panel=1)

dev.off()

#png(paste0(bedfile,'_cov.png'), width=2400, height=3000, res=150)

#kp <- plotKaryotype(genome="hg38", cex=1.6, plot.type=2)
#kpDataBackground(kp, color = "#FFFFFFAA")


#dev.off()
