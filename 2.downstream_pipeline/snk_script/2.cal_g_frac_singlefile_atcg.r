library(dplyr)
library(readr)
library(ggseqlogo)
library(ggplot2)
library(ggprism)
library(argparse)
library(stringr)
library(tidyverse)
library(ggsci)
parser <- ArgumentParser()

parser$add_argument('-d','--workdir')

args <- parser$parse_args()
thedir <- args$workdir
thelist <- read_tsv(paste0(thedir,'/list.txt'), col_names=c('file'))
thelist <- thelist$file[!grepl('inpeak', thelist$file)]

dftreat <- list()
dfctrl <- list()

dfggtreat <- data.frame(label = c(), ptrn = c(), con = c())
dfggctrl <- data.frame(label = c(), ptrn = c(), con = c())


for(dfname in thelist){
	tryCatch(
	expr = {
		dff <- read_tsv(paste0(thedir, '/' ,dfname), col_names=c('site','seq'))
		if(nrow(dff) <= 1){
			next
		}
		label <- str_extract(dfname, "^.*(?=_stop)")  
			# dff$label <- str_extract(dfname, "^.*(?=_stop)")  
		con <- str_extract(dfname, "(ctrl|treat)")
		if(con == 'treat'){
			dftreat[[label]] <- dff$seq |> as.vector()
			dfggtreat <- dfggtreat |> bind_rows(
												data.frame(
														label = label,
														ptrn = substr(dff$seq, 10, 11),
														con = 'treat'

														)
												)
			
		} else {
			dfctrl[[label]] <- dff$seq |> as.vector()
			dfggctrl <- dfggctrl |> bind_rows(
												data.frame(
														label = label,
														ptrn = substr(dff$seq, 10, 11),
														con = 'ctrl'
														)
												)
		}
	},
	error = function(e){
		print(paste0("Error!", dfname))
	}
    )
}

tryCatch(
	expr = {
		figtreat <<- ggseqlogo(
								lapply(dftreat, toupper),
								ncol = 4,
								seq_type = 'dna',
								method = 'prob'
							) +
			annotate('rect', xmin=9.5, xmax=11.5, ymin=0, ymax=1, col='black', fill='yellow', alpha=.1) +
			scale_x_continuous(breaks=seq(1,21,2),labels=seq(-10,10,2))

	},
	error = function(e){
		figtreat <<- ggplot()
	}
)

tryCatch(
	expr = {
		figctrl <<- ggseqlogo(
								lapply(dfctrl, toupper),
								ncol = 4,
								seq_type = 'dna',
								method = 'prob'
							) +
			annotate('rect', xmin=9.5, xmax=11.5, ymin=0, ymax=1, col='black', fill='yellow', alpha=.1) +
			scale_x_continuous(breaks=seq(1,21,2),labels=seq(-10,10,2))
	},
	error = function(e){
		figctrl <<- ggplot()
	}
)

dfgg <- bind_rows(dfggtreat, dfggctrl)
dfgg$ptrn <- toupper(dfgg$ptrn)

dfgg$ptrn <- gsub(
                  "(?!AG|GG)[ACGTN]{2}",
                  "Others",
                  dfgg$ptrn,
                  perl=TRUE
                  )

labels <- dfgg$label |> unique()

dfgg_nonctrl <- dfgg |> filter(con != 'ctrl')

tryCatch(
	expr={
		dfgg_ctrl1 <<- dfgg |> filter(label %in% c(labels[1]), con == 'ctrl')
		dfgg_ctrl1$label <- 'ctrl_fwd' #cannot be <<-
	},
	error = function(e){
		dfgg_ctrl1 <<- data.frame(label = c(), ptrn = c(), con = c())
	}
)

tryCatch(
		 expr = {
			dfgg_ctrl2 <<- dfgg |> filter(label %in% c(labels[2]), con == 'ctrl')
			dfgg_ctrl2$label <- 'ctrl_rvs' #cannot be <<-
		 },
		 error = function(e){
		dfgg_ctrl1 <<- data.frame(label = c(), ptrn = c(), con = c())
		 }
)



dfgg_data <- bind_rows(dfgg_nonctrl, dfgg_ctrl1) |> bind_rows(dfgg_ctrl2)

#get gg fraction

ggcount <- dfgg_data |> group_by(label, ptrn) |> summarise(count=n())

total_counts <- ggcount |> group_by(label) %>%
  summarise(total_count = sum(count))

gg_fraction <- ggcount %>%
  filter(ptrn == "GG") %>%
  inner_join(total_counts, by = "label") %>%
    mutate(fraction = count / total_count)

write.table(gg_fraction, paste0(thedir,'/gg_fraction.txt'), quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

figgg <- ggplot(dfgg_data, aes(x=factor(label, levels=c('ctrl_fwd','ctrl_rvs',labels)), fill=factor(ptrn, levels=c('Others','AG','GG')))) +
    geom_bar(stat='count', position='fill', color='black') + 
    geom_text(data=gg_fraction, aes(x=factor(label, levels=c('ctrl_fwd','ctrl_rvs',labels)), y=1.05, label=round(fraction,2)), angle=45, vjust=1) +
    scale_fill_manual(values=c("GG"="#008397", "AG"="#E1DDCA", "Others"="#D8D8D8")) +
    theme_prism(axis_text_angle=45) +
    labs(x="", y='')

ggsave(filename=paste0(thedir,'/treat_seqlogo.pdf'), plot = figtreat, height=4.5,width=16)
ggsave(filename=paste0(thedir,'/ctrl_seqlogo.pdf'), plot = figctrl, height=4.5, width=16)
ggsave(filename=paste0(thedir,'/pattern_fraction.pdf'), plot = figgg, height=9, width=16)
