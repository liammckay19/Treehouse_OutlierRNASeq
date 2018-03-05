# p95_TPM_scatter.r

options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))


setwd("~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/ckcc_rsem_genes_results")

sample_file_list=list.files(, "TH0")

rawTPMDf<-lapply(sample_file_list, function(x) {
		read_tsv(x, col_types=cols()) %>%
		add_column(sampleID=gsub("_rsem_genes", "", x))
		})  %>%
		bind_rows()

percentileOfEachTPMSampleDf <- rawTPMDf %>%
	group_by(sampleID) %>% 
	summarize(p95 = quantile(TPM, 0.95))

dfNotZerosOrNotZeros <- rawTPMDf %>%
	select(sampleID, TPM) %>%
	group_by(sampleID) %>%
	count(TPM == 0) 

dfNotZeros <- dfNotZerosOrNotZeros %>%
	group_by(sampleID) %>%
	filter(`TPM == 0` == F)


dfNotZeros$p95<-percentileOfEachTPMSampleDf$p95
dfNotZeros$zeroCount <- dfNotZeros$n

dfScatter <- dfNotZeros 
ggplot(dfScatter, aes(p95, zeroCount/1000)) + 
	scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
	scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
	geom_point() +
	ylab('Expressed Genes (Thousands)') + xlab('95th Percentile per Sample (TPM)') +
	ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile (TPM)') +
	geom_smooth(method = 'lm')+
	annotate(
	"text",
	x = 30,
	y = 20 ,
	label = paste0(
		"correlation : ",
		round(cor(dfScatter$n,dfScatter$p95),3)
		
	)
	)
