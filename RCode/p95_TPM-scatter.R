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
dfScatter$shortSampleID <- gsub('.results','',dfScatter$sampleID)

dfScatter$shortSampleID <- gsub('[_][0-9S]+','',dfScatter$shortSampleID)
dfScatter$Method <- gsub('[TH02-9]+[^TH01]', 'PolyA', dfScatter$shortSampleID)
dfScatter$Method <- gsub('TH01', 'RiboD', dfScatter$Method)
i<-1
corlist<-list(0)
for(sampleCenter in dfScatter$shortSampleID){
  print(sampleCenter)
  dfab <- dfScatter %>% filter(shortSampleID == sampleCenter)
  corlist[[i]] <- round(cor(dfab$n, dfab$p95),4)
  i <- i+ 1
}

correlations<-data.frame(sampleID = unique(dfScatter$shortSampleID))
correlations$cor <- unique(corlist)
resultsCorrelations <- paste((paste0(correlations$sampleID, ": ", correlations$cor, "; \n")), collapse = '')

i<-1
corlistM<-list(0)
for(method in unique(dfScatter$Method)){
  print(method)
  dfab <- dfScatter %>% filter(Method == method)
  corlistM[[i]] <- round(cor(dfab$n, dfab$p95),4)
  i <- i+ 1
}

correlationsM<-data.frame(method = unique(dfScatter$Method))
correlationsM$corMethod <- corlistM
resultsCorrelationsM <- paste((paste0(correlationsM$method, ": ", correlationsM$corMethod, "; \n")), collapse = '')


ggplot(dfScatter, aes(p95, zeroCount/1000, color=shortSampleID)) + 
	scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
	scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
	geom_point() +
	ylab('Expressed Genes (Thousands)') + xlab('95th Percentile per Sample (TPM)') +
	ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile (TPM)') +
	geom_smooth(method = 'lm')+
	annotate("text",x = 15,y = 35 ,label = paste0(	"correlation : \n",	resultsCorrelations	))
ggplot(dfScatter, aes(p95, zeroCount/1000, color=Method)) + 
	scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
	scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
	geom_point() +
	ylab('Expressed Genes (Thousands)') + xlab('95th Percentile per Sample (TPM)') +
	ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile (TPM)') +
	geom_smooth(method = 'lm')+
	annotate("text",x = 15,y = 35 ,label = paste0(	"correlation : \n",	resultsCorrelationsM	))
