# expectedCountScatter.r

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


	

setwd(paste0(getwd(), "comp4.3_tert8.ckcc.outlier_results"))

up_outlier_files = list.files(, "outlier_results_")

outlierResults <- lapply(up_outlier_files, function(x) {
	read_tsv(x, col_types = cols()) %>%
		add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
	bind_rows()


nfpDF <-
	outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))

fifteenth = quantile(nfpDF$nfp, 0.15)
worst15pctSamples <-
	nfpDF %>% filter(nfp < fifteenth) %>% arrange(desc(nfp))


getwd() <-
	

setwd(paste0(getwd()))

mostVariableGenes <- read_tsv("JVB TH_log2tpm_variable_genes.txt", col_types = cols(), col_names=c("Gene")) 

percentileOfEachTPMSampleDf <- rawTPMDf %>%
	group_by(sampleID) %>% 
	summarize(p95q = quantile(TPM, 0.95), p95e_c = quantile(expected_count, 0.95))


rawTPMDf %>% gsub(".[0-9]","",gene_id) %>% filter(gene_id %in% conversionEnsmblGeneID$To)

rawTPMDf <- gsub("[0-9].[0-9]", "[0-9]",rawTPMDf$gene_id)

ggplot(percentileOfEachTPMSampleDf, aes(p95q, p95e_c/1000)) + geom_point() + geom_smooth(method = 'lm')  +
	ggtitle("Upper Normalized Measured Read Counts vs. Upper Expected Read Counts") +
	xlab("Sample's 95th Percentile of TPM Values") + ylab("Sample's 95th Percentile of Expected Count (Thousands)")+
	annotate(
		"text",
		x = 30,
		y = 10,
		label = paste0(
			"correlation: ",
			round(cor(percentileOfEachTPMSampleDf$p95q, percentileOfEachTPMSampleDf$p95e_c),3)

		))

percentileOfEachTPMSampleDf$shortSampleID <- gsub('.results','',percentileOfEachTPMSampleDf$sampleID)

percentileOfEachTPMSampleDf$shortSampleID <- gsub('[_][0-9S]+','',percentileOfEachTPMSampleDf$shortSampleID)
percentileOfEachTPMSampleDf$Method <- gsub('[TH02-9]+[^TH01]', 'PolyA', percentileOfEachTPMSampleDf$shortSampleID)
percentileOfEachTPMSampleDf$Method <- gsub('TH01', 'RiboD', percentileOfEachTPMSampleDf$Method)
i<-1
corlist<-list(0)
for(sampleCenter in unique(percentileOfEachTPMSampleDf$shortSampleID)){
  print(sampleCenter)
  dfab <- percentileOfEachTPMSampleDf %>% filter(shortSampleID == sampleCenter)
  corlist[[i]] <- round(cor(dfab$p95q, dfab$p95e_c),4)
  i <- i+ 1
}

correlations<-data.frame(sampleID = unique(percentileOfEachTPMSampleDf$shortSampleID))
correlations$cor <- unique(corlist)
resultsCorrelations <- paste((paste0(correlations$sampleID, ": ", correlations$cor, "; \n")), collapse = '')


ggplot(percentileOfEachTPMSampleDf, aes(p95q, p95e_c/1000, color=shortSampleID)) + geom_point() + geom_smooth(method = 'lm')  +
	ggtitle("95th Pctl Expected Count of Samples vs. 95th Pctl TPM of Samples") +
	xlab("Sample's 95th Percentile of TPM Values") + ylab("Sample's 95th Percentile of Expected Count (Thousands)")+
	annotate(
		"text",
		x = 12,
		y = 8,
		label = paste0(
			"correlation: \n",
			resultsCorrelations

		))


i<-1
corlistM<-list(0)
for(method in unique(percentileOfEachTPMSampleDf$Method)){
  print(method)
  dfab <- percentileOfEachTPMSampleDf %>% filter(Method == method)
  corlistM[[i]] <- round(cor(dfab$p95q, dfab$p95e_c),4)
  i <- i+ 1
}

correlationsM<-data.frame(method = unique(percentileOfEachTPMSampleDf$Method))
correlationsM$corMethod <- corlistM
resultsCorrelationsM <- paste((paste0(correlationsM$method, ": ", correlationsM$corMethod, "; \n")), collapse = '')
ggplot(percentileOfEachTPMSampleDf, aes(p95q, p95e_c/1000, color=Method)) + geom_point() + geom_smooth(method = 'lm')  +
	ggtitle("95th Pctl Expected Count of Samples vs. 95th Pctl TPM of Samples") +
	xlab("Sample's 95th Percentile of TPM Values") + ylab("Sample's 95th Percentile of Expected Count (Thousands)")+
	annotate(
		"text",
		x = 12,
		y = 8,
		label = paste0(
			"correlation: \n",
			resultsCorrelationsM

		))


#  variance
dfGeneVar <- rawTPMDf %>%
	group_by(gene_id) %>%
	summarize(variation = var(TPM)) %>%
	arrange(desc(variation))

geneList <- dfGeneVar %>% filter(variation > quantile(dfGeneVar$variation, 0.95))
# get names of genes p95 of variation and up

dfPercentile <- rawTPMDf %>%
	group_by(sampleID) %>%
	summarize(p95 = quantile(TPM, 0.95)) %>%
	arrange(desc(p95))

dfSamples <-
		rawTPMDf %>% group_by(gene_id) %>% filter(gene_id %in% geneList$gene_id)
	# match names to all of their ESN001, ESN0032 etc...


x <- list()
for (thisSample in dfPercentile$sampleID) {
	print(thisSample)
	dfi <- rawTPMDf %>% filter(sampleID == thisSample)

	percentileDf <- dfPercentile %>% filter(sampleID == thisSample)
	p<-ggplot(dfi,aes(expected_count/1000)) +
		geom_histogram(binwidth=0.1)+
		ggtitle(paste0(thisSample, " p95: (TPM) ", round(percentileDf$p95,4))) +
		scale_y_continuous(limits = c(0,100)) +
		scale_x_continuous(limits = c(0,200)) +
		xlab("Expected Count (Thousands)") + ylab("Frequency")

	ggsave(filename = paste0("Expected-count-",round(percentileDf$p95,4),thisSample,".png"), p, device='png', paste0(getwd(), "Batch-expectedCount"))

}