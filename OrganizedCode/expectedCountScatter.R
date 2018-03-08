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

liamsWorkingDir <-
	"~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/"

setwd(paste0(liamsWorkingDir, "comp4.3_tert8.ckcc.outlier_results"))

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


liamsWorkingDir <-
	"~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/"

setwd(paste0(liamsWorkingDir))

mostVariableGenes <- read_tsv("JVB TH_log2tpm_variable_genes.txt", col_types = cols(), col_names=c("Gene")) 

percentileOfEachTPMSampleDf <- rawTPMDf %>%
	group_by(sampleID) %>% 
	summarize(p95q = quantile(TPM, 0.95), p95e_c = quantile(expected_count, 0.95))

conversionEnsmblGeneID<- read_tsv("gene-DAVID-EnsembleID-Convert.txt", col_types = cols(), col_names=TRUE)

rawTPM

rawTPMDf %>% gsub(".[0-9]","",gene_id) %>% filter(gene_id %in% conversionEnsmblGeneID$To)

rawTPMDf <- gsub("[0-9].[0-9]", "[0-9]",rawTPMDf$gene_id)

thisSample <- NULL
order <- 0
sumBump <- 0
countBump <- 0 
averageBump <- 0 


# https://stackoverflow.com/questions/28543517/how-can-i-convert-ensembl-id-to-gene-symbol-in-r
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library('biomaRt')


ggplot(percentileOfEachTPMSampleDf, aes(p95q, p95e_c/1000)) + geom_point() + geom_smooth(method = 'lm')  +
	ggtitle("95th Pctl Expected Count of Samples vs. 95th Pctl TPM of Samples") +
	xlab("Sample's 95th Percentile of TPM Values") + ylab("Sample's 95th Percentile of Expected Count (Thousands)")+
	annotate(
		"text",
		x = 30,
		y = 10,
		label = paste0(
			"correlation: ",
			round(cor(percentileOfEachTPMSampleDf$p95q, percentileOfEachTPMSampleDf$p95e_c),3)

		))

# NOOO lookup gene id in tpm df and use tpm to get 1000 most var genes

	# dfVarTPM <- rawTPMDf %>% mutate(TPMlog2 = log2(TPM+1))

	# samp <- data.frame(outlierResults %>% select(sample, sampleID))

	# dfVarTPM <- left_join(dfVarTPM, samp, by="sampleID")

	# dfSamples %>% arrange((global95))

	# dfSamples$TH01 <- grepl("TH01", dfSamples$sampleID)

	# # dfVarSamples <- dfVarTPM %>% group_by(sampleID) %>% filter(TPMlog2 %in% outlierResults$sample)
	# dfCollected<- data.frame(dfVarTPM,outlierResults$sample)

#  variance
dfGeneVar <- rawTPMDf %>%
	group_by(gene_id) %>%
	summarize(variation = var(TPM)) %>%
	arrange(desc(variation))

geneList <- dfGeneVar %>% filter(variation > quantile(dfGeneVar$variation, 0.95))
# get names of genes p95 of variation and up

dfPercentile <- rawTPMDf %>%
	group_by(gene_id) %>%
	summarize(p95 = quantile(TPM, 0.95))

dfSamples <-
		rawTPMDf %>% group_by(gene_id) %>% filter(gene_id %in% geneList$gene_id)
	# match names to all of their ESN001, ESN0032 etc...


x <- list()
for (thisSample in worst15pctSamples$sampleID) {
	print(thisSample)
	dfi <- rawTPMDf %>% filter(sampleID == paste0(thisSample,'.results'))


	p<-ggplot(dfi,aes(expected_count/1000)) +
		geom_histogram(binwidth=1)+
		ggtitle(thisSample) +
		scale_x_continuous(limits=c(0,5000))+
		scale_y_continuous(limits=c(0,100)) +
		xlab("Expected Count (Thousands)") + ylab("Frequency")

	ggsave(filename = paste0("Expected-count-",thisSample,".png"), p, device='png', paste0(liamsWorkingDir, "Batch-expectedCount"))

}