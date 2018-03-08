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


rawTPMDf %>% filter(Gene == mostVariableGenes)

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
dfVarSamples <- rawTPMDf %>% group_by(sampleID) %>% filter(Gene %in% geneList$Gene)
x <- list()
for (thisSample in worst15pctSamples$sampleID) {
	print(thisSample)
	dfi <- rawTPMDf %>% filter(sampleID == paste0(thisSample,'.results'))
	mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
	genes <- dfi$gene_id
	dfi<-dfi[,-4]
	G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)
	merge(dfi,G_list,by.x="gene",by.y="ensembl_peptide_id")

	
	p<- ggplot(dfi) +
		geom_histogram(aes(expected_count), binwidth=1)+
		ggtitle(thisSample)
	ggsave(filename = paste0("Expected-count-",thisSample,".png"), p,
       width = 2, height = 2, dpi = 150, units = "in", device='png', paste0(liamsWorkingDir, "Batch-expectedCount"))

}