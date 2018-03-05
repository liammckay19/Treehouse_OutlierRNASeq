# just plot genes in the top 5th of each samples, consider the shape of those distnrutions

# top5thGenesPerSamples.R


options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)


liamsWorkingDir <-
  "~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/"

setwd(paste0(liamsWorkingDir, "comp4.3_tert8.ckcc.outlier_results"))

up_outlier_files=list.files(, "outlier_results_")

outlierResults<-lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types=cols()) %>%
    add_column(sampleID=gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()


dfTopFiveGenes <- outlierResults %>%
	filter(sample > quantile(sample, 0.95))


sampleList <- dfTopFiveGenes %>%
	select(sampleID,sample) %>%
	group_by(sampleID) %>%
	summarize()

for (thisSample in sampleList$sampleID) {
	print(thisSample)
	dfi <- dfTopFiveGenes %>% filter(sampleID == thisSample)

	maxGene <- max(dfi$sample)
	varSample <- var(dfi$sample)
	maxVarGene<-dfi[which.max(dfi$sample),]$Gene
	variationOfMax <- round(maxGene- mean(dfi$sample),3)*1000

	p <- ggplot(dfTopFiveGenes %>% filter(sampleID == thisSample)) +
		geom_histogram(aes(sample), binwidth = 0.1) +
		ggtitle(
			paste0("mGene:",
			maxVarGene,
			" var:",
			varSample,
			" ",
			thisSample
		)) +
		xlab("log2(TPM+1)") + 
		ylab("Gene Expression") +
		scale_x_continuous(limits = c(5,15)) +
		scale_y_continuous(limits = c(0,300))

	ggsave(
		paste0("Var_",
			varSample,
			"-",
			thisSample,
			".png"
		),
		plot = p,
		"png",
		paste0(liamsWorkingDir, "Batch-top5-sample-histograms")
	)

}

