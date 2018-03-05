# most variable genes.r


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

dfGeneMean <- outlierResults %>%
	group_by(Gene) %>%
	summarize(mean = mean(sample))

mean(dfGeneMean$mean)
# the overall mean is 1.019 sample 

dfGeneVar <- outlierResults %>%
	group_by(Gene) %>%
	summarize(variation = var(sample)) %>%
	arrange(desc(variation))

mean(dfGeneVar$variation) # 0.56
# so most genes differ from the norm on average by 0.56 

sd(dfGeneVar$variation) # standard deviation = 1.15

summary(dfGeneVar)
 #     Gene             variation        
 # Length:58581       Min.   : 0.000000  
 # Class :character   1st Qu.: 0.001579  
 # Mode  :character   Median : 0.069611  
 #                    Mean   : 0.565726  
 #                    3rd Qu.: 0.707318  
 #                    Max.   :28.579563  


quantile(dfGeneVar$variation, 0.95)
# > 95% of the data variates from the mean by 2.56



geneList <- dfGeneVar %>% filter(variation > quantile(dfGeneVar$variation, 0.95))
# get names of genes p95 of variation and up


dfSamples <-
		outlierResults %>% group_by(sampleID) %>% filter(Gene %in% geneList$Gene)
	# match names to all of their th01 th02 etc...

sampleList <- dfSamples %>%
	select(sampleID,sample) %>%
	group_by(sampleID) %>%
	summarize()

for (thisSample in sampleList$sampleID) {
	print(thisSample)
	dfi <- dfSamples %>% filter(sampleID == thisSample)
	p <- ggplot(dfSamples %>% filter(sampleID == thisSample)) +
		geom_histogram(aes(sample), binwidth = 0.1) +
		ggtitle(thisSample) +
		xlab("log2(TPM+1)") + ylab("Gene Expression")
	maxGene <- max(dfi$sample)
	maxVarGene<-dfi[which.max(dfi$sample),]$Gene
	variationOfMax <- round(maxGene- mean(dfi$sample),3) 
	ggsave(
		paste0(
			maxVarGene,
			"-",
			variationOfMax,
			"-",
			thisSample,
			".png"
		),
		plot = p,
		"png",
		paste0(liamsWorkingDir, "MostVariantGenes")
	)

}


# saves plots and bumps for 22 sample files on the low end of the 95th percentile
{
	nfpDF <-
		dfSamples %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
	
	fifteenth = quantile(nfpDF$nfp, 0.15)
	worst15pctSamples <-
		nfpDF %>% filter(nfp < fifteenth) %>% arrange(desc(nfp))
	
	thisSample <- NULL
	order <- 0
	sumBump <- 0
	countBump <- 0 
	averageBump <- 0 

	x <- list()
	for (thisSample in worst15pctSamples$sampleID) {
		print(thisSample)
		df <-
			dfSamples %>% filter(sampleID == thisSample) %>% filter(sample > 1.8)
		dfn <- count(df, sample=round(sample,1))
		dfn$index <- seq(1,length(dfn$n))

		maxGene <- max(df$sample)
		maxVarGene<-df[which.max(df$sample),]$Gene
		variationOfMax <- round(maxGene- mean(df$sample),3) 

		order = order + 1
		p <- ggplot(dfSamples %>% filter(sampleID == thisSample)) +
			geom_histogram(aes(sample), binwidth = 0.1) +
			ggtitle(paste0(thisSample,"\n maxGene: ", maxVarGene)) +
			scale_x_continuous(limits = c(0,15))+
			scale_y_continuous(limits = c(0,125))+
			xlab("log2(TPM+1)") + ylab("Gene Expression")

		if(dfn[which.max(dfn$n),]$sample > 2.1) {
			p = p + annotate(
					"text",
					x = dfn[which.max(dfn$n),]$sample+3,
					y = 1000,
					label = paste0(
						"bump: ",
						dfn[which.max(dfn$n),]$sample

					)
				) + geom_vline(xintercept = dfn[which.max(dfn$n),]$sample)
			sumBump = sumBump + dfn[which.max(dfn$n),]$sample
			countBump = countBump + 1
		}
	

		ggsave(
			paste0(
				order + 100,
				"_",
				round(worst15pctSamples[[2]][order], digits = 3),
				"_",
				thisSample,
				".png"
			),
			plot = p,
			"png",
			paste0(liamsWorkingDir, "BatchPlotsMostVar-Below-p15")
		)
		
	}
	averageBump = sumBump / countBump
	print(averageBump)
}


# saves plots and bumps for 22 sample files on the high end of the 95th percentile
{
	nfpDF <-
		dfSamples %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
	
	eightyFive = quantile(nfpDF$nfp, 0.85)
	best85pctSamples <-
		nfpDF %>% filter(nfp > eightyFive) %>% arrange(desc(nfp))
	
	thisSample <- NULL
	order <- 0
	sumBump <- 0
	countBump <- 0 
	averageBump <- 0 

	x <- list()
	for (thisSample in best95pctSamples$sampleID) {
		print(thisSample)
		df <-
			dfSamples %>% filter(sampleID == thisSample) %>% filter(sample > 1.8)
		dfn <- count(df, sample=round(sample,1))
		dfn$index <- seq(1,length(dfn$n))

		maxGene <- max(df$sample)
		maxVarGene<-df[which.max(df$sample),]$Gene
		variationOfMax <- round(maxGene- mean(df$sample),3) 

		order = order + 1
		p <- ggplot(dfSamples %>% filter(sampleID == thisSample)) +
			geom_histogram(aes(sample), binwidth = 0.1) +
			ggtitle(paste0(thisSample,"> p85\n maxGene: ", maxVarGene)) +
			scale_x_continuous(limits = c(0,20))+
			scale_y_continuous(limits = c(0,300))+
			xlab("log2(TPM+1)") + ylab("Gene Expression")

		if(dfn[which.max(dfn$n),]$sample > 2.5) {
			p = p + annotate(
					"text",
					x = dfn[which.max(dfn$n),]$sample+3,
					y = 1000,
					label = paste0(
						"bump: ",
						dfn[which.max(dfn$n),]$sample

					)
				) + geom_vline(xintercept = dfn[which.max(dfn$n),]$sample)
			sumBump = sumBump + dfn[which.max(dfn$n),]$sample
			countBump = countBump + 1
		}
	

		ggsave(
			paste0(
				order + 100,
				"_",
				round(best95pctSamples[[2]][order], digits = 3),
				"_",
				thisSample,
				".png"
			),
			plot = p,
			"png",
			paste0(liamsWorkingDir, "BatchPlotsMostVar-Above-p95")
		)
		
	}
	averageBump = sumBump / countBump
	print(averageBump)
}

maxGene <- max(dfSamples$sample)
maxVarGene<-dfSamples[which.max(dfSamples$sample),]$Gene
variationOfMax <- round(maxGene- mean(dfSamples$sample),3) 

ggplot(dfSamples %>% filter(sampleID == best85pctSamples$sampleID), aes(sample)) +
	geom_histogram(binwidth=0.1) +
	ggtitle(paste0("Samples > p85 | maxGene: ", maxVarGene, " | Distance From Mean: ", variationOfMax)) +
	xlab("log2(TPM+1)") + ylab("Gene Expression")+
	facet_wrap(~ sampleID)
