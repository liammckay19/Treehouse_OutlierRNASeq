# most variable genes.r


options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)


getwd() <-
  "~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/"

setwd(paste0(getwd(), "comp4.3_tert8.ckcc.outlier_results"))

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

dfPercentile <- outlierResults %>%
	group_by(sampleID) %>%
	summarize(p95 = quantile(sample, 0.95))

dfSamples <-
		outlierResults %>% group_by(sampleID) %>% filter(Gene %in% geneList$Gene)
	# match names to all of their th01 th02 etc...

dfSamples$TH01 <- grepl(pattern = 'TH01', dfSamples$sampleID)
dfSamples$TH01 <- gsub('TRUE', 'blue',dfSamples$TH01)
dfSamples$TH01 <- gsub('FALSE', 'red',dfSamples$TH01)



sampleList <- dfSamples %>%
	select(sampleID,sample) %>%
	group_by(sampleID) %>%
	summarize()


for (thisSample in sampleList$sampleID) {
	print(thisSample)
	dfi <- dfSamples %>% filter(sampleID == thisSample)
	dfp <- dfPercentile %>% filter(sampleID == thisSample)
	dfc <- dfSamples %>% group_by(TH01) %>% filter(sampleID == thisSample) %>% summarize()
	print(dfc[1][[1]])
	ggplot(dfSamples %>% filter(sampleID == thisSample)) +
		geom_histogram(aes(sample), binwidth = 0.1) +
		scale_fill_manual(values = c('red')) +
		ggtitle(thisSample) +
		xlab("log2(TPM+1)") + ylab("Gene Expression") +
		scale_x_continuous(limits = c(0,20)) +
		scale_y_continuous(limits = c(0,200)) +
		geom_vline(xintercept = dfp$p95) +
		annotate(
					"text",
					x = round(dfp$p95,4)+4.3,
					y = 150,
					label = paste0(
						"glp95: ",
						round(dfp$p95,4)
					)
				) 

	maxGene <- max(dfi$sample)
	maxVarGene<-dfi[which.max(dfi$sample),]$Gene
	variationOfMax <- round(maxGene- mean(dfi$sample),3) 
	ggsave(
		paste0(
			"pctl=",
			format(round(dfp$p95,4),nsmall = 4),
			"-var=",
			variationOfMax,
			"-",
			thisSample,
			".png"
		),
		plot = p,
		"png",
		paste0(getwd(), "Batch-MostVariantGenesSorted-by-p95")
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
			paste0(getwd(), "BatchPlotsMostVar-Below-p15")
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
			xlab("log2(TPM+1)") + ylab("Gene Expression")+
			geom_vline(xintercept=)

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
			paste0(getwd(), "BatchPlotsMostVar-Above-p95")
		)
		
	}
	averageBump = sumBump / countBump
	print(averageBump)
}
nfpDF <-
		dfSamples %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
	
eightyFive = quantile(nfpDF$nfp, 0.85)
best85pctSamples <-
	nfpDF %>% filter(nfp > eightyFive) %>% arrange(desc(nfp))
fifteenth = quantile(nfpDF$nfp, 0.15)
worst15pctSamples <-
	nfpDF %>% filter(nfp < fifteenth) %>% arrange(desc(nfp))

maxGene <- max(dfSamples$sample)
maxVarGene<-dfSamples[which.max(dfSamples$sample),]$Gene
variationOfMax <- round(maxGene- mean(dfSamples$sample),3) 

pctl <- data.frame(outlierResults %>% group_by(sampleID) %>% summarize(global95 = quantile(sample, 0.95)))

dfSamples <- left_join(dfSamples, pctl, by="sampleID")

dfSamples %>% arrange((global95))

dfSamples$TH01 <- grepl("TH01", dfSamples$sampleID)

facetBigPlot <- ggplot(dfSamples, aes(sample, fill= TH01)) + geom_histogram(binwidth=0.1) +
	ggtitle(paste0("Samples All | maxGene: ", maxVarGene, " | Distance From Mean: ", variationOfMax)) +
	xlab("log2(TPM+1)") + ylab("Gene Expression")+
	scale_x_continuous(limits = c(0,20)) +
	scale_y_continuous(limits = c(0,100)) +
	scale_fill_brewer(palette = "Set1") +
	facet_wrap(~ global95)


ggsave(filename = "facetWrapColored3.png", facetBigPlot,
       width = 20, height = 20, dpi = 150, units = "in", device='png', paste0(getwd()))

ggplot(dfSamples %>% filter(sampleID == best85pctSamples$sampleID), aes(sample)) +
	geom_histogram(binwidth=0.1) +
	ggtitle(paste0("Highest 22 p95s | maxVarGene: ", maxVarGene, " | Distance From Mean: ", variationOfMax)) +
	xlab("log2(TPM+1)") + ylab("Gene Expression")+	
	scale_x_continuous(limits = c(0,17)) +
	scale_y_continuous(limits = c(0,20)) +
	facet_wrap(~ sampleID)

dfBadSamples <- dfSamples %>% filter(sampleID == worst15pctSamples$sampleID)
maxGene <- max(dfBadSamples$sample)
maxVarGene<-dfBadSamples[which.max(dfBadSamples$sample),]$Gene
variationOfMax <- round(maxGene- mean(dfBadSamples$sample),3) 

ggplot(dfBadSamples, aes(sample)) +
	geom_histogram(binwidth=0.1) +
	ggtitle(paste0("Lowest 22 p95s | maxVarGene: ", maxVarGene, " | Distance From Mean: ", variationOfMax)) +
	xlab("log2(TPM+1)") + ylab("Gene Expression")+
	scale_x_continuous(limits = c(0,17)) +
	scale_y_continuous(limits = c(0,20)) +
	facet_wrap(~ sampleID)
