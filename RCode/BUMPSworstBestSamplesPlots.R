

# IDEA
# plot the line on the second maximum of the distribution to get numerical data of the bump

options(stringsAsFactors = FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))


liamsWorkingDir <-
	"~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/"

setwd(paste0(liamsWorkingDir, "comp4.3_tert8.ckcc.outlier_results"))

up_outlier_files = list.files(, "outlier_results_")

outlierResults <- lapply(up_outlier_files, function(x) {
	read_tsv(x, col_types = cols()) %>%
		add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
	bind_rows()


# saves plots and bumps for 22 sample files on the low end of the 95th percentile
{
	nfpDF <-
		outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
	
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
			outlierResults %>% filter(sampleID == thisSample) %>% filter(sample > 1.8)
		dfn <- count(df, sample=round(sample,1))
		dfn$index <- seq(1,length(dfn$n))

		order = order + 1
		p <- ggplot(outlierResults %>% filter(sampleID == thisSample)) +
			geom_histogram(aes(sample), binwidth = 0.1) +
			ggtitle(thisSample) +
			scale_x_continuous(limits = c(0, 20)) +
			scale_y_continuous(limits = c(0, 2000)) 

		if(dfn[which.max(dfn$n),]$sample > 1.9) {
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
			paste0(liamsWorkingDir, "BatchPlots")
		)
		
	}
	averageBump = sumBump / countBump
	print(averageBump)
}
	# average bump = 2.638889

ggplot(outlierResults %>% filter(sampleID == worst15pctSamples$sampleID), aes(sample)) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~ sampleID) + 
	scale_y_continuous(limits = c(0,50)) + 
	scale_x_continuous(limits = c(0,10)) +
	ylab("Frequency of Sample") + xlab("log2(TPM+1)") + ggtitle("22 Crummy Sample's TPM Frequency")
ggsave("22 Crummy Sample's TPM Frequency", plot=p1, "png", width = 10, height = 10, path = liamsWorkingDir)
# saves plots and bumps for 22 sample files on the high end of the 95th percentile
{
	nfpDF <-
		outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
	
	p85 = quantile(nfpDF$nfp, 0.85)
	sumBump <- 0
	countBump <- 0
	averageBump <- 0
	best85pctSamples <-
		nfpDF %>% filter(nfp > p85) %>% arrange(desc(nfp))
	order <- 0
	thisSample <- NULL
	for (thisSample in best85pctSamples$sampleID) {
		print(thisSample)
		df <- outlierResults %>% filter(sampleID == thisSample) %>% filter(sample > 1.8)
		dfn <- count(df, sample=round(sample,1))
		dfn$index <- seq(1,length(dfn$n))

		p <- ggplot(outlierResults %>% filter(sampleID == thisSample)) +
			geom_histogram(aes(sample), binwidth = 0.1) +
			ggtitle(thisSample) +
			scale_x_continuous(limits = c(0, 20)) +
			scale_y_continuous(limits = c(0, 2000)) 

		if(dfn[which.max(dfn$n),]$sample > 1.9) {
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

		order = order + 1
		ggsave(
			paste0(
				order + 100,
				"_",
				round(best85pctSamples[[2]][order], digits = 3),
				"_",
				thisSample,
				".png"
			),
			plot = p,
			"png",
			paste0(liamsWorkingDir, "BestPercentilePlots")
		)
	}
	averageBump = sumBump / countBump
	print(averageBump)
}
	# average bump is 4.104545

ggplot(outlierResults %>% filter(sampleID == best85pctSamples$sampleID), aes(sample)) +
	geom_histogram(binwidth=0.1) +
	facet_wrap(~ sampleID) + 
	scale_y_continuous(limits = c(0,50)) + 
	scale_x_continuous(limits = c(0,10)) +
	ylab("Frequency of Sample") + xlab("log2(TPM+1)") + ggtitle("22 Best Sample's TPM Frequency")

# every sample in the dataset's bump plotted and saved: 
{
	nfpDF <-
		outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
	
	p85 = quantile(nfpDF$nfp, 0.85)
	sumBump <- 0
	countBump <- 0
	averageBump <- 0
	allSamples <-
		nfpDF %>% arrange(desc(nfp))
	order <- 0
	thisSample <- NULL
	for (thisSample in allSamples$sampleID) {
		print(thisSample)
		df <- outlierResults %>% filter(sampleID == thisSample) %>% filter(sample > 1.8)
		dfn <- count(df, sample=round(sample,1))
		dfn$index <- seq(1,length(dfn$n))

		p <- ggplot(outlierResults %>% filter(sampleID == thisSample)) +
			geom_histogram(aes(sample), binwidth = 0.1) +
			ggtitle(thisSample) + xlab("log2(TPM+1)")


		if(dfn[which.max(dfn$n),]$sample > 1.9) {
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
	

		order = order + 1
		ggsave(
			paste0(
				order + 1000,
				"_",
				round(allSamples[[2]][order], digits = 3),
				"_",
				thisSample,
				".png"
			),
			plot = p,
			"png",
			paste0(liamsWorkingDir, "BatchUnexpressionPlots")
		)
	}
	averageBump = sumBump / countBump
	print(averageBump)
}
	# average bump = 3.344928