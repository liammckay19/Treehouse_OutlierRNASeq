# Plots.r
# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------


ggplot(dfScatter, aes(n, p95)) + geom_point() + 
  xlab('Unexpressed Genes') + ylab('95th Percentile of Sample') +
  ggtitle('Each Sample\'s Count of Unexpressed Genes and its 95th Percentile') +
  geom_smooth(method = 'lm')

ggplot(TPM053_sample, aes(log2(TPM+1))) +geom_histogram()
ggplot(countLog2TPMDf, aes(vars, n)) + geom_point() # weird binning was observed
# it does spread in a different way

ggplot(TPMDf, aes(TPM)) + # with zeros
  geom_histogram() # plot of bell curve centered at 3 

ggplot(log2(no_zero_TPMDf-1), aes(TPM)) + # no zeros 
  geom_histogram() + 
  geom_vline(xintercept = 0) # plot of bell curve centered at 3 

# binwidth of 1 for 95th percentile of TPM across all genes
ggplot(percentileOfEachTPMDf, aes(p95)) + geom_histogram(binwidth = 1)

# binwidth of 0.1 for 75th percentile of TPM across all genes
ggplot(percentileOfEachTPMDf, aes(p75)) + geom_histogram(binwidth = 0.1)


ggplot(singleSample, aes(sample)) + geom_histogram(binwidth = 0.1) + 
  ggtitle("Histogram of Sample TH01_0053_S01; binwidth = 0.1")


# make file batch of plots
nfp <-
	outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))

fifteenth = quantile(nfp$nfp, 0.15)
worst15pctSamples <- nfp %>% filter(nfp < fifteenth)

thisSample <- NULL
for (thisSample in worst15pctSamples$sampleID) {
	print(thisSample)
	df <- outlierResults %>% filter(sampleID == thisSample)
	
	p <- ggplot(df) +
		geom_histogram(aes(sample), binwidth = 0.1) +
		ggtitle(thisSample)
	
	ggsave(
		paste0(thisSample, ".png"),
		plot = p,
		"png",
		"~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/"
	)
	
}
