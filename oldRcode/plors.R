# Plots.r
# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------

# 95th 75th pctl
percentileOfEachSampleDf <- outlierResults %>%
  group_by(sampleID) %>%
  summarise(p95 = quantile(sample, c(0.95)),p75 = quantile(sample, c(0.75))) %>%
  arrange(desc(p95))

}
bothPercentilesOne <- createBoundedData(percentileOfEachSampleDf$p95, 
										percentileOfEachSampleDf$p75, 
										"percentile", '95th', '75th')
ggplot(percentileOfEachSampleDf, aes(p95)) + 
  geom_histogram(alpha = 0.5, position = 'identity', aes(y=..density..)) + 
  ggtitle("95th Percentiles of All Samples") +
  geom_density(alpha = 0, position = 'identity', aes(y=..density..)) +
  xlab("sample")

ggplot(bothPercentilesOne, aes(values, fill = percentile)) + 
  geom_histogram(alpha = 0.5, position = 'identity', aes(y=..density..)) + 
  ggtitle("95th and 75th Percentiles of All Samples") +
  geom_density(alpha = 0, position = 'identity', aes(y=..density..)) +
  xlab("sample")


# comparezerosofeachsample

dfCombined95thZeroCount <- createBoundedData(eachSampleZeroCount$n, eachSample95thCount$nn, "count", 'zero', '95th')
dfCombined95thZeroCount$sampleID <- c(eachSampleZeroCount$sampleID, eachSample95thCount$sampleID)



ggplot(dfCombined95thZeroCount, aes(values, sampleID, color = count)) + 
  geom_point() + 
  ggtitle("95th and Zero Counts of Each Sample") 


#oulierResults...
paFifty<- ggplot(bothPercentilesThree, aes(values, fill = percentile)) + 
	geom_histogram(binwidth = 0.1, alpha = 0.5, position = 'identity', aes(y=..density..)) + 
	ggtitle("95th and 50th Percentiles of All Samples") +
	geom_density(alpha = 0, position = 'identity', aes(y=..density..))
pbSixty<- ggplot(bothPercentilesFour, aes(values, fill = percentile)) + 
	geom_histogram(binwidth = 0.1, alpha = 0.5, position = 'identity', aes(y=..density..)) + 
	ggtitle("95th and 60th Percentiles of All Samples") +
	geom_density(alpha = 0, position = 'identity', aes(y=..density..))
  # ---> PLOT 95th 75th percentile histograms of total samples
  pcSeventyFive<- ggplot(bothPercentilesOne, aes(values, fill = percentile)) + 
	geom_histogram(binwidth = 0.1, alpha = 0.5, position = 'identity', aes(y=..density..)) + 
	ggtitle("95th and 75th Percentiles of All Samples") +
	geom_density(alpha = 0, position = 'identity', aes(y=..density..)) +
	xlab("sample")
  
  
  
  pdEighty<- ggplot(bothPercentilesTwo, aes(values, fill = percentile)) + 
	geom_histogram(binwidth = 0.1, alpha = 0.5, position = 'identity', aes(y=..density..)) + 
	ggtitle("95th and 80th Percentiles of All Samples") +
	geom_density(alpha = 0, position = 'identity', aes(y=..density..))
  
  peEightyFive <- ggplot(binwidth = 0.1, bothPercentilesFive, aes(values, fill = percentile)) + 
	geom_histogram(binwidth = 0.1, alpha = 0.5, position = 'identity', aes(y=..density..)) + 
	ggtitle("95th and 85th Percentiles of All Samples") +
	geom_density(alpha = 0, position = 'identity', aes(y=..density..))

grob <- arrangeGrob(grobs = list(paFifty,pbSixty,pcSeventyFive,pdEighty,peEightyFive), nrow=5, ncol=1 )
grid.arrange(grob)
thisSample = "paFifty"
p <- paFifty
ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/")
thisSample = "pbSixty"
p <- pbSixty
ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/")
thisSample = "pcSeventyFive"
p <- pcSeventyFive
ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/")
thisSample = "pdEighty"
p <- pdEighty
ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/")
thisSample = "peEightyFive"
p <- peEightyFive
ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/")

ggplot(bothPercentilesOverall_Bad, aes(values, fill = sampleFile)) + 
    geom_histogram(alpha = 0.5,  position = 'identity', aes(y=..density..)) + 
    ggtitle("95th Percentiles of All Samples and the worst sample") +
    geom_density(alpha = 0, position = 'identity', aes(y=..density..))

ggplot(bothPercentilesGood_Bad, aes(values, fill = sampleFile)) + 
    geom_histogram(alpha = 0.5,  position = 'identity', aes(y=..density..)) + 
    ggtitle("Expression of mean sample and the worst sample") +
    geom_density(alpha = 0, position = 'identity', aes(y=..density..))
  

p1 <- ggplot(bothPercentilesBest_Bad, aes(values, fill = sampleFile)) + 
	geom_histogram(alpha = 0.5,  position = 'identity', aes(y=..density..), binwidth = 0.1) + 
	ggtitle("Expression of best Sample and the worst sample") +
	geom_density(alpha = 0, position = 'identity', aes(y=..density..)) +
	xlab("sample") +
	scale_fill_manual(breaks = bothPercentilesBest_Bad$sampleFile, 
	                  values = c("#FF2D00","#03BDC0")) 


p2 <- ggplot(bothPercentilesGood_Bad, aes(values, fill = sampleFile)) + 
	geom_histogram(alpha = 0.5,  position = 'identity', aes(y=..density..), binwidth = 0.1) + 
	ggtitle("Expression  of mean Sample and the worst sample") +
	geom_density(alpha = 0, position = 'identity', aes(y=..density..))+
	xlab("sample") + 
	scale_fill_manual(breaks = bothPercentilesGood_Bad$sampleFile, values = c("#03BDC0","#FF2D00")) 

# plot both graphs side by side
grob <- arrangeGrob(p1,p2)
grid.arrange(p1,p2) 


