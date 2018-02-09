
options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))


setwd("~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/comp4.3_tert8.ckcc.outlier_results")

up_outlier_files=list.files(, "outlier_results_")

outlierResults<-lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types=cols()) %>%
    add_column(sampleID=gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()

worstSamplesNames <- percentileOfEachSampleDf %>%
  select(p95, sampleID) %>%
  filter(p95 <= quantile(p95, c(0.15))) %>%
  mutate(names = gsub("T", "outlier_results_T",sampleID))
  
worstSamples<-lapply(worstSamplesNames$names, function(x) {
  read_tsv(x, col_types=cols()) %>%
    add_column(sampleID=gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()

ggplot(worstSamples, aes(sample, color = sampleID)) + 
  geom_histogram(binwidth = 0.1) +
  ggtitle("Bump In Worst Samples") +
  ylim(0, 2e+04)
# interesting graph of all the worst samples in a total histogram

splitDF <- split(worstSamples, worstSamples$sampleID)
plotTest<-data.frame(splitDF[1])
ggplot(plotTest, aes(sample)) + 
  geom_histogram(alpha = 0.5,  position = 'identity', aes(y=..density..)) + 
  ggtitle("95th Percentiles of All Samples and the worst sample") +
  geom_density(alpha = 0, position = 'identity', aes(y=..density..))

#for (this sample) in splitDF$sampleID{

}

for (x in splitDF) {
  plotdfSample <- data.frame(splitDF[x])
  p<-ggplot(plotdfSample, aes(sample)) + 
    geom_histogram() + 
    geom_density(alpha = 0, position = 'identity', aes(y=..density..))
}


# saves plots for all sample files on the low end of the 95th percentile
{
  nfp<-outlierResults %>% group_by(sampleID) %>% summarize(nfp=quantile(sample, 0.95))
  
  fifteenth=quantile(nfp$nfp, 0.15)
  worst15pctSamples<-nfp %>% filter(nfp<fifteenth)
  
  thisSample <- NULL
  for (thisSample in worst15pctSamples$sampleID){
    print(thisSample)
    df <- outlierResults %>% filter(sampleID == thisSample) 
    
    p <- ggplot(df) +
      geom_histogram(aes(sample), binwidth = 0.1) +
      ggtitle(thisSample) +
      scale_x_continuous(limits = c(0, 20)) +
      scale_y_continuous(limits = c(0, 2500))
    
    ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/")
    
  }
}


# saves plots for all sample files on the high end of the 95th percentile
{
nfp<-outlierResults %>% group_by(sampleID) %>% summarize(p95=quantile(sample, 0.95))

p85=quantile(nfp$p95, 0.85)
best85pctSamples<-nfp %>% filter(p95>p85)

thisSample <- NULL
for (thisSample in best85pctSamples$sampleID){
  print(thisSample)
  df <- outlierResults %>% filter(sampleID == thisSample) 
  
  p <- ggplot(df) +
    geom_histogram(aes(sample), binwidth = 0.1) +
    ggtitle(thisSample) +
    scale_x_continuous(limits = c(0, 20)) +
    scale_y_continuous(limits = c(0, 2500))
  
  ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/BestPercentilePlots/")
  
}
}
