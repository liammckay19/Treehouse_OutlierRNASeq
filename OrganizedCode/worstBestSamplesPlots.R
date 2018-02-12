
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

worstSamplesNames <- percentileOfEachSampleDf %>%
  select(p95, sampleID) %>%
  filter(p95 <= quantile(p95, c(0.15))) %>%
  mutate(names = gsub("T", "outlier_results_T", sampleID))

worstSamples <- lapply(worstSamplesNames$names, function(x) {
  read_tsv(x, col_types = cols()) %>%
    add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()

ggplot(worstSamples, aes(sample, color = sampleID)) +
  geom_histogram(binwidth = 0.1) +
  ggtitle("Bump In Worst Samples") +
  ylim(0, 2e+04)
# interesting but hard to read graph of all the worst samples in a total histogram



# saves plots for all sample files on the low end of the 95th percentile
{
  nfpDF <-
    outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
  
  fifteenth = quantile(nfpDF$nfp, 0.15)
  worst15pctSamples <-
    nfpDF %>% filter(nfp < fifteenth) %>% arrange(desc(nfp))
  
  thisSample <- NULL
  order <- 0
  for (thisSample in worst15pctSamples$sampleID) {
    print(thisSample)
    df <- outlierResults %>% filter(sampleID == thisSample)
    order = order + 1
    p <- ggplot(df) +
      geom_histogram(aes(sample), binwidth = 0.1) +
      ggtitle(thisSample) +
      scale_x_continuous(limits = c(0, 20)) +
      scale_y_continuous(limits = c(0, 2000))
    
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
      paste0(liamsWorkingDir, "WorstPercentilePlots")
    )
    
  }
}




# saves plots for all sample files on the high end of the 95th percentile
{
  nfpDF <-
    outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
  
  p85 = quantile(nfpDF$nfp, 0.85)
  
  best85pctSamples <-
    nfpDF %>% filter(nfp > p85) %>% arrange(desc(nfp))
  order <- 0
  thisSample <- NULL
  for (thisSample in best85pctSamples$sampleID) {
    print(thisSample)
    thisSample <- best85pctSamples[[1]][1]
    df <- outlierResults %>% filter(sampleID == thisSample)
    
    
    p <- ggplot(df) +
      geom_histogram(aes(sample), binwidth = 0.1) +
      ggtitle(thisSample) +
      scale_x_continuous(limits = c(0, 20)) +
      scale_y_continuous(limits = c(0, 2000)) 
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
}
