# for next week, for each sample, find the number of genes with zero counts and compare 
# the numbers to the 95th pctl (e.g. scatterplot, #genes with zero counts vs 95th percentile, 
# one point per sample)


options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))


setwd("~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/comp4.3_tert8.ckcc.outlier_results")

outlierResults<-lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types=cols()) %>%
    add_column(sampleID=gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()

outlierResults_sample <- outlierResults %>%
  select(sample)

# find the number of genes with zero counts
# number of samples that are 0 and not zero
count(outlierResults_sample, sample == 0)
#   `sample == 0`       n
#   <lgl>           <int>
#   1 F             4124131
#   2 T             4428695

# get dataframe overall to compute for each sample file
eachSampleZeroCount <- outlierResults %>%
  group_by(sampleID) %>%
  count(sample) %>%
  filter(sample == 0) 


# the numbers to the 95th pctl
eachSample95thCount <- outlierResults %>%
  group_by(sampleID) %>%
  count(sample) %>%
  filter(sample >= quantile(sample, c(0.95))) %>%
  tally(n)
  

# 4428695 is the number of samples with 0 count

# scatterplot
ggplot(eachSampleZeroCount, aes(sampleID, n)) + 
  geom_point() +
  theme(axis.text.x = element_blank()) 
  
#genes with zero counts vs 95th percentile, 
# one point per sample)

createBoundedData <- function(col1, col2, nameOfComparison, name1, name2) {
  ## very useful tutorial for plotting two histograms together
  # https://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r/3557042
  # p95df <- data.frame(percentileOfEachSampleDf$p95)
  # p95df <- rename(p95df, 'value' = percentileOfEachSampleDf.p95)
  # 
  # p75df <- data.frame(percentileOfEachSampleDf$p75) 
  # p75df <- rename(p75df, 'value' = percentileOfEachSampleDf.p75)
  # 
  # p95df$percentile <- '95th'
  # p75df$percentile <- '75th'
  # 
  # bothPercentiles <- rbind(p95df, p75df)
  
  col1df <- data.frame(col1)
  col1df <- rename(col1df, 'value' = col1)
  
  
  col2df <- data.frame(col2)
  col2df <- rename(col2df, 'value' = col2)
  
  col1df$nameOfComparison <- name1
  col2df$nameOfComparison <- name2
  
  bothBoundedData <- rbind(col1df, col2df)
  colnames (bothBoundedData) <- c("values", nameOfComparison)
  return(bothBoundedData)
}

dfCombined95thZeroCount <- createBoundedData(eachSampleZeroCount$n, eachSample95thCount$nn, "count", 'zero', '95th')
dfCombined95thZeroCount$sampleID <- c(eachSampleZeroCount$sampleID, eachSample95thCount$sampleID)



ggplot(dfCombined95thZeroCount, aes(values, sampleID, color = count)) + 
  geom_point() + 
  ggtitle("95th and Zero Counts of Each Sample") 
