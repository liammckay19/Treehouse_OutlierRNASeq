# Unecessary code at this point feb8/2018

# creates 95th vs 75th percentiles of all samples in one histogram
# creates 95th percentile of all samples in a histogram

options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))
import::here(createBoundedData, .from = "createBoundedData.R")

setwd("~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/comp4.3_tert8.ckcc.outlier_results")

up_outlier_files=list.files(, "outlier_results_")
outlierResults<-lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types=cols()) %>%
    add_column(sampleID=gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()
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


