# scatterplot 

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

dfZerosOrNotZeros <- outlierResults %>%
  select(sampleID, sample) %>%
  group_by(sampleID) %>%
  count(sample == 0) 

dfZeros <- dfZerosOrNotZeros %>%
  group_by(sampleID) %>%
  filter(`sample == 0` == T)

p95df <- outlierResults %>% group_by(sampleID) %>% summarize(p95 = quantile(sample, 0.95))

dfZeros$p95 = p95df$p95

dfScatter <- dfZeros
ggplot(dfScatter, aes(n/1000, p95)) + geom_point() + 
  xlab('Unexpressed Genes per Thousand') + ylab('95th Percentile of Sample') +
  ggtitle('Each Sample\'s Count of Unexpressed Genes and its 95th Percentile') +
  geom_smooth(method = 'lm')


