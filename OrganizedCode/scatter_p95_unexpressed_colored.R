# plot 95th percentile per sample vs .# unexpressed genes and color the points by whether the smaple ID starts with "TH01".

# scatterplot 

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

dfZerosOrNotZeros <- outlierResults %>%
  select(sampleID, sample) %>%
  group_by(sampleID) %>%
  count(sample == 0) 

dfZeros <- dfZerosOrNotZeros %>%
  group_by(sampleID) %>%
  filter(`sample == 0` == T)

p95df <- outlierResults %>% group_by(sampleID) %>% summarize(p95 = quantile(sample, 0.95))

dfZeros$p95 = p95df$p95

dfZeros$TH01 = grepl("TH01", p95df$sampleID)
dfScatter <- dfZeros 
ggplot(dfScatter, aes(p95,n/1000)) + geom_point() + 
  ylab('Unexpressed Genes per Thousand') + xlab('95th Percentile per Sample') +
  ggtitle('Each Sample\'s Count of Unexpressed Genes and its 95th Percentile') +
  geom_smooth(method = 'lm')


