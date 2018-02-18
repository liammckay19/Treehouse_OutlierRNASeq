# plot 95th percentile per sample vs .# Expressed genes and color the points by whether the smaple ID starts with "TH01".

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

dfNotZerosOrNotZeros <- outlierResults %>%
  select(sampleID, sample) %>%
  group_by(sampleID) %>%
  count(sample == 0) 

dfNotZeros <- dfNotZerosOrNotZeros %>%
  group_by(sampleID) %>%
  filter(`sample == 0` == F)

p95df <- outlierResults %>% group_by(sampleID) %>% summarize(p95 = quantile(sample, 0.95))

dfNotZeros$p95 = p95df$p95

dfNotZeros$TH01 = grepl("TH01", p95df$sampleID)
dfTH01s <- dfNotZeros %>% filter(TH01 == T)
dfNotTH01s <- dfNotZeros %>% filter(TH01 == F)
dfScatter <- dfNotZeros 
ggplot(dfScatter,aes(p95,n/1000, fill = TH01)) + geom_point(aes(p95,n/1000, color = TH01)) + 
  ylab('Expressed Genes per Thousand') + xlab('95th Percentile per Sample') +
  ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile') +
  geom_smooth(method = 'lm', aes(color = TH01))+
  annotate(
    "text",
    x = 3,
    y = 35 ,
    label = paste0(
      "correlation TH01: ",
      round(cor(dfTH01s$n,dfTH01s$p95),3),
      "\ncorrelation Not TH01: ",
      round(cor(dfNotTH01s$n,dfNotTH01s$p95),3)
      
    )
  )

