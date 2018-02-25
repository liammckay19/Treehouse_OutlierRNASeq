# scatterplot 

options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))


setwd("~/Documents/UCSC/Junior/Treehouse/Treehouse_outlierRNASeq/comp4.3_tert8.ckcc.outlier_results")

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
  filter(`sample == 0` == F)

p95df <- outlierResults %>% group_by(sampleID) %>% summarize(p95 = quantile(sample, 0.95))

dfZeros$p95 = p95df$p95

dfScatter <- dfZeros
ggplot(dfScatter, aes( p95,n/1000)) + geom_point() + 
  ylab('Expressed Genes per Thousand') + xlab('95th Percentile of Sample') +
  ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile') +
  geom_smooth(method = 'lm')+
  annotate(
    "text",
    x = 3,
    y = 35 ,
    label = paste0(
      "correlation: ",
      round(cor(dfScatter$n,dfScatter$p95),3)
    )
  )

cor(dfScatter$n,dfScatter$p95)
# 0.08755389
# about 8.8% correlation

# REFERENCE
x <- seq(1,5)
y <- x
df = data.frame(x,y)

cor(x,y)
ggplot(df, aes(x,y)) + geom_point() + ggtitle("Correlation of 1")
# correlation of 1 (when x = y)