
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

singleSample <- outlierResults %>%
  filter(sampleID == "TH01_0054_S01")

ggplot(singleSample, aes(sample)) + 
  geom_histogram(binwidth = 0.1, aes(y = ..count..)) + 
  ggtitle("Histogram of Sample TH01_0054_S01; binwidth = 0.1") +
  scale_x_continuous(name = "Log2(TPM+1) value") +
  scale_y_continuous(name = "Count (thousands)",
                       labels = function(y) y / 1000)
