# holly's code
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



nfp<-outlierResults %>% group_by(sampleID) %>% summarize(nfp=quantile(sample, 0.95))

fifteenth=quantile(nfp$nfp, 0.15)
worst15pctSamples<-nfp %>% filter(nfp<fifteenth)

thisSample <- NULL
for (thisSample in worst15pctSamples$sampleID){
  print(thisSample)
  df <- outlierResults %>% filter(sampleID == thisSample) 
  
  p<-ggplot(df) +
    geom_histogram(aes(sample), binwidth = 0.1) +
    ggtitle(thisSample)

  ggsave(paste0(thisSample, ".png"), plot = p, "png", "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/WorstPercentilePlots/")

}
