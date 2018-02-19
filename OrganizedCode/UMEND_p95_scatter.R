
options(stringsAsFactors = FALSE) # for compatibile code between us

library(tidyverse)


liamsWorkingDir <-
  "~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/"
setwd(liamsWorkingDir)

umend_file = list.files(, "ckccSampleUMEND_reads")

umendResults <- lapply(umend_file, function(x) {
  read_tsv(x, col_types = cols()) %>%
    add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows() %>%
  arrange(th_sampleid)


# ''' plot 95 percentile vs UMEND reads per sample '''

setwd(paste0(liamsWorkingDir, "comp4.3_tert8.ckcc.outlier_results"))

up_outlier_files = list.files(, "outlier_results_")

outlierResults <- lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types = cols()) %>%
    add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()



UMENDp95DF <- outlierResults %>%
  group_by(sampleID) %>%
  summarise(p95 = quantile(sample, c(0.95))) %>%
  arrange(sampleID) %>%
  filter(sampleID %in% umendResults$th_sampleid) %>%
  add_column(rawUMEND = umendResults$umendCountRaw, umendID = umendResults$th_sampleid)

ggplot(UMENDp95DF, aes(p95, rawUMEND)) + geom_point() +
  ylab("Raw UMEND Count") + xlab("95th Percentile Per Sample") + 
  geom_smooth(method = 'lm', level = 0.95, se = TRUE) +
  ggtitle("95th Percentile verses UMEND Reads Per Sample") + 
  annotate("text",
    x = 3,
    y = 7e+07,
    label = paste0("correlation: ",round(cor(UMENDp95DF$rawUMEND,UMENDp95DF$p95),3),"\nconfidence = 0.95"))

