# plot # Expressed genes vs UMEND reads per sample;
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


# plot # Expressed genes vs UMEND reads per sample;

setwd(paste0(liamsWorkingDir, "comp4.3_tert8.ckcc.outlier_results"))

up_outlier_files = list.files(, "outlier_results_")

outlierResults <- lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types = cols()) %>%
    add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()

dfNonZeros <- outlierResults %>%
  select(sampleID, sample) %>%
  group_by(sampleID) %>%
  count(sample == 0) %>%
  filter(`sample == 0` == F)


umendNonZerosDF <- dfNonZeros %>%
  arrange(sampleID) %>%
  filter(sampleID %in% umendResults$th_sampleid) %>%
  add_column(rawUMEND = umendResults$umendCountRaw, umendID = umendResults$th_sampleid)

ggplot(umendNonZerosDF, aes(n/1000, rawUMEND)) + geom_point() +
  xlab("Thousands of Expressed Genes Per Sample") + ylab("Raw UMEND Count") +
  geom_smooth(method = 'lm') +
  ggtitle("Number of Expressed Genes verses UMEND Reads Per Sample") +
  annotate(
    "text",
    x = 23.5,
    y = 7e+07,
    label = paste0(
      "correlation: ",
      round(cor(umendNonZerosDF$rawUMEND,umendNonZerosDF$n),3)
    )
  )
