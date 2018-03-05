
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



UMENDp95DF$Sample = grepl("TH01", UMENDp95DF$umendID)

UMENDp95DFTH01 <- UMENDp95DF %>% filter(Sample == T)
UMENDp95DFOther <- UMENDp95DF %>% filter(Sample == F)

UMENDp95DF$Sample <- gsub("TRUE", "TH01_[Ribo]...",UMENDp95DF$Sample)
UMENDp95DF$Sample <- gsub("FALSE", "Not TH01_[PolyA]...",UMENDp95DF$Sample)

ggplot(UMENDp95DF, aes(p95, rawUMEND, color = Sample)) + geom_point() +
  ylab("Raw UMEND Count") + xlab("95th Percentile Per Sample") + 
  geom_smooth(method = 'lm', level = 0.95, se = TRUE) +
  ggtitle("[Ribo-D / PolyA-S] 95th Percentile vs. UMEND Reads Per Sample") + 
  annotate(
    "text",
    x = 3,
    y = 7e+07,
    label = paste0(
      "correlation TH01 [Ribo-D]: ",
      round(cor(UMENDp95DFTH01$rawUMEND,UMENDp95DFTH01$p95),3),
      "\ncorrelation Not TH01 [PolyA-S]: ",
      round(cor(UMENDp95DFOther$rawUMEND,UMENDp95DFOther$p95),3)
      
    )
  )


