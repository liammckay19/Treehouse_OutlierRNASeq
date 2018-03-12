
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

UMENDp95DF$shortSampleID <- gsub('[_][0-9S]+','',UMENDp95DF$sampleID)
i<-1
corlist<-list(0)
for(sampleCenter in UMENDp95DF$shortSampleID){
  print(sampleCenter)
  dfab <- UMENDp95DF %>% filter(shortSampleID == sampleCenter)
  corlist[[i]] <- round(cor(dfab$rawUMEND, dfab$p95),4)
  i <- i+ 1
}

correlations<-data.frame(sampleID = unique(UMENDp95DF$shortSampleID))
correlations$cor <- unique(corlist)
resultsCorrelations <- paste((paste0(correlations$sampleID, ": ", correlations$cor, "; \n")), collapse = '')


ggplot(UMENDp95DF, aes(p95, rawUMEND, color = shortSampleID)) + geom_point() +
  ylab("Raw UMEND Count") + xlab("95th Percentile Per Sample") + 
  geom_smooth(method = 'lm', level = 0.95, se = TRUE) +
  ggtitle("95th Percentile vs. UMEND Reads Per Sample") + 
  annotate(
    "text",
    x = 3,
    y = 7e+07,
    label = paste0("correlations:\n", resultsCorrelations)
      
    )
  )


