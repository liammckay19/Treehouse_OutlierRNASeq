# ==================================================================================================
# 
# 		ALL SCATTER PLOTS FROM PRESENTATION
# 
# ==================================================================================================
# LOAD DATA FROM comp4.3_tert8.ckcc.outlier_results
# ==================================================================================================
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

# ==================================================================================================
# Number of Expressed Genes vs. p95
# ==================================================================================================
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

dfNotZeros$TH01 <- gsub("TRUE", "RiboMinus",dfNotZeros$TH01)
dfNotZeros$TH01 <- gsub("FALSE", "PolyA",dfNotZeros$TH01)


dfScatter <- dfNotZeros 
colnames(dfScatter)[which(names(dfScatter) == "TH01")] <- "Sample_Name"
dfScatter$Has_Bump <- allSamples[[3]]
dfScatter$Bump_Value  <- allSamples[[4]]
ggplot(dfScatter,aes(n/1000,p95, color = Sample_Name)) + 
  scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
  scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
  geom_point(aes(shape=Has_Bump)) +
  scale_size_continuous(range = c(1,6))+
  xlab('Expressed Genes (Thousands)') + ylab('95th Percentile per Sample') +
  ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile') +
  geom_smooth(method = 'lm')+
  annotate(
    "text",
    x = 30,
    y = 1 ,
    label = paste0(
      "correlation TH01: ",
      round(cor(dfTH01s$n,dfTH01s$p95),3),
      "\ncorrelation Not TH01: ",
      round(cor(dfNotTH01s$n,dfNotTH01s$p95),3)
      
    )
  )

dfScatter$shortSampleID <- gsub('[_][0-9S]+','',dfScatter$sampleID)
i<-1
corlist<-list(0)
for(sampleCenter in dfScatter$shortSampleID){
	print(sampleCenter)
	dfab <- dfScatter %>% filter(shortSampleID == sampleCenter)
	corlist[[i]] <- round(cor(dfab$n, dfab$p95),4)
	i <- i+ 1
}

correlations<-data.frame(sampleID = unique(dfScatter$shortSampleID))
correlations$cor <- unique(corlist)

for (i in seq(1,4)){
	paste0()
}


resultsCorrelations <- paste((paste0(correlations$sampleID, ": ", correlations$cor, "; \n")), collapse = '')
ggplot(dfScatter,aes(n/1000,p95, color = shortSampleID)) + 
  scale_colour_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  geom_point() +
  scale_size_continuous(range = c(1,6))+
  geom_smooth(method = 'lm')+
  annotate(
    "text",
    x = 30,
    y = 2.5 ,
    label = paste0("correlations\n",resultsCorrelations)
    )+
  xlab('Expressed Genes (Thousands)') + ylab('95th Percentile per Sample') +
  ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile') 





# ==================================================================================================
# UMEND count vs p95
# ==================================================================================================
# LOAD DATA FROM ckccSampleUMEND_reads and comp4.3_tert8.ckcc.outlier_results
# ==================================================================================================

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





# ==================================================================================================
# TPM vs. p95
# ==================================================================================================
options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))


setwd("~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/ckcc_rsem_genes_results")

sample_file_list=list.files(, "TH0")

rawTPMDf<-lapply(sample_file_list, function(x) {
		read_tsv(x, col_types=cols()) %>%
		add_column(sampleID=gsub("_rsem_genes", "", x))
		})  %>%
		bind_rows()

percentileOfEachTPMSampleDf <- rawTPMDf %>%
	group_by(sampleID) %>% 
	summarize(p95 = quantile(TPM, 0.95))

dfNotZerosOrNotZeros <- rawTPMDf %>%
	select(sampleID, TPM) %>%
	group_by(sampleID) %>%
	count(TPM == 0) 

dfNotZeros <- dfNotZerosOrNotZeros %>%
	group_by(sampleID) %>%
	filter(`TPM == 0` == F)


dfNotZeros$p95<-percentileOfEachTPMSampleDf$p95
dfNotZeros$zeroCount <- dfNotZeros$n

dfScatter <- dfNotZeros 
ggplot(dfScatter, aes(p95, zeroCount/1000)) + 
	scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
	scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
	geom_point() +
	ylab('Expressed Genes (Thousands)') + xlab('95th Percentile per Sample (TPM)') +
	ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile (TPM)') +
	geom_smooth(method = 'lm')+
	annotate(
	"text",
	x = 30,
	y = 20 ,
	label = paste0(
		"correlation : ",
		round(cor(dfScatter$n,dfScatter$p95),3)
		
	)
	)
dfScatter$shortSampleID <- gsub('.results','',dfScatter$sampleID)

dfScatter$shortSampleID <- gsub('[_][0-9S]+','',dfScatter$shortSampleID)
dfScatter$Method <- gsub('[TH02-9]+[^TH01]', 'PolyA', dfScatter$shortSampleID)
dfScatter$Method <- gsub('TH01', 'RiboD', dfScatter$Method)
i<-1
corlist<-list(0)
for(sampleCenter in dfScatter$shortSampleID){
  print(sampleCenter)
  dfab <- dfScatter %>% filter(shortSampleID == sampleCenter)
  corlist[[i]] <- round(cor(dfab$n, dfab$p95),4)
  i <- i+ 1
}

correlations<-data.frame(sampleID = unique(dfScatter$shortSampleID))
correlations$cor <- unique(corlist)
resultsCorrelations <- paste((paste0(correlations$sampleID, ": ", correlations$cor, "; \n")), collapse = '')

i<-1
corlistM<-list(0)
for(method in unique(dfScatter$Method)){
  print(method)
  dfab <- dfScatter %>% filter(Method == method)
  corlistM[[i]] <- round(cor(dfab$n, dfab$p95),4)
  i <- i+ 1
}

correlationsM<-data.frame(method = unique(dfScatter$Method))
correlationsM$corMethod <- corlistM
resultsCorrelationsM <- paste((paste0(correlationsM$method, ": ", correlationsM$corMethod, "; \n")), collapse = '')


ggplot(dfScatter, aes(p95, zeroCount/1000, color=shortSampleID)) + 
	scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
	scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
	geom_point() +
	ylab('Expressed Genes (Thousands)') + xlab('95th Percentile per Sample (TPM)') +
	ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile (TPM)') +
	geom_smooth(method = 'lm')+
	annotate("text",x = 15,y = 35 ,label = paste0(	"correlation : \n",	resultsCorrelations	))
ggplot(dfScatter, aes(p95, zeroCount/1000, color=Method)) + 
	scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
	scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
	geom_point() +
	ylab('Expressed Genes (Thousands)') + xlab('95th Percentile per Sample (TPM)') +
	ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile (TPM)') +
	geom_smooth(method = 'lm')+
	annotate("text",x = 15,y = 35 ,label = paste0(	"correlation : \n",	resultsCorrelationsM	))


# ==================================================================================================
# UMEND vs. Expressed Genes per Sample
# ==================================================================================================
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
  geom_smooth(method = 'lm', level = 0) +
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




# ==================================================================================================
# Correlation of 1
# ==================================================================================================

# REFERENCE
x <- seq(1,5)
y <- x
df = data.frame(x,y)

cor(x,y)
ggplot(df, aes(x,y)) + geom_point() + ggtitle("Correlation of 1")
# correlation of 1 (when x = y)


