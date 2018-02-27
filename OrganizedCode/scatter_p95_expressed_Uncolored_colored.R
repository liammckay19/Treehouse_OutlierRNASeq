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


# every sample in the dataset's bump plotted and saved: 
{
	nfpDF <-
		outlierResults %>% group_by(sampleID) %>% summarize(nfp = quantile(sample, 0.95))
	
	p85 = quantile(nfpDF$nfp, 0.85)
	sumBump <- 0
	countBump <- 0
	averageBump <- 0
	allSamples <-
		nfpDF %>% arrange(desc(nfp))
	order <- 0
	thisSample <- NULL
	allSamples$bump <- TRUE
	allSamples$bumpVal <- 0
	i <- 0
	for (thisSample in allSamples$sampleID) {
		print(thisSample)
		df <- outlierResults %>% filter(sampleID == thisSample) %>% filter(sample > 1.8)
		dfn <- count(df, sample=round(sample,1))
		dfn$index <- seq(1,length(dfn$n))

		# p <- ggplot(outlierResults %>% filter(sampleID == thisSample)) +
		# 	geom_histogram(aes(sample), binwidth = 0.1) +
		# 	ggtitle(thisSample) +
		# 	scale_x_continuous(limits = c(0, 20)) +
		# 	scale_y_continuous(limits = c(0, 2000)) 


		if(dfn[which.max(dfn$n),]$sample > 1.9) { # gets the highest count
			allSamples[[3]][i] <- TRUE
			allSamples[[4]][i] <- dfn[which.max(dfn$n),]$sample
			# p = p + annotate(
			# 		"text",
			# 		x = dfn[which.max(dfn$n),]$sample+3,
			# 		y = 1000,
			# 		label = paste0(
			# 			"bump: ",
			# 			dfn[which.max(dfn$n),]$sample
			# 		)
			# 	) + geom_vline(xintercept = dfn[which.max(dfn$n),]$sample)
			sumBump = sumBump + dfn[which.max(dfn$n),]$sample
			countBump = countBump + 1
		} else {
			allSamples[[3]][i] <- FALSE
			allSamples[[4]][i] <- 0

		}
		i <- i + 1
	}
	averageBump = sumBump / countBump
	print(averageBump)
}
	# average bump = 3.344928

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

dfNotZeros$TH01 <- gsub("TRUE", "TH01_...",dfNotZeros$TH01)
dfNotZeros$TH01 <- gsub("FALSE", "Not TH01_...",dfNotZeros$TH01)


dfScatter <- dfNotZeros 
colnames(dfScatter)[which(names(dfScatter) == "TH01")] <- "Sample_Name"
dfScatter$Has_Bump <- allSamples[[3]]
dfScatter$Bump_Value  <- allSamples[[4]]
ggplot(dfScatter,aes(n/1000,p95, color = Sample_Name)) + 
  scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
  scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
  geom_point(aes(shape=Has_Bump, alpha=Bump_Value)) +
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


ggplot(dfScatter,aes(n/1000,p95)) + 
  scale_colour_brewer(type = "seq", palette = "Set1", direction = 1)+
  scale_fill_brewer(type = "seq", palette = "Set1", direction = 1)+
  geom_point() +
  scale_size_continuous(range = c(1,6))+
  geom_smooth(method = 'lm')+
  annotate(
    "text",
    x = 30,
    y = 2.5 ,
    label = paste0(
      "correlation : ",
      round(cor(dfScatter$n,dfScatter$p95),3)
      
    )
  ) +
  ylab('Expressed Genes (Thousands)') + xlab('95th Percentile per Sample') +
  ggtitle('Each Sample\'s Count of Expressed Genes and its 95th Percentile') 



