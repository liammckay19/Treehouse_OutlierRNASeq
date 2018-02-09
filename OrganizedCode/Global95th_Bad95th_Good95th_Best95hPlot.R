


options(stringsAsFactors = FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))
# Create plots of 95th vs 50,60,75,80,85-93 percentile

setwd(
  "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/comp4.3_tert8.ckcc.outlier_results"
)

createBoundedData <-
  function(col1, col2, nameOfComparison, name1, name2) {
    # input: column_name1, column_name2, nameOfSimilarColumn, group_label1, group_label2
    
    # example bothPercentilesOne <- createBoundedData(percentileOfEachSampleDf$p95,
    # percentileOfEachSampleDf$p75,
    # "percentile", '95th', '75th')
    
    # output: dataframe with two columns with labels for each comparison
    # use ggplot(dataframe, aes(values, fill = nameOfSimilarColumn))
    
    ## created my own function to concatenate two columns and compare them in ggplot with a
    # comparison needed
    # (basically generalizes what is written from the tutorial)
    
    ## very useful tutorial for plotting two histograms together
    # https://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r/3557042
    
    # example:
    # p95df <- data.frame(percentileOfEachSampleDf$p95)
    # p95df <- rename(p95df, 'value' = percentileOfEachSampleDf.p95)
    #
    # p75df <- data.frame(percentileOfEachSampleDf$p75)
    # p75df <- rename(p75df, 'value' = percentileOfEachSampleDf.p75)
    #
    # p95df$percentile <- '95th'
    # p75df$percentile <- '75th'
    #
    # bothPercentiles <- rbind(p95df, p75df)
    
    col1df <- data.frame(col1)
    col1df <- rename(col1df, 'value' = col1)
    
    
    col2df <- data.frame(col2)
    col2df <- rename(col2df, 'value' = col2)
    
    col1df$nameOfComparison <- name1
    col2df$nameOfComparison <- name2
    
    bothBoundedData <- rbind(col1df, col2df)
    colnames (bothBoundedData) <- c("values", nameOfComparison)
    return(bothBoundedData)
  }

up_outlier_files=list.files(, "outlier_results_")

outlierResults<-lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types=cols()) %>%
    add_column(sampleID=gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()


# Single bad sample comparison 95% to global samples

# dataframes
{
  # taking the histogram of a bad sample and comparing it to all samples
  # TH01_0069_S01
  up_outlier_file_bad = list.files(, "outlier_results_TH01_0069_S01")
  
  bad_sample <- lapply(up_outlier_file_bad, function(x) {
    read_tsv(x, col_types = cols()) %>%
      add_column(sampleID = gsub("outlier_results_TH01_0069_S01", "", x))
  }) 	%>%
    bind_rows()
  
  percentileOfEachSampleDf <- outlierResults %>%
    group_by(sampleID) %>%
    summarise( p95 = quantile(sample, c(0.95))) %>%
    arrange(desc(p95))
  
  percentileOfEachSampleDf_bad_sample <- bad_sample %>%
    summarise(p95 = quantile(sample, 0.95))
  
  bothPercentilesOverall_Bad <-
    createBoundedData(
      percentileOfEachSampleDf$p95,
      percentileOfEachSampleDf_bad_sample$p95,
      "sampleFile",
      'all samples',
      'Low outlier'
    )
  
  # taking the histogram of a good sample and comparing it to bad sample
  sample_file_good = list.files(, "outlier_results_TH01_0062_S01")
  
  good_sample <- lapply(sample_file_good, function(x) {
    read_tsv(x, col_types = cols()) %>%
      add_column(sampleID = gsub("outlier_results_TH01_0062_S01", "", x))
  }) 	%>%
    bind_rows()
  
  percentileOfEachSampleDf_good_sample <- good_sample %>%
    summarise(p95 = quantile(sample, c(0.95)))
  
  
  
  bothPercentilesGood_Bad <-
    createBoundedData(
      percentileOfEachSampleDf_good_sample$p95,
      percentileOfEachSampleDf_bad_sample$p95,
      "sampleFile",
      'Mean',
      'Low outlier'
    )
  # taking the histogram of a best sample and comparing it to bad sample
  
  
  sample_file_best = list.files(, "outlier_results_TH03_0008_S01")
  
  best_sample <- lapply(sample_file_best, function(x) {
    read_tsv(x, col_types = cols()) %>%
      add_column(sampleID = gsub("outlier_results_TH03_0008_S01", "", x))
  }) 	%>%
    bind_rows()
  
  percentileOfEachSampleDf_best_sample <- best_sample %>%
    summarise(p95 = quantile(sample, c(0.95)))
  
  
  
  bothPercentilesBest_Bad <-
    createBoundedData(
      percentileOfEachSampleDf_best_sample$p95,
      percentileOfEachSampleDf_bad_sample$p95,
      "sampleFile",
      'Best Outlier',
      'Low Outlier'
    )
  }


ggplot(bothPercentilesOverall_Bad, aes(values, fill = sampleFile)) +
  geom_histogram(alpha = 0.5,  position = 'identity', aes(y = ..density..)) +
  ggtitle("95th Percentiles of All Samples and the worst sample") +
  geom_density(alpha = 0, position = 'identity', aes(y = ..density..))


# Single good sample comparison 95% to bad sample 95%

ggplot(bothPercentilesGood_Bad, aes(values, fill = sampleFile)) +
  geom_histogram(alpha = 0.5,  position = 'identity', aes(y = ..density..)) +
  ggtitle("Expression of mean sample and the worst sample") +
  geom_density(alpha = 0, position = 'identity', aes(y = ..density..))


# Single best sample comparison 95% to bad sample 95%


p1 <-
  ggplot(bothPercentilesBest_Bad, aes(values, fill = sampleFile)) +
  geom_histogram(
    alpha = 0.5,
    position = 'identity',
    aes(y = ..density..),
    binwidth = 0.1
  ) +
  ggtitle("Expression of best Sample and the worst sample") +
  geom_density(alpha = 0, position = 'identity', aes(y = ..density..)) +
  xlab("sample") +
  scale_fill_manual(breaks = bothPercentilesBest_Bad$sampleFile,
                    values = c("#FF2D00", "#03BDC0"))


p2 <-
  ggplot(bothPercentilesGood_Bad, aes(values, fill = sampleFile)) +
  geom_histogram(
    alpha = 0.5,
    position = 'identity',
    aes(y = ..density..),
    binwidth = 0.1
  ) +
  ggtitle("Expression  of mean Sample and the worst sample") +
  geom_density(alpha = 0, position = 'identity', aes(y = ..density..)) +
  xlab("sample") +
  scale_fill_manual(breaks = bothPercentilesGood_Bad$sampleFile,
                    values = c("#03BDC0", "#FF2D00"))

# plot both graphs side by side
grob <- arrangeGrob(p1, p2)
grid.arrange(p1, p2)
