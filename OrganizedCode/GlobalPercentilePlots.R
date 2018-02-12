


options(stringsAsFactors = FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))
# Create plots of 95th vs 50,60,75,80,85-93 percentile

liamsWorkingDir <- "~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/comp4.3_tert8.ckcc.outlier_results"
setwd(
  liamsWorkingDir
)

up_outlier_files = list.files(, "outlier_results_")

outlierResults <- lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types = cols()) %>%
    add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()



percentileOfEachSampleDf <- outlierResults %>%
  group_by(sampleID) %>%
  summarise(
    p95 = quantile(sample, c(0.95)),
    p75 = quantile(sample, c(0.75))
    ,
    p80 = quantile(sample, c(0.80)),
    p60 = quantile(sample, c(0.60))
    ,
    p50 = quantile(sample, c(0.50)),
    p85 = quantile(sample, c(0.85))
    ,
    p86 = quantile(sample, c(0.86)),
    p90 = quantile(sample, c(0.90))
    ,
    p87 = quantile(sample, c(0.87)),
    p91 = quantile(sample, c(0.91))
    ,
    p88 = quantile(sample, c(0.88)),
    p92 = quantile(sample, c(0.92))
    ,
    p89 = quantile(sample, c(0.89)),
    p93 = quantile(sample, c(0.93))
  ) %>%
  arrange(desc(p95))


#dataframes
{
  p95p50 <-
    gather(percentileOfEachSampleDf, percentile, sample, p95, p50)
  p95p60 <-
    gather(percentileOfEachSampleDf, percentile, sample, p95, p60)
  p95p75 <-
    gather(percentileOfEachSampleDf, percentile, sample, p95, p75)
  p95p80 <-
    gather(percentileOfEachSampleDf, percentile, sample, p95, p80)
  p95p85 <-
    gather(percentileOfEachSampleDf, percentile, sample, p95, p85)
  p95p92 <-
    gather(percentileOfEachSampleDf, percentile, sample, p95, p92)
  }
# # Multiple plot of 50-92

{
  nfp = list(p95p50, p95p60, p95p75, p95p85, p95p92)
  names = list("50th and 95th",
               "60th and 95th",
               "75th and 95th",
               "85th and 95th",
               "92nd and 95th")
  
  for (i in seq(1, 5)) {
    p1 <- ggplot(nfp[[i]], aes(sample, fill = percentile)) +
      geom_histogram(
        alpha = 0.5,
        position = 'identity',
        aes(y = ..density.., color = percentile),
        binwidth = 0.1
      ) +
      ggtitle(names[i]) +
      geom_density(alpha = 0,
                   position = 'identity',
                   aes(y = ..density..)) +
      xlab("log2(TPM+1) Pctl Cutoff \nof Each Sample") +
      ylab("Density of Measurements") +
      scale_y_continuous(limits = c(0, 1.5)) +
      scale_x_continuous(limits = c(0, 6.5))
    ggsave(
      paste0(names[i], "pctlPlot.png"),
      plot = p1,
      "png",
      "~/Documents/UCSC/Junior/Treehouse/Treehouse_OutlierRNASeq/GlobalPercentilePlots/",
    )
    
  }
}


ggplot(p95p75, aes(sample, fill = percentile)) +
  geom_histogram(
    alpha = 0.5,
    position = 'identity',
    aes(y = ..density.., color = percentile),
    binwidth = 0.1
  ) +
  ggtitle("Frequency of Global 75th and 95th Percentile of Each Sample") +
  geom_density(alpha = 0,
               position = 'identity',
               aes(y = ..density..)) +
  xlab("log2(TPM+1) Pctl Cutoff \nof Each Sample") +
  ylab("Density of Measurements") +
  scale_y_continuous(limits = c(0, 1.5)) +
  scale_x_continuous(limits = c(0, 6.5)) +
  geom_vline(xintercept = 5.24) + annotate("text", x = 4.4, y = 1, label = "log2(TPM+1) = 5.24")

ggplot(percentileOfEachSampleDf, aes(p95)) +
  geom_histogram(
    alpha = 0.5,
    position = 'identity',
    aes(y = ..density..),
    binwidth = 0.1
  ) +
  ggtitle("Frequency of Global 95th Percentile in Each Sample") +
  geom_density(alpha = 0,
               position = 'identity',
               aes(y = ..density..)) +
  xlab("95th Percentile") +
  ylab("Density of Measurements")


m <- ggplot(percentileOfEachSampleDf, aes(p95))
m <- m + geom_density()
p <- print(m)

head(p$data[[1]]$y)
sum(p$data[[1]]$y)
# sum of data is 173.9702
# I want to divide the points of the density by the sum of all to get a normalized curve that adds to one

normalized <- data.frame(p$data[[1]]$y,p$data[[1]]$x)
sum(p$data[[1]]$y/sum(p$data[[1]]$y))                                
# this sum is 1

ggplot(normalized, aes(p$data[[1]]$x, p$data[[1]]$y/sum(p$data[[1]]$y))) + geom_line() +
  ggtitle("Normalized Density Curve Sum All y = 1") +
  ylab("Density of Observations") +
  xlab("log2(TPM+1)")
