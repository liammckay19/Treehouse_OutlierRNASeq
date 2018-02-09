

# Single bad sample comparison 95% to global samples 95%
# taking the histogram of a bad sample and comparing it to all samples
# TH01_0069_S01
up_outlier_file_bad = list.files(, "outlier_results_TH01_0069_S01")

bad_sample <- lapply(up_outlier_file_bad, function(x) {
  read_tsv(x, col_types = cols()) %>%
    add_column(sampleID = gsub("outlier_results_TH01_0069_S01", "", x))
}) 	%>%
  bind_rows()

percentileOfEachSampleDf_bad_sample <- bad_sample
summarise(p95 = quantile(sample, c(0.95)))

bothPercentilesOverall_Bad <-
  createBoundedData(
    percentileOfEachSampleDf$p95,
    percentileOfEachSampleDf_bad_sample$p95,
    "sampleFile",
    'all samples',
    'Low outlier'
  )

ggplot(bothPercentilesOverall_Bad, aes(values, fill = sampleFile)) +
	geom_histogram(alpha = 0.5,  position = 'identity', aes(y = ..density..)) +
	ggtitle("95th Percentiles of All Samples and the worst sample") +
	geom_density(alpha = 0, position = 'identity', aes(y = ..density..))
