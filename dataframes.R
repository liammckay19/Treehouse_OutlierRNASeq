# dataframes.r
# ------------------------------------------------------------------------------
# dataframes
# ------------------------------------------------------------------------------

# Log2(TPM+1) data set
setwd(
  "~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/comp4.3_tert8.ckcc.outlier_results"
)

up_outlier_files = list.files(, "outlier_results_")

outlierResults <- lapply(up_outlier_files, function(x) {
  read_tsv(x, col_types = cols()) %>%
    add_column(sampleID = gsub("outlier_results_", "", x))
}) 	%>%
  bind_rows()

outlierResults_sample <- outlierResults %>%
  select(sample)
dfZerosOrNotZeros <- outlierResults %>%
  select(sampleID, sample) %>%
  group_by(sampleID) %>%
  count(sample == 0) 

dfZeros <- dfZerosOrNotZeros %>%
  group_by(sampleID) %>%
  filter(`sample == 0` == T)

p95df <- outlierResults %>% group_by(sampleID) %>% summarize(p95 = quantile(sample, 0.95))

dfZeros$p95 = p95df$p95

dfScatter <- dfZeros

singleSample <- outlierResults %>%
  select(sampleID, sample, Gene) %>%
  filter(sampleID == "TH01_0053_S01")




# TPM Data set
setwd("~/Documents/UCSC/Junior/Treehouse/OutlierRNAseq_Treehouse_Repo/ckcc_rsem_genes_results")

sample_file_list=list.files(, "TH0")

rawTPMDf<-lapply(sample_file_list, function(x) {
      read_tsv(x, col_types=cols()) %>%
      add_column(sampleID=gsub("outlier_results_", "", x))
      }) 	%>%
      bind_rows()

percentileOfEachTPMDf <- rawTPMDf %>%
  group_by(sampleID) %>%
  summarise(p95 = quantile(TPM, c(0.95)), p75 = quantile(TPM, c(0.75))) %>%
  arrange(desc(p95))

TPMDf <- rawTPMDf %>%
  select(TPM)

no_zero_TPMDf <- TPMDf %>%
  filter(TPM > 0)
summary(TPMDf)
  

Log2TPMDf <- log2(TPMDf+1)
countLog2TPMDf <- count(Log2TPMDf, vars = TPM)


TPM053_sample <- subset(TPMsampleIDDf, sampleID == "TH01_0053_S01_rsem_genes.results")

sum(TPMDf)

sum_TPMDf <- rawTPMDf %>%
  group_by(sampleID) %>%
  summarise(sum = sum(TPM))
# = 1.46e+08 total TPM

countOfTPMDf <- count(TPMDf, vars = TPM)

count(TPMDf, TPM > 0)
# `TPM > 0`       n
# <lgl>       <int>
#   1 F         4690960
#   2 T         4141748


head(countOfTPMDf)
# confirms above works  (the vars = 0; n = 4690960)

sum_TPMDf

# back up reference within 1 of one million for each sample
min(sum_TPMDf$sum)
# [1] 999997.7
max(sum_TPMDf$sum)
# [1] 1000001
