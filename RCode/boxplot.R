# Boxplot
# do ribo depleted samples generally have more expressed genes? 
# make a boxplot and do a statistics (t?)

options(stringsAsFactors=FALSE) # for compatibile code between us

library(tidyverse)
library(gridExtra) # easy for putting graphs onto the same page (just use ggarrange(graph1, graph2, ncol = # of display
# columns, nrow = #row))


setwd("~/Documents/UCSC/Junior/Treehouse/Treehouse_outlierRNASeq/comp4.3_tert8.ckcc.outlier_results")

up_outlier_files=list.files(, "outlier_results_")

outlierResults<-lapply(up_outlier_files, function(x) {
	read_tsv(x, col_types=cols()) %>%
		add_column(sampleID=gsub("outlier_results_", "", x))
}) 	%>%
	bind_rows()

dfZerosOrNotZeros <- outlierResults %>%
	select(sampleID, sample) %>%
	group_by(sampleID) %>%
	count(sample == 0) 

dfNotZeros <- dfZerosOrNotZeros %>%
	group_by(sampleID) %>%
	filter(`sample == 0` == F)

p95df <- outlierResults %>% group_by(sampleID) %>% summarize(p95 = quantile(sample, 0.95))

dfNotZeros$p95 = p95df$p95

dfBox <- dfNotZeros

dfBox$Method = grepl("TH01", p95df$sampleID)
dfBox$Method <- gsub("TRUE", "Ribo-depletion(TH01)",dfBox$Method)
dfBox$Method <- gsub("FALSE", "PolyA-selection(Not TH01)",dfBox$Method)
riboDepleted <- filter(dfBox, Method=="Ribo-depletion(TH01)")
polyaSelected <- filter(dfBox, Method=="PolyA-selection(Not TH01)")

ggplot(dfBox, aes(x=Method, y=p95)) + geom_boxplot()+
	ylab('95th Percentile Per Sample') + xlab('Method') +
	ggtitle('Ribo-depletion and PolyA-selection 95th Percentile Values') 

t.test(riboDepleted$p95, polyaSelected$p95,
       alternative = "two.sided",
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)

# Welch Two Sample t-test

# data:  riboDepleted$p95 and polyaSelected$p95
# t = -5.3986, df = 24.118, p-value = 1.496e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -1.2779641 -0.5712024
# sample estimates:
# mean of x mean of y 
#  4.374965  5.299548 

ggplot(dfBox, aes(x=Method, y=n/1000)) + geom_boxplot()+
	ylab('Expressed Genes (Thousands)') + xlab('Method') +
	ggtitle('Ribo-depletion and PolyA-selection Measured Expression') +
	annotate("text",x = -1,y = 30 , hjust =0, label = paste0(  "Welch Two Sample t-test
		data:  riboDepleted$n and polyaSelected$n
		t = 5.121, df = 49.389, p-value = 5.016e-06
		alternative hypothesis: 
			true difference in means is not equal to 0
		95 percent confidence interval:
		 	1352.013 3097.878
		sample estimates: 
			mean of x 30121.91  
			mean of y 27896.97 
		"
		)
	
	)

ggplot(dfBox, aes(x=Method, y=n/1000)) + geom_boxplot()+
	ylab('Expressed Genes (Thousands)') + xlab('Method') +
	ggtitle('Ribo-depletion and PolyA-selection Measured Expression') 
	annotate("text",x = -1,y = 30 , hjust =0, label = paste0(  "Welch Two Sample t-test
		data:  riboDepleted$n and polyaSelected$n
		t = 5.121, df = 49.389, p-value = 5.016e-06
		alternative hypothesis: 
			true difference in means is not equal to 0
		95 percent confidence interval:
		 	1352.013 3097.878
		sample estimates: 
			mean of x 30121.91  
			mean of y 27896.97 
		"
		)
	
	)


t.test(riboDepleted$n, polyaSelected$n,
       alternative = "two.sided",
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95)
# 	Welch Two Sample t-test

# data:  riboDepleted$n and polyaSelected$n
# t = 5.121, df = 49.389, p-value = 5.016e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  1352.013 3097.878
# sample estimates:
# mean of x mean of y 
#  30121.91  27896.97 

t.test(riboDepleted$n, polyaSelected$n,
       alternative = "two.sided",
       mu = 0, paired = FALSE, var.equal = TRUE,
       conf.level = 0.95)

# Two Sample t-test

# data:  riboDepleted$n and polyaSelected$n
# t = 3.6244, df = 144, p-value = 0.000401
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  1011.554 3438.337
# sample estimates:
# mean of x mean of y 
#  30121.91  27896.97

t.test(riboDepleted$n, polyaSelected$n,
       alternative = "two.sided",
       mu = 5, paired = FALSE, var.equal = TRUE,
       conf.level = 0.95)




# REFERENCE
x <- seq(1,5)
y <- x
df = data.frame(x,y)

cor(x,y)
ggplot(df, aes(x,y)) + geom_point() + ggtitle("Correlation of 1")
# correlation of 1 (when x = y)