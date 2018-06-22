# createBoundedData.R
# input: column_name1, column_name2, nameOfSimilarColumn, group_label1, group_label2

# example bothPercentilesOne <- createBoundedData(percentileOfEachSampleDf$p95, 
                                        # percentileOfEachSampleDf$p75, 
                                        # "percentile", '95th', '75th')
                                        
# output: dataframe with two columns with labels for each comparison
# use ggplot(dataframe, aes(values, fill = nameOfSimilarColumn))

## created my own function to concatenate two columns and compare them in ggplot with a 
# comparison needed 
# (basically generalizes what is written from the tutorial)
createBoundedData <- function(col1, col2, nameOfComparison, name1, name2) {
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