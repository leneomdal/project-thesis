library(tidyverse)
path = "\\\\fil.nice.ntnu.no\\nice\\p489\\miRNA 10 days\\data\\copies\\"

#dir.exists(path)
cpm.data = read.csv(paste0(path, "cpm_countMatrix_new data_subset.txt"), sep = "\t", dec= ",", header = TRUE, row.names = 1)
View(cpm.data)
colnames(cpm.data) = str_remove(colnames(cpm.data), "[X]")
ncol(cpm.data)


#Remove last two collumns?
# this is only for full data set, not for the subset
#cpm.data[,c(ncol(cpm.data)-1,ncol(cpm.data))]
#mirna.data = subset(cpm.data, select = -c(ncol(cpm.data)-1, ncol(cpm.data)))
#View(mirna.data)


# For subset, just run the following
mirna.data = cpm.data


meta.data = read.csv(paste0(path, "Breastmilk_overview_10days.csv"), sep = ";", dec = ",", header = TRUE)
class(meta.data)
View(meta.data)
?data.frame


#extract 10 days samples
ten.days.sample.nr = meta.data$bm_sample_no1[seq(1, length(meta.data$bm_sample_no1)-2)]
ten.days.sample.nr = ten.days.sample.nr[toString(ten.days.sample.nr)]


#Find only the ones included in miRNA sequencing project for 10 day samples, mirna = 1 in overview file

ten.days.sample.included.nr = meta.data$bm_sample_no1[meta.data$mirna1 == 1]
n.samples.10.days = length(ten.days.sample.included.nr)


# Check that all samples are present
for (i in seq(from = 1,to = length(ten.days.sample.included.nr))) {
  if(toString(ten.days.sample.included.nr[i] %in% colnames(mirna.data))){
    print("im good")
  }
  else {print(ten.days.sample.included.nr[i])}
}





# Define data frame containing only the samples we want
ten.days.sample.df =  mirna.data[, as.character(ten.days.sample.included.nr)]

View(ten.days.sample.df)
ncol(ten.days.sample.df)


#check for missing values
apply(is.na(ten.days.sample.df), 2, which)
sum(is.na(ten.days.sample.df))

#data frame for only needed metadata

ten.days.meta.data = meta.data[meta.data$mirna1 == 1,]
View(ten.days.meta.data)



#Explore data
mirna.means = rowMeans(ten.days.sample.df)
mirna.means[1]
plot(mirna.means)
hist(mirna.means)

#Function for plotting histogram of expression levels in one miRNA

class(t(ten.days.sample.df))

hist_expression_level = function(dataframe, mirna){
  hist(as.numeric(dataframe[mirna,]))
}

hist_expression_level2 = function(dataframe, mirna){
  dataframe = as.data.frame(t(dataframe))
  col = sym(colnames(dataframe)[mirna])
  p = ggplot(data = dataframe, aes(x = !!col, y = ..density..)) + geom_histogram()
  return(p)
}
hist_expression_level2(ten.days.sample.df, 1)

# installing limma
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")

