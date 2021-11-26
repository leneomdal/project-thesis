library(limma)

path = "C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\new-data\\"

#Read the data
count.data <- read.csv(paste0(path, "MatureMatrix.csv"),sep="\t", header =T,row.names=1)
sample.sheet<-read.table(paste0(path,"SampleSheet.txt"),sep="\t", header=T)

#Remove X from colnames
colnames(count.data) = str_remove(colnames(count.data), "[X]")

# check for missing values
sum(is.na(sample.sheet) == TRUE)             # Two samples miss value for sens2yrs, row 16 and 64



#Find only the ones included in miRNA sequencing project for 10 day samples, mirna = 1 in overview file
ten.days.sample.included.nr = sample.sheet$bm_sample_no1[sample.sheet$day10 == 1]
n.samples.10.days = length(ten.days.sample.included.nr)


# Check that all samples are present
for (i in seq(from = 1,to = length(ten.days.sample.included.nr))) {
  if(toString(ten.days.sample.included.nr[i] %in% colnames(count.data))){
    print("im good")
  }
  else {print(ten.days.sample.included.nr[i])}
}



# Define data frame containing only the samples we want
ten.days.sample.df =  count.data[, as.character(ten.days.sample.included.nr)]



#check for missing values in these extrated samples
apply(is.na(ten.days.sample.df), 2, which)
sum(is.na(ten.days.sample.df))


#data frame for only needed metadata
ten.days.meta.data = sample.sheet[sample.sheet$day10 ==1,]
#nrow(ten.days.meta.data)

count.df = ten.days.sample.df
metadata.df = ten.days.meta.data
