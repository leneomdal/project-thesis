source("load_data.R")
library(svglite)
#Explore data

#Path to store plots
path_plots = "plots"

mirna.means = rowMeans(ten.days.sample.df)
mirna.means[1]
plot(mirna.means)
hist(mirna.means)

#Function for plotting histogram of expression levels in one miRNA
hist_expression_level = function(dataframe, mirna){
  hist(as.numeric(dataframe[mirna,]))
}

#Same function using ggplot
hist_expression_level_gg = function(dataframe, mirna){
  dataframe = as.data.frame(t(dataframe))
  col = sym(colnames(dataframe)[mirna])
  p = ggplot(data = dataframe, aes(x = !!col, y = ..density..)) + geom_histogram()
  return(p)
}
hist_expression_level2(ten.days.sample.df, 1)



# BOXPLOT OF MOST ABUNDANT miRNAs

sum.reads = apply(ten.days.sample.df, 1, sum)

sorted.sum.reads = sum.reads[order(sum.reads, decreasing = TRUE)]

n.mirna.to.include = 20


most.abundant.mirna.names = names(sorted.sum.reads[seq(1, n.mirna.to.include)])

most.abundant.mirnas = ten.days.sample.df[most.abundant.mirna.names,]


abundant.mirna.long = pivot_longer(data = as.data.frame(t(most.abundant.mirnas)), cols = rownames(most.abundant.mirnas), names_to = "miRNA", values_to = "cpmValue")




abundant.mirna.long$miRNA = str_remove(abundant.mirna.long$miRNA, "hsa-")
abundant.mirna.long$miRNA = factor(abundant.mirna.long$miRNA, levels = str_remove(most.abundant.mirna.names, "hsa-"))




ggplot(data = abundant.mirna.long, aes(x = miRNA, y = cpmValue )) + geom_boxplot() + ylab("cpm value") + ggtitle("20 most abundant miRNAs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major.x = element_blank()) 
ggsave("abundant-miRNAs-boxplot.svg", path = path_plots, height = 4.5, width = 5)


# SEPARATED BOXPLOT

abundant.mirna.long$probiotic.bool = rep(ten.days.meta.data$probiotic, each = n.mirna.to.include)

abundant.mirna.long$Probiotic = ifelse(abundant.mirna.long$probiotic.bool == 1, "Yes", "No" )




ggplot(data = abundant.mirna.long, aes(x = miRNA, y = cpmValue, fill = Probiotic)) + geom_boxplot() + 
  ylab("cpm value") + ggtitle("20 most abundant miRNAs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major.x = element_blank())
ggsave("separate-abundant-miRNAs-boxplot.svg", path = path_plots, height = 4.5, width = 6)
