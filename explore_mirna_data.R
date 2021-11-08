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



#boxplot of most abundant miRNAs


sum.reads = apply(ten.days.sample.df, 1, sum)

sorted.sum.reads = sum.reads[order(sum.reads, decreasing = TRUE)]

n.mirna.to.include = 20


most.abundant.mirna.names = names(sorted.sum.reads[seq(1, n.mirna.to.include)])

most.abundant.mirnas = ten.days.sample.df[most.abundant.mirna.names,]
class(as.data.frame(t(most.abundant.mirnas)))

abundant.mirna.long = pivot_longer(data = as.data.frame(t(most.abundant.mirnas)), cols = rownames(most.abundant.mirnas), names_to = "miRNA", values_to = "cpmValue")

abundant.mirna.long$miRNA = factor(abundant.mirna.long$miRNA, levels = most.abundant.mirna.names)


ggplot(data = abundant.mirna.long, aes(x = miRNA, y = cpmValue )) + geom_boxplot() + ylab("cpm value") + ggtitle("20 most abundant miRNAs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
