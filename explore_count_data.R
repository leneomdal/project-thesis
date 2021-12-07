source("load_count_data.R")


# BARPLOT OF EXPRESSION VALUES OF ONE miRNA


df  = data.frame(x =colnames(count.df), y = as.numeric(count.df[8,])) 

ggplot(data = df, aes(x = x, y = y)) + geom_bar(stat = "identity")

# HOW MANY IN ECAH GROUP
sum(metadata.df$probiotic ==0)
# 30 in each


#TEST FOR GENDER DIFFERENCES

sum(metadata.df$sex == 0)
# 27 males overall
sum(metadata.df$sex == 1 & metadata.df$probiotic == 1)
#16 males in P group

sum(metadata.df$sex == 0 & metadata.df$probiotic == 0)
# 11 males in nP group

#HOW MANY MOTHERS have history of atopic dermatitis

sum(metadata.df$probiotic == 0 & metadata.df$matatopy == 1)
# 9 in P group



# BARPLOT OF LIBRARY SIZES
library.sizes = apply(count.df, 2,sum)
barplot(library.sizes)

ggplot(data = data.frame(  lib.s = library.sizes, library = names(library.sizes)), 
       aes(x =library , y = lib.s)) + geom_bar(stat = "identity") + ylab("Library size") +
  xlab("Library") + ggtitle("Barplot of library sizes") + theme(axis.text.x=element_blank())
ggsave("barplot-library-sizes.svg", path = path_plots, height = 4.5, width = 6)




# BOXPLOT OF MOST ABUNDANT miRNAs

sum.reads = apply(count.df, 1, sum)

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
