library(edgeR)
source("load_count_data.R")


# Make dge object
dge <- DGEList(counts=count.df)

View(metadata.df)

treatment = ifelse(metadata.df$probiotic ==1, "P", "nP")
outcome = ifelse(metadata.df$ad == 1, "AD", "nAD")

groups_p = as.factor(treatment)
groups = as.factor(paste0(treatment,"_", outcome))
groups_ad = as.factor(outcome)

dge$samples$group = groups

# convert to cpm and log cpm ( voom recomputes its own log-CPM values internally with a smaller prior count)
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
summary(lcpm)


# TABLE: not expressed miRNAs
table(rowSums(dge$counts==0)==60)



# filter out low expressed miRNA
keep.exprs <- filterByExpr(dge, group = groups)
dge1 <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge1$counts)


View(metadata.df)

groups

# design matrix
design = model.matrix(~probiotic + ad, data = metadata.df)
head(design)





# second way of filtering????????
keep <- filterByExpr(dge, design = design)
dge2 <- dge[keep,,keep.lib.sizes=FALSE]
dim(dge2$counts)


# check that design matrix gives same result as groups
tret = ifelse(design[,2] == 1, "P", "nP")
outc = ifelse(design[,3] == 1, "AD", "nAD")
groups_design = as.factor(paste0(treatment,"_", outcome))
groups == groups_design



# check which miRNA are different between the two filterings
diff.index = c()
for(i in 1:nrow(dge1$counts)){
  if(rownames(dge1$counts)[i] %in% rownames(dge2$counts)){
    next
  }
  else{diff.index = c(diff.index, which(rownames(dge1$counts) == rownames(dge1$counts)[i]))}
}
diff.index
# ER DET RART AT ALLE DISSE BLIR TATT BORT?
dge1$counts[diff.index,]


#Check if some a filtered from dge1 and not dge2
dim(dge1$counts[diff.index,])
dim(dge1$counts)[1] - dim(dge2$counts)

for(i in 1:nrow(dge2$counts)){
  if(rownames(dge2$counts)[i] %in% rownames(dge1$counts)){
    next
  }
  else{print(rownames(dge2$counts)[i])}
}

# eventuell filtrering?
View(dge$counts)
summary(dge)

keep.modified <- filterByExpr(dge, design = design, min.count = 0, min.total.count = 60*3 )
dge.filtered <- dge[keep.modified,,keep.lib.sizes=FALSE]
dim(dge.filtered$counts)
dim(dge2)

# CHOOSE TO CONTINIUE WITH dge1

dge = dge1
dge


#Optional? normalization
# dge = calcNormFactors(dge, method = "TMM")


# Plot 
library(RColorBrewer)


lcpm = cpm(dge, log=TRUE)
col.group = groups_ad       # or= groups     # or = groups_p
levels(col.group) =  brewer.pal(nlevels(col.group), "Set1")
col.group = as.character(col.group)
plotMDS(lcpm, labels=groups, col=col.group)
title(main="MDS plot")



# Define design matrix

design.matrix = model.matrix(~probiotic + ad , data = metadata.df)

View(design.matrix)

# design matrix using
#design.matrix = model.matrix(~0+groups)


#Make contrast matrix
colnames(design.matrix) <- gsub("group", "", colnames(design.matrix))

contr.matrix <- makeContrasts(
  PvsnP = probiotic,
  ADvsnAD = ad,
  nPADvsP = (intercept + ad) - (intercept + probiotic + 0.5* ad),
  PADvsnPnAD = (intercept + probiotic + ad) - (intercept),
  PnADvsnPAD = (probiotic + intercept) - (intercept + ad),
  PnADvsAll = (probiotic + intercept) - (intercept + intercept + ad + intercept + ad + probiotic)/3,
  
  levels = colnames(design.matrix))

contr.matrix

colnames(design.matrix)[1] = "intercept"

# Voom method

voom.weights <- voom(dge, design.matrix, plot=TRUE)

#save plot
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\mean-variance-trend.svg")
v = voom(dge, design.matrix, plot=TRUE)
# Close the graphics device
dev.off() 

voom.fit <- lmFit(voom.weights, design.matrix)
voom.fit <- contrasts.fit(voom.fit, contrasts=contr.matrix)
eB.fit <- eBayes(voom.fit)
#save plot
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\final-MV-trend.svg")
plotSA(eB.fit, main="Final model: Mean-variance trend")
dev.off()


# table for quick comparison, adjusted p-value cutoff is set at 5% by default
summary(decideTests(eB.fit))

#stricter cut off by treat
tfit <- treat(voom.fit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


#Examining individual miRNA

voom.PvsnP = topTable(fit = eB.fit, coef = 5)
voom.PvsnP
