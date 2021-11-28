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

dge.g = dge
dge.g$samples$group = groups


# convert to cpm and log cpm ( voom recomputes its own log-CPM values internally with a smaller prior count)
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
summary(lcpm)


# TABLE: not expressed miRNAs
table(rowSums(dge$counts==0)==60)



# filter out low expressed miRNA
keep.exprs <- filterByExpr(dge, group = groups)
dge.g <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge.g$counts)



# design matrix
design = model.matrix(~probiotic*ad, data = metadata.df)
head(design)


# second way of filtering????????
keep <- filterByExpr(dge, design = design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dim(dge$counts)


# check that design matrix gives same result as groups
treat = ifelse(design[,2] == 1, "P", "nP")
outc = ifelse(design[,3] == 1, "AD", "nAD")
groups_design = as.factor(paste0(treatment,"_", outcome))
groups == groups_design



# check which miRNA are different between the two filterings
diff.index = c()
for(i in 1:nrow(dge.g$counts)){
  if(rownames(dge.g$counts)[i] %in% rownames(dge$counts)){
    next
  }
  else{diff.index = c(diff.index, which(rownames(dge.g$counts) == rownames(dge.g$counts)[i]))}
}
diff.index
# ER DET RART AT ALLE DISSE BLIR TATT BORT?
dge.g$counts[diff.index,]


#Check if some a filtered from dge.g and not dge
dim(dge.g$counts[diff.index,])
dim(dge.g$counts)[1] - dim(dge$counts)

for(i in 1:nrow(dge$counts)){
  if(rownames(dge$counts)[i] %in% rownames(dge.g$counts)){
    next
  }
  else{print(rownames(dge$counts)[i])}
}

# eventuell filtrering?

#keep.modified <- filterByExpr(dge, design = design, min.count = 0, min.total.count = 60*3 )
#dge.filtered <- dge[keep.modified,,keep.lib.sizes=FALSE]
#dim(dge.filtered$counts)
#dim(dge)



#Optional? normalization
# dge = calcNormFactors(dge, method = "TMM")


# Plot 
library(RColorBrewer)


col.group = groups       # or= groups     # or = groups_p
levels(col.group) =  brewer.pal(nlevels(col.group), "Set1")
col.group = as.character(col.group)
plotMDS(lcpm, labels=groups, col=col.group)
title(main="MDS plot")



# Define design matrix

design.matrix = model.matrix(~probiotic + ad , data = metadata.df)



# design matrix using the groups
design.matrix.g = model.matrix(~0+groups)
colnames(design.matrix.g) <- gsub("groups", "", colnames(design.matrix.g))


#Make contrast matrix
colnames(design.matrix)[1] = "intercept"
contr.matrix <- makeContrasts(
  PvsnP = probiotic,
  ADvsnAD = ad,
  nPADvsP = (intercept + ad) - (intercept + probiotic + 0.5* ad),
  PADvsnPnAD = (intercept + probiotic + ad) - (intercept),
  PnADvsnPAD = (probiotic + intercept) - (intercept + ad),
  PnADvsAll = (probiotic + intercept) - (intercept + intercept + ad + intercept + ad + probiotic)/3,
  
  levels = colnames(design.matrix))
contr.matrix


#Alternative contrast matrix
contr.matrix.g <- makeContrasts(
  PnADvsnPAD = P_nAD - nP_AD,
  PvsnP = 0.5*(P_AD+P_nAD)-0.5*(nP_AD+nP_nAD),
  levels = colnames(design.matrix.g))

contr.matrix.g



# Voom method

voom.weights = voom(dge, design.matrix, plot=TRUE)
voom.weights.g = voom(dge.g, design.matrix.g, plot = TRUE) 

#save plot
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\mean-variance-trend.svg")
v = voom(dge, design.matrix, plot=TRUE)
# Close the graphics device
dev.off() 

voom.fit = lmFit(voom.weights, design.matrix)
voom.fit = contrasts.fit(voom.fit, contrasts=contr.matrix)
eB.fit = eBayes(voom.fit)

voom.fit.g = lmFit(voom.weights.g, design.matrix.g)
voom.fit.g = contrasts.fit(voom.fit.g, contrasts = contr.matrix.g)
eB.fit.g = eBayes(voom.fit.g)
#save plot
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\final-MV-trend.svg")
plotSA(eB.fit.g, main="Final model: Mean-variance trend")
dev.off()


# table for quick comparison, adjusted p-value cutoff is set at 5% by default
summary(decideTests(eB.fit.g))

#stricter cut off by treat
tfit <- treat(voom.fit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


#Examining individual miRNA

voom.PvsnP = topTable(fit = eB.fit, coef = 1)
voom.PvsnP

voom.PvsnP.g = topTable(fit = eB.fit.g, coef = 2)
voom.PvsnP.g

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))



