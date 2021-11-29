library(edgeR)
library(xtable)
library(RColorBrewer)
source("load_count_data.R")


# Make dge object
dge <- DGEList(counts=count.df)

# Make groups
treatment = ifelse(metadata.df$probiotic ==1, "P", "nP")
outcome = ifelse(metadata.df$ad == 1, "AD", "nAD")

groups_p = as.factor(treatment)
groups = as.factor(paste0(treatment,"_", outcome))
groups_ad = as.factor(outcome)


# Assign groups to dge.g
dge.g = dge
dge.g$samples$group = groups


# convert to cpm and log cpm ( voom recomputes its own log-CPM values internally with a smaller prior count)
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
summary(lcpm)


# TABLE: not expressed miRNAs
table.zero.expression = table(rowSums(dge$counts==0)==60)



# filter out low expressed miRNA
keep.exprs <- filterByExpr(dge, group = groups)
dge.g <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge.g$counts)


# second way of filtering, using design matrix
design = model.matrix(~probiotic*ad, data = metadata.df)
head(design)
keep <- filterByExpr(dge, design = design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dim(dge$counts)



# check that design matrix gives same result as groups
treat = ifelse(design[,2] == 1, "P", "nP")
outc = ifelse(design[,3] == 1, "AD", "nAD")
groups_design = as.factor(paste0(treatment,"_", outcome))
groups == groups_design



#Optional? normalization. Not used here
# dge = calcNormFactors(dge, method = "TMM")


#MDS Plot
dge.g$samples$group = groups_p
lcpm = cpm(dge.g, log =  TRUE) 

col.group = groups_p       # or= groups     # or = groups_p
levels(col.group) =  brewer.pal(nlevels(col.group), "Set1")
col.group = as.character(col.group)
plotMDS(lcpm, labels=groups_p, col=col.group)
title(main="MDS plot")

#save plot
#svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\MDS-plot.svg")
#plotMDS(lcpm, labels=groups_p, col=col.group)
#title(main="MDS plot")
#legend("topleft", as.character(unique(dge.g$samples$group)),col = unique(col.group) , pch = 20)
#dev.off() 





# Define design matrix
design.matrix = model.matrix(~probiotic + ad , data = metadata.df)


# design matrix using the groups
design.matrix.g = model.matrix(~0+groups)
colnames(design.matrix.g) <- gsub("groups", "", colnames(design.matrix.g))


#Make contrast matrix, denne er feil nå som interaction er med. MÅ ENDRES
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
  ADvsnAD = 0.5*(P_AD + nP_AD) - 0.5*(P_nAD + nP_nAD),
  PnADvsAll = P_nAD - (P_AD + nP_nAD + nP_AD)/3,
  levels = colnames(design.matrix.g))

contr.matrix.g



# VOOM METHOD

#voom.weights = voom(dge, design.matrix, plot=TRUE)
voom.weights.g = voom(dge.g, design.matrix.g, plot = TRUE) 

#save plot
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\mean-variance-trend.svg")
v = voom(dge.g, design.matrix.g, plot=TRUE)
# Close the graphics device
dev.off() 

voom.fit = lmFit(voom.weights, design.matrix)
voom.fit = contrasts.fit(voom.fit, contrasts=contr.matrix)
eB.fit = eBayes(voom.fit)

voom.fit.g = lmFit(voom.weights.g, design.matrix.g)
voom.fit.g = contrasts.fit(voom.fit.g, contrasts = contr.matrix.g)
eB.fit.g = eBayes(voom.fit.g, trend = FALSE)

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


# Results by Voom
voom.PvsnP = topTable(fit = eB.fit, coef = 1)
voom.PvsnP

voom.PnADvsnPAD.g = topTable(fit = eB.fit.g, coef = 1)
voom.PnADvsnPAD.g
xtable(voom.PnADvsnPAD.g %>% slice_head(n = 10))

voom.PnADvsAll = topTable(fit = eB.fit.g, coef = 4)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))


# LIMMA-TREND

lpcm = cpm(dge.g, log = TRUE)
limma.fit = lmFit(lcpm, design.matrix.g)
limma.fit.contr = contrasts.fit(limma.fit, contrasts = contr.matrix.g)
limma.fit.contr = eBayes(limma.fit.contr, trend = TRUE)
limma.PnADvsnPAD = topTable(limma.fit.contr, coef = 1)
limma.PnADvsnPAD

xtable(limma.PnADvsnPAD %>% slice_head(n = 10))

limma.PnADvsAll = topTable(limma.fit.contr, coef = 4)


#Including more covariates

design.matrix.full = model.matrix(~0+groups + sex + sib + matatopy, data = metadata.df)
design.matrix.full

voom.full = voom(dge.g, design.matrix.full, plot = TRUE) 
voom.fit.full = lmFit(voom.full, design.matrix.g)
voom.fit.full = contrasts.fit(voom.fit.full, contrasts = contr.matrix.g)
eB.fit.full = eBayes(voom.fit.full, trend = FALSE)

plotSA(eB.fit.full, main="Final model: Mean-variance trend")
summary(decideTests(eB.fit.full))

topTable(eB.fit.full, coef = 1)



#SJEKK COVARIATENE OG DERES SIGNIFIKANS
