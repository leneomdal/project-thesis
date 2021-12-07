rm(list=ls()) 
library(edgeR)
library(xtable)
library(RColorBrewer)
library(stats)
source("load_count_data.R")


# Make dge object
dge <- DGEList(counts=count.df)

# Make groups
treatment = ifelse(metadata.df$probiotic ==1, "P", "nP")
outcome = ifelse(metadata.df$ad == 1, "AD", "nAD")

groups_p = as.factor(treatment)
#groups = as.factor(paste0(treatment,"_", outcome))
#groups_ad = as.factor(outcome)


# Assign groups to dge.g
dge.g = dge
#dge.g$samples$group = groups


# convert to cpm and log cpm ( voom recomputes its own log-CPM values internally with a smaller prior count)
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
summary(lcpm)


# TABLE: not expressed miRNAs
table.zero.expression = table(rowSums(dge$counts==0)==60)

# Filtering out miRNAs with cpm value less than 500 in 10% or less of the samples

keep = rowSums(cpm(dge)>500)>=6
dge$counts = dge$counts[keep,]
dim(dge)

#keep.g = rowSums(cpm(dge.g)>500)>=6
#dge.g$counts = dge.g$counts[keep,]
#dim(dge.g)



# filter out low expressed miRNA
#keep.exprs <- filterByExpr(dge, group = groups)
#dge.g <- dge[keep.exprs,, keep.lib.sizes=FALSE]
#dim(dge.g$counts)


# second way of filtering, using design matrix
#design = model.matrix(~probiotic*ad, data = metadata.df)
#keep <- filterByExpr(dge, design = design)
#dge <- dge[keep,,keep.lib.sizes=FALSE]
#dim(dge$counts)



# check that design matrix gives same result as groups
#treat = ifelse(design[,2] == 1, "P", "nP")
#outc = ifelse(design[,3] == 1, "AD", "nAD")
#groups_design = as.factor(paste0(treatment,"_", outcome))
#groups == groups_design



#Optional? normalization. Not used here
# dge = calcNormFactors(dge, method = "TMM")


#MDS Plot
dge$samples$group = groups_p
lcpm = cpm(dge, log =  TRUE) 

col.group = groups_p       
levels(col.group) =  brewer.pal(nlevels(col.group), "Set1")
col.group = as.character(col.group)
plotMDS(lcpm, labels=groups_p, col=col.group)
title(main="MDS plot")

#save plot
#svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\MDS-plot.svg")
#plotMDS(lcpm, labels=groups_p, col=col.group)
#title(main="MDS plot")
#legend(x = -1.25, y = 0.2, as.character(unique(dge.g$samples$group)),col = unique(col.group) , pch = 20)
#dev.off() 





# Define design matrix
design.matrix = model.matrix(~probiotic + ad +matatopy + sex +sib, data = metadata.df)
design.matrix

# design matrix using the groups
#design.matrix.g = model.matrix(~0+groups)
#colnames(design.matrix.g) <- gsub("groups", "", colnames(design.matrix.g))

# Design matrix including interaction
design.intr = model.matrix(~probiotic*ad +matatopy + sex +sib, data = metadata.df)
design.intr


# Fit model to examine interaction term (limma-trend)

lpcm = cpm(dge, log = TRUE)
fit = lmFit(lcpm, design.intr)
fit = eBayes(fit, trend = TRUE)
interaction.t = topTable(fit, coef = 7)
interaction.t

xtable(interaction.t %>% slice_head(n = 10))




#Make contrast matrix, denne må endres om alle covariater er med, bare slett de siste
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

# contrast matrix including interaction, denne må også endres, KAN EGT SLETTE HELE
contr.matrix <- makeContrasts(
  PvsnP = probiotic,
  ADvsnAD = ad,
  nPADvsP = (intercept + ad) - (intercept + probiotic + 0.5* ad),
  PADvsnPnAD = (intercept + probiotic + ad) - (intercept),
  PnADvsnPAD = (probiotic + intercept) - (intercept + ad),
  PnADvsAll = (probiotic + intercept) - (intercept + intercept + ad + intercept + ad + probiotic)/3,
  levels = colnames(design.matrix))
contr.matrix

#Alternative contrast matrix using groups
#contr.matrix.g <- makeContrasts(
#  PnADvsnPAD = P_nAD - nP_AD,
#  PvsnP = 0.5*(P_AD+P_nAD)-0.5*(nP_AD+nP_nAD),
#  ADvsnAD = 0.5*(P_AD + nP_AD) - 0.5*(P_nAD + nP_nAD),
#  PnADvsAll = P_nAD - (P_AD + nP_nAD + nP_AD)/3,
#  levels = colnames(design.matrix.g))



# vanilla limma
lcpm = cpm(dge, log = TRUE)
limma.fit = lmFit(lcpm, design.matrix)
fit.vanilla = eBayes(limma.fit, trend = FALSE)
fit.vanilla$s2.prior

xtable(topTable(fit.vanilla, coef = 2, adjust.method = "BH") %>% slice_head(n = 10))



# LIMMA-TREND
limma.fit = lmFit(lcpm, design.matrix)
fit.eB = eBayes(limma.fit, trend = TRUE)
#table for p vs nP
topTable(fit.eB, coef = 2)
topTable(fit.eB, coef = 2, adjust.method = "bonferroni")
xtable(topTable(fit.eB, coef = 2) %>% slice_head(n = 10))



# HVORDAN FUNGERER LIMMA TREND?
means.prior.df = data.frame(means = limma.fit$Amean, s2.prior = fit.eB$s2.prior, s2.posterior = fit.eB$s2.post)
plot(means.prior.df$means, limma.fit$sigma )
plot(lowess(limma.fit$Amean, sqrt(limma.fit$sigma)))
points(limma.fit$Amean, sqrt(means.prior.df$s2.post))


#save plot VANILLA LIMMA
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\prior-variance-vanilla-limma.svg")
plotSA(fit.vanilla, main="Mean variance points and prior variance")
dev.off()

#save plot LIMMA TREND
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\final-MV-trend-limma-trend.svg")
plotSA(fit.eB, main="limma-trend: Mean-variance trend", ylab = "Sqrt( standard deviation)")
legend("topright", as.character(c("limma-trend curve")),col = unique(col.group) , pch = "-", cex = 1)
dev.off()



# mean-variance trend in limma trend
ggplot(data = means.prior.df, aes(x = means, y = (s2.prior)^(1/4))) + geom_point() + 
  ggtitle("Estimated prior against mean expression for each miRNA") +
  xlab("Average log2-expression levels") + ylab("s^(1/2) prior")
ggplot(data = means.prior.df, aes(x = means, y = s2.posterior)) + geom_point() + 
  ggtitle("Estimated posterior against mean expression for each miRNA") +
  xlab("Average log2-expression levels") + ylab("s^2 posterior")









# VOOM METHOD

voom.weights = voom(dge, design.matrix, plot=TRUE)
#voom.weights.g = voom(dge.g, design.matrix.g, plot = TRUE) 



#save plot
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\mean-variance-trend-by-voom.svg")
v = voom(dge, design.matrix, plot=TRUE)
# Close the graphics device
dev.off() 

voom.fit = lmFit(voom.weights, design.matrix)
#voom.fit = contrasts.fit(voom.fit, contrasts=contr.matrix)
eB.voom.fit = eBayes(voom.fit)

topTable(eB.voom.fit, coef = 2)
xtable(topTable(eB.voom.fit, coef = 2) %>% slice_head(n = 10))

#voom.fit.g = lmFit(voom.weights.g, design.matrix.g)
#voom.fit.g = contrasts.fit(voom.fit.g, contrasts = contr.matrix.g)
#eB.fit.g = eBayes(voom.fit.g, trend = FALSE)

#save plot
svg("C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\plots\\final-MV-trend-voom.svg")
plotSA(eB.fit, main="Final model voom: Mean-variance trend")
dev.off()



# table for quick comparison, adjusted p-value cutoff is set at 5% by default
summary(decideTests(eB.fit.g))



#SJEKK COVARIATENE OG DERES SIGNIFIKANS
