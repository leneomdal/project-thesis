source("explore_mirna_data.R")

library(edgeR)

dge = DGEList(counts = ten.days.sample.df)
View(dge$counts)

dge = calcNormFactors(dge)


logCPM = log(ten.days.sample.df)
View(logCPM)


#Remove this
#logCPM = DGEList(counts = logCPM)

# design matrix

design = model.matrix(~probiotic, data = ten.days.meta.data)
head(design)

# Limma trend
fit = lmFit(logCPM, design)
fit = eBayes(fit, trend = TRUE)
topTable(fit, coef = "probiotic")                                            #basically top of summary

#including AD as a covariate

design.include.AD = model.matrix(~probiotic + ad_ukwp2yr, data = ten.days.meta.data)
head(design.include.AD)

fit.AD = lmFit(logCPM, design.include.AD)
fit.AD = eBayes(fit.AD, trend = TRUE)
topTable(fit.AD, coef = "probiotic")


# include interaction term, noe rart her med design matrix? hvordan skal den se ut her?

design.interaction = model.matrix(~probiotic*ad_ukwp2yr, data = ten.days.meta.data)
head(design.interaction, n = 20)

fit.interaction = lmFit(logCPM, design = design.interaction)
fit.interaction = eBayes(fit.interaction, trend = TRUE)
topTable(fit.interaction, coef = "probiotic")
